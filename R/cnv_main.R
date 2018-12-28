# License Part ------------------------------------------------------------

# MIT License
#
# Copyright (c) [2018] [Geoffrey Macintyre]
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#     The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#------------------------------------------------------------------
#' @title  Read copy number as a list of data.frame from data.frame or files
#' @description This function is used to read copy number profile for preparing CNV signature
#' analysis.
#' @param input a `data.frame` or a file or a directory contains copy number profile.
#' @param is_dir logical. Is `input` a directory?
#' @param pattern an optional regular expression used to select part of files if input is a directory, more detail please see `list.files` function.
#' @param ignore_case logical. Should pattern-matching be case-insensitive?
#' @param sep the field separator character for input file(s) if `input` is not `data.frame`.
#' @param cols four characters used to specify chromosome, start position,
#'  end position and copy number value in input, respectively.
#'  Default use names from ABSOLUTE calling result.
#' @param have_sampleCol logical. Does input have sample column?
#' This parameter must be `TRUE` and `sample_col` also
#' properly be assigned when `input` is a file or a \code{data.frame} (i.e. not a directory).
#' @param sample_col a character used to specify the sample column name.
#' @author  Shixiang Wang <w_shixiang@163.com>
#' @return a `list` contains absolute copy-number profile for multiple samples.
#' @importFrom utils read.csv
#' @export
#' @family CNV analysis functions
#' @seealso [cnv_derivefeatures()] for deriving CNV features, [cnv_getLengthFraction()] for calculating
#' CNV length fraction (normalized to arm), [cnv_plotDistributionProfile()] for plotting profile of
#' CNV distribution.
cnv_readprofile = function(input,
                           is_dir = FALSE,
                           pattern = NULL,
                           ignore_case = FALSE,
                           sep = "\t",
                           cols = c("Chromosome", "Start.bp", "End.bp", "modal_cn"),
                           have_sampleCol = TRUE,
                           sample_col = "sample") {
    stopifnot(
        is.logical(is_dir),
        is.logical(have_sampleCol),
        is.character(sample_col),
        length(sample_col) == 1
    )
    if (is_dir) {
        message("Treat input as a directory...")
        if (length(input) != 1) {
            stop("Only can take one directory as input!")
        }
        # get files and exclude directories
        all.files <- list.files(
            path = input,
            pattern = pattern,
            all.files = FALSE,
            recursive = FALSE,
            ignore.case = ignore_case
        )
        files = all.files[!file.info(file.path(input, all.files))$isdir]
        if (length(files) == 0)
            stop("No files exist, please check!")
        files_path = file.path(input, files)
        files_list = list()
        for (i in seq_along(files_path)) {
            temp = read.csv(
                file = files_path[i],
                sep = sep,
                comment.char = "#",
                stringsAsFactors = FALSE
            )
            if (!all(cols %in% colnames(temp)))
                stop("not all cols are in file, please check.")
            if (have_sampleCol) {
                tempName = unique(temp[, sample_col])
                if (length(tempName) > 1) {
                    stop("When input is a directory, a file contains only one sample.")
                }
                temp = temp[, cols]
                colnames(temp) = c("chromosome", "start", "end", "segVal")
                if (nrow(temp) <= 22) {
                    warning("Sample ",
                            tempName,
                            " is discarded because of few segments (<=22)")
                } else {
                    files_list[[tempName]] = temp
                }

            } else {
                message("Select file names as sample names.")
                temp = temp[, cols]
                colnames(temp) = c("chromosome", "start", "end", "segVal")
                if (nrow(temp) <= 22) {
                    warning("File ",
                            files[i],
                            " is discarded because of few segments (<=22)")
                } else {
                    files_list[[files[i]]] = temp
                }

            }
        }
        return(files_list)
    } else if (all(is.character(input))) {
        message("Treat input as a file...")
        if (length(input) > 1) {
            stop("Muliple files are not a valid input, please use directory as input.")
        }

        if (!file.exists(input))
            stop("input file not exists")
        if (!have_sampleCol)
            stop("When input is a file, sample column must set.")
        input = read.csv(
            file = input,
            sep = sep,
            comment.char = "#",
            stringsAsFactors = FALSE
        )
    }

    if (!sample_col %in% colnames(input))
        stop("sample column user set not exists in input file.")
    if (!all(cols %in% colnames(input)))
        stop("not all cols are in file, please check.")
    samples = unique(input[, sample_col])

    res_list = list()
    for (i in seq_along(samples)) {
        tempDF = input[input[, sample_col] == samples[i],]
        tempDF = tempDF[, cols]
        colnames(tempDF) = c("chromosome", "start", "end", "segVal")

        if (nrow(tempDF) <= 22) {
            warning("Sample ",
                    samples[i],
                    " is discarded because of few segments (<=22)")
        } else {
            res_list[[samples[i]]] = tempDF
        }
    }

    return(res_list)
}

#-------------------------------------------

#' @title  Derive copy number feature distributions
#' @description This function summarise each copy-number profile using a number of different
#' feature distributions: sigment size, breakpoint number (per ten megabase), change-point copy-number,
#' segment copy-number, breakpoint number (per chromosome arm), length of segments with oscilating
#' copy-number.
#' @param CN_data a `QDNAseqCopyNumbers` object or a `list` contains multiple `data.frame`s (recommended),
#' each `data.frame` stores copy-number profile for one sample with 'chromosome', 'start', 'end' and
#' 'segVal' these four necessary columns. Of note, 'segVal' column shoule be absolute copy number values.
#' @param cores number of compute cores to run this task.
#' You can use [parallel::detectCores()] function to check how
#' many cores you can use. If you are using [cnv_pipe()] feature,
#' please do not use maximal number of
#' cores in your computer, it may cause some unexpected problems.
#' @param genome_build genome build version, must be one of 'hg19' or 'hg38'.
#' @author Geoffrey Macintyre, Shixiang Wang
#' @return a `list` contains six copy number feature distributions.
#' @import foreach doParallel
#' @export
#' @family CNV analysis functions
#' @seealso [cnv_plotFeatureDistribution()] for plotting feature distributions.
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = cnv_derivefeatures(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#'}
cnv_derivefeatures = function(CN_data,
                              cores = 1,
                              genome_build = c("hg19", "hg38")) {
    genome_build = match.arg(genome_build)
    # get chromosome lengths and centromere locations
    if (genome_build == "hg19") {
        data("chromsize.hg19",
             package = "VSHunter",
             envir = environment())
        data("centromeres.hg19",
             package = "VSHunter",
             envir = environment())
        chrlen = chromsize.hg19
        centromeres = centromeres.hg19
    } else {
        data("chromsize.hg38",
             package = "VSHunter",
             envir = environment())
        data("centromeres.hg38",
             package = "VSHunter",
             envir = environment())
        chrlen = chromsize.hg38
        centromeres = centromeres.hg38
    }

    # only keep 1:22 and x, y
    chrlen = chrlen[chrlen$chrom %in% centromeres$chrom,]
    if (cores > 1) {
        #require(foreach)
        requireNamespace("foreach", quietly = TRUE)
        #doMC::registerDoMC(cores)
        doParallel::registerDoParallel(cores = cores)

        temp_list = foreach::foreach(i = 1:6) %dopar% {
            if (i == 1) {
                list(segsize = getSegsize(CN_data))
            } else if (i == 2) {
                list(bp10MB = getBPnum(CN_data, chrlen))
            } else if (i == 3) {
                list(osCN = getOscilation(CN_data))
            } else if (i == 4) {
                list(bpchrarm = getCentromereDistCounts(CN_data, centromeres, chrlen))
            } else if (i == 5) {
                list(changepoint = getChangepointCN(CN_data))
            } else {
                list(copynumber = getCN(CN_data))
            }

        }
        unlist(temp_list, recursive = FALSE)
    } else {
        segsize <- getSegsize(CN_data)
        bp10MB <- getBPnum(CN_data, chrlen)
        osCN <- getOscilation(CN_data)
        bpchrarm <-
            getCentromereDistCounts(CN_data, centromeres, chrlen)
        changepoint <- getChangepointCN(CN_data)
        copynumber <- getCN(CN_data)

        list(
            segsize = segsize,
            bp10MB = bp10MB,
            osCN = osCN,
            bpchrarm = bpchrarm,
            changepoint = changepoint,
            copynumber = copynumber
        )
    }

}


#------------------------------------------------
#
#' @title Fit optimal number of mixture model components
#' @description Apply mixture modelling to breakdown each feature distribution into mixtures
#' of Gaussian or mixtures of Poison distributions using the **flexmix** package.
#'
#' @param CN_features a `list` generate from [cnv_derivefeatures()] function.
#' @param seed seed number.
#' @param min_comp minimal number of components to fit, default is 2.
#' @param max_comp maximal number of components to fit, default is 10.
#' @param min_prior minimal prior value, default is 0.001.
#' Details about custom setting please refer to **flexmix** package.
#' @param model_selection model selection strategy, default is 'BIC'.
#' Details about custom setting please refer to **flexmix** package.
#' @param nrep number of run times for each value of component,
#' keep only the solution with maximum likelihood.
#' @param niter maximal number of iteration to achive converge.
#' @inheritParams cnv_derivefeatures
#' @param featsToFit integer vector used for task assignment in parallel computation.
#' **Do not change it!**
#' @author Geoffrey Macintyre, Shixiang Wang
#' @return a `list` contain `flexmix` object of copy-number features.
#' @import flexmix
#' @export
#' @family CNV analysis functions
#' @seealso [cnv_plotMixComponents()] for plotting mixture component models.
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = cnv_derivefeatures(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' ## fit mixture model  (this will take some time)
#' tcga_components = cnv_fitMixModels(CN_features = tcga_features, cores = 1)
#' }
cnv_fitMixModels = function(CN_features,
                            seed = 123456,
                            min_comp = 2,
                            max_comp = 10,
                            min_prior = 0.001,
                            model_selection = "BIC",
                            nrep = 1,
                            niter = 1000,
                            cores = 1,
                            featsToFit = seq(1, 6)) {
    if (cores > 1) {
        #require(foreach)
        requireNamespace("foreach", quietly = TRUE)
        #doMC::registerDoMC(cores)
        doParallel::registerDoParallel(cores = cores)

        temp_list = foreach(i = 1:6) %dopar% {
            if (i == 1 & i %in% featsToFit) {
                message("Fit feature: Segment size")
                dat <- as.numeric(CN_features[["segsize"]][, 2])
                list(
                    segsize = fitComponent(
                        dat,
                        seed = seed,
                        model_selection = model_selection,
                        min_prior = min_prior,
                        niter = niter,
                        nrep = nrep,
                        min_comp = min_comp,
                        max_comp = max_comp
                    )
                )

            } else if (i == 2 & i %in% featsToFit) {
                message("Fit feature: Breakpoint count per 10 Mb")
                dat <- as.numeric(CN_features[["bp10MB"]][, 2])
                list(
                    bp10MB = fitComponent(
                        dat,
                        dist = "pois",
                        seed = seed,
                        model_selection = model_selection,
                        min_prior = min_prior,
                        niter = niter,
                        nrep = nrep,
                        min_comp = min_comp,
                        max_comp = max_comp
                    )
                )

            } else if (i == 3 & i %in% featsToFit) {
                message("Fit feature: Length of oscillating copy-number chain")
                dat <- as.numeric(CN_features[["osCN"]][, 2])
                list(
                    osCN = fitComponent(
                        dat,
                        dist = "pois",
                        seed = seed,
                        model_selection = model_selection,
                        min_prior = min_prior,
                        niter = niter,
                        nrep = nrep,
                        min_comp = min_comp,
                        max_comp = max_comp
                    )
                )

            } else if (i == 4 & i %in% featsToFit) {
                message("Fit feature: Breakpoint count per arm")
                dat <- as.numeric(CN_features[["bpchrarm"]][, 2])
                list(
                    bpchrarm = fitComponent(
                        dat,
                        dist = "pois",
                        seed = seed,
                        model_selection = model_selection,
                        min_prior = min_prior,
                        niter = niter,
                        nrep = nrep,
                        min_comp = min_comp,
                        max_comp = max_comp
                    )
                )

            } else if (i == 5 & i %in% featsToFit) {
                message("Fit feature: Copy number change")
                dat <- as.numeric(CN_features[["changepoint"]][, 2])
                list(
                    changepoint = fitComponent(
                        dat,
                        seed = seed,
                        model_selection = model_selection,
                        min_prior = min_prior,
                        niter = niter,
                        nrep = nrep,
                        min_comp = min_comp,
                        max_comp = max_comp
                    )
                )

            } else if (i == 6 & i %in% featsToFit) {
                message("Fit feature: Absolute copy number")
                dat <- as.numeric(CN_features[["copynumber"]][, 2])
                list(
                    copynumber = fitComponent(
                        dat,
                        seed = seed,
                        model_selection = model_selection,
                        nrep = nrep,
                        min_comp = min_comp,
                        max_comp = max_comp,
                        min_prior = 0.005,
                        niter = 2000
                    )
                )

            }

        }
        unlist(temp_list, recursive = FALSE)
    } else {
        dat <- as.numeric(CN_features[["segsize"]][, 2])
        message("Fit feature: Segment size")
        segsize_mm <-
            fitComponent(
                dat,
                seed = seed,
                model_selection = model_selection,
                min_prior = min_prior,
                niter = niter,
                nrep = nrep,
                min_comp = min_comp,
                max_comp = max_comp
            )

        dat <- as.numeric(CN_features[["bp10MB"]][, 2])
        message("Fit feature: Breakpoint count per 10 Mb")
        bp10MB_mm <-
            fitComponent(
                dat,
                dist = "pois",
                seed = seed,
                model_selection = model_selection,
                min_prior = min_prior,
                niter = niter,
                nrep = nrep,
                min_comp = min_comp,
                max_comp = max_comp
            )

        dat <- as.numeric(CN_features[["osCN"]][, 2])
        message("Fit feature: Length of oscillating copy-number chain")
        osCN_mm <-
            fitComponent(
                dat,
                dist = "pois",
                seed = seed,
                model_selection = model_selection,
                min_prior = min_prior,
                niter = niter,
                nrep = nrep,
                min_comp = min_comp,
                max_comp = max_comp
            )

        dat <- as.numeric(CN_features[["bpchrarm"]][, 2])
        message("Fit feature: Breakpoint count per arm")
        bpchrarm_mm <-
            fitComponent(
                dat,
                dist = "pois",
                seed = seed,
                model_selection = model_selection,
                min_prior = min_prior,
                niter = niter,
                nrep = nrep,
                min_comp = min_comp,
                max_comp = max_comp
            )

        dat <- as.numeric(CN_features[["changepoint"]][, 2])
        message("Fit feature: Copy number change")
        changepoint_mm <-
            fitComponent(
                dat,
                seed = seed,
                model_selection = model_selection,
                min_prior = min_prior,
                niter = niter,
                nrep = nrep,
                min_comp = min_comp,
                max_comp = max_comp
            )

        dat <- as.numeric(CN_features[["copynumber"]][, 2])
        message("Fit feature: Absolute copy number")
        copynumber_mm <-
            fitComponent(
                dat,
                seed = seed,
                model_selection = model_selection,
                nrep = nrep,
                min_comp = min_comp,
                max_comp = max_comp,
                min_prior = 0.005,
                niter = 2000
            )

        list(
            segsize = segsize_mm,
            bp10MB = bp10MB_mm,
            osCN = osCN_mm,
            bpchrarm = bpchrarm_mm,
            changepoint = changepoint_mm,
            copynumber = copynumber_mm
        )
    }
}
#------------------------------------
#' @title Generate a sample-by-component matrix
#' @description This function generate a sample-by-component matrix representing the sum of
#' posterior probabilities of each copy-number event being assigned to each component.
#' @param CN_features a `list` contains six copy number feature distributions,
#' obtain this from [cnv_derivefeatures()] function.
#' @param all_components a `list` contain `flexmix` object of copy-number features, obtain this
#' from [cnv_fitMixModels] function or use pre-compiled components data which come from CNV signature paper
#' https://www.nature.com/articles/s41588-018-0179-8 (set this parameter as `NULL`).
#' @inheritParams cnv_derivefeatures
#' @param rowIter step size of iteration for rows of ech CNV feature \code{data.frame}.
#' @author Geoffrey Macintyre, Shixiang Wang
#' @import doParallel
#' @return a numeric sample-by-component `matrix`
#' @export
#' @family CNV analysis functions
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = cnv_derivefeatures(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' ## fit mixture model  (this will take some time)
#' tcga_components = cnv_fitMixModels(CN_features = tcga_features, cores = 1)
#' ## generate a sample-by-component matrix
#' tcga_sample_component_matrix = cnv_generateSbCMatrix(tcga_features, tcga_components, cores = 1)
#' }
cnv_generateSbCMatrix = function(CN_features,
                                 all_components = NULL,
                                 cores = 1,
                                 rowIter = 1000)
{
    if (is.null(all_components))
    {
        # all_components <-
        #     readRDS(paste(this_path, "data/component_parameters.rds", sep = "/"))
        message(
            "About reference components\n   more detail please see https://github.com/ShixiangWang/absoluteCNVdata"
        )
        if (!file.exists("Nat_Gen_component_parameters.rds")) {
            message(
                "Nat_Gen_component_parameters.rds doesnot exist, will download reference components."
            )
            download.file(url = "https://github.com/ShixiangWang/absoluteCNVdata/raw/master/component_parameters.rds",
                          destfile = "Nat_Gen_component_parameters.rds")
        }
        all_components <-
            readRDS("Nat_Gen_component_parameters.rds")
    }



    full_mat <- cbind(
        calculateSumOfPosteriors(CN_features[["segsize"]], all_components[["segsize"]], "segsize", cores = cores),
        calculateSumOfPosteriors(CN_features[["bp10MB"]], all_components[["bp10MB"]], "bp10MB", cores = cores),
        calculateSumOfPosteriors(CN_features[["osCN"]], all_components[["osCN"]], "osCN", cores = cores),
        calculateSumOfPosteriors(CN_features[["changepoint"]], all_components[["changepoint"]], "changepoint", cores = cores),
        calculateSumOfPosteriors(CN_features[["copynumber"]], all_components[["copynumber"]], "copynumber", cores = cores),
        calculateSumOfPosteriors(CN_features[["bpchrarm"]], all_components[["bpchrarm"]], "bpchrarm", cores = cores)
    )

    # if (cores > 1) {
    #     require(foreach)
    #     requireNamespace("foreach", quietly = TRUE)
    #     feats = c("segsize",
    #               "bp10MB",
    #               "osCN",
    #               "changepoint",
    #               "copynumber",
    #               "bpchrarm")
    #     doMC::registerDoMC(cores)
    #
    #     full_mat = foreach(feat = feats, .combine = cbind) %dopar% {
    #         calculateSumOfPosteriors(
    #             CN_features[[feat]],
    #             all_components[[feat]],
    #             feat,
    #             rowIter = rowIter,
    #             cores = cores
    #         )
    #     }
    # } else {
    #     full_mat <- cbind(
    #         calculateSumOfPosteriors(CN_features[["segsize"]], all_components[["segsize"]], "segsize"),
    #         calculateSumOfPosteriors(CN_features[["bp10MB"]], all_components[["bp10MB"]], "bp10MB"),
    #         calculateSumOfPosteriors(CN_features[["osCN"]], all_components[["osCN"]], "osCN"),
    #         calculateSumOfPosteriors(CN_features[["changepoint"]], all_components[["changepoint"]], "changepoint"),
    #         calculateSumOfPosteriors(CN_features[["copynumber"]], all_components[["copynumber"]], "copynumber"),
    #         calculateSumOfPosteriors(CN_features[["bpchrarm"]], all_components[["bpchrarm"]], "bpchrarm")
    #     )
    # }

    rownames(full_mat) <- unique(CN_features[["segsize"]][, 1])
    full_mat[is.na(full_mat)] <- 0
    full_mat
}


#------------------------------------
#' @title Choose optimal number of signatures
#' @description This function use **NMF** package to evaluate the optimal number of signatures.
#' The most common approach is to choose the smallest rank for which cophenetic correlation coefficient
#' starts decreasing (Used by this function). Another approach is to choose the rank for which the plot
#' of the residual sum of squares (RSS) between the input matrix and its estimate shows an inflection point.
#' @param sample_by_component a sample-by-component `matrix`, generate from [cnv_generateSbCMatrix] function.
#' @param nTry the maximal tried number of signatures, default is 12.
#' Of note, this value should far less than number of features or samples.
#' @param nrun the number of run to perform for each value in range of 2 to `nTry`, default is 10.
#' According to **NMF** package documentation, `nrun` set to 30~50 is enough to achieve robust result.
#' @inheritParams cnv_derivefeatures
#' @param seed seed number.
#' @param plot logical. If `TRUE`, plot rank survey.
#' @param consensusmap_name a character, basename of consensus map output path.
#' @param testRandom Should generate random data from input to test measurements. Default is `TRUE`.
#' @param nmfalg specification of the NMF algorithm.
#' @author Geoffrey Macintyre, Shixiang Wang
#' @importFrom NMF nmfEstimateRank consensusmap summary
#' @importFrom grDevices pdf dev.off
#' @return a `list` contains information of NMF run and rank survey.
#' @export
#' @family CNV analysis functions
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = cnv_derivefeatures(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' ## fit mixture model  (this will take some time)
#' tcga_components = cnv_fitMixModels(CN_features = tcga_features, cores = 1)
#' ## generate a sample-by-component matrix
#' tcga_sample_component_matrix = cnv_generateSbCMatrix(tcga_features, tcga_components, cores = 1)
#' ## optimal rank survey
#'  tcga_sig_choose = cnv_chooseSigNumber(tcga_sample_component_matrix,
#'  nrun = 10, cores = 1, plot = FALSE)
#' }
#'
cnv_chooseSigNumber <-
    function(sample_by_component,
             nTry = 12,
             nrun = 10,
             cores = 1,
             seed = 123456,
             plot = TRUE,
             consensusmap_name = "nmf_consensus",
             testRandom = TRUE,
             nmfalg = "brunet")
    {
        message('Estimating best rank..')
        #nmfalg <- "brunet"

        #suppressMessages(library(NMF))
        estim.r <-
            NMF::nmfEstimateRank(
                t(sample_by_component),
                seq(2, nTry),
                seed = seed,
                nrun = nrun,
                verbose = TRUE,
                method = nmfalg,
                .opt = list(shared.memory = FALSE, paste0("p", cores))
            )

        pdf(
            paste0(consensusmap_name, ".pdf"),
            bg = 'white',
            pointsize = 9,
            width = 12,
            height = 12,
            paper = "special"
        )
        NMF::consensusmap(estim.r)
        dev.off()
        message('created ', paste0(consensusmap_name, ".pdf"))

        #--- copy from maftools and modified ---#
        nmf.sum = NMF::summary(estim.r) # Get summary of estimates
        print(nmf.sum)
        nmf.sum$diff = c(0, diff(nmf.sum$cophenetic))
        bestFit = nmf.sum$rank[which(nmf.sum$diff < 0)][1]
        #bestFit = nmf.sum[diff < 0, rank][1] #First point where cophenetic correlation coefficient starts decreasing

        # https://blog.csdn.net/YJJ18636810884/article/details/83214566
        # https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-367
        message(
            paste(
                'Using ',
                bestFit,
                ' as a best-fit rank based on decreasing cophenetic correlation coefficient.',
                sep = ''
            )
        )
        n = as.numeric(bestFit)

        if (testRandom) {
            message("Generating random matrix and run NMF...")
            V.random <- NMF::randomize(t(sample_by_component))
            estim.r.random <-
                NMF::nmfEstimateRank(
                    V.random,
                    seq(2, nTry),
                    seed = seed,
                    nrun = nrun,
                    verbose = TRUE,
                    method = nmfalg,
                    .opt = list(shared.memory = FALSE, paste0("p", cores))
                )
        }

        if (plot) {
            message('Creating nmf rank survey plot...')

            if (testRandom) {
                p <- NMF::plot(
                    estim.r,
                    estim.r.random,
                    what = c("cophenetic",
                             "dispersion",
                             "sparseness",
                             #"silhouette",
                             #"residuals",
                             "rss"),
                    xname = "Observed",
                    yname = "Randomised",
                    main = "NMF Rank Survey"
                )
            } else {
                p <- NMF::plot(
                    estim.r,
                    what = c("cophenetic",
                             "dispersion",
                             "sparseness",
                             # "silhouette",
                             # "residuals",
                             "rss"),
                    main = "NMF Rank Survey"
                )
            }

            print(p)

        }

        if (!plot)
            p = NULL
        if (!testRandom)
            estim.r.random = NULL

        return(
            list(
                nmfEstimate = estim.r,
                nmfEstimate.random = estim.r.random,
                bestRank = n,
                survey = nmf.sum,
                survey_plot = p,
                seed = seed
            )
        )
    }

#--------------------------
# extract signatures
#' @title Extract signature based on specified rank value
#' @inheritParams cnv_chooseSigNumber
#' @param nsig specification of the factorization rank.
#' @author Geoffrey Macintyre, Shixiang Wang
#' @importFrom NMF nmf
#' @return a object of \code{NMF} run.
#' @export
#' @family CNV analysis functions
#' @seealso [cnv_plotSignatures()] for plot signatures and their contributions.
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = cnv_derivefeatures(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' ## fit mixture model  (this will take some time)
#' tcga_components = cnv_fitMixModels(CN_features = tcga_features, cores = 1)
#' ## generate a sample-by-component matrix
#' tcga_sample_component_matrix = cnv_generateSbCMatrix(tcga_features, tcga_components, cores = 1)
#' ## optimal rank survey
#'  tcga_sig_choose = cnv_chooseSigNumber(tcga_sample_component_matrix, nrun = 10,
#'  cores = 1, plot = FALSE)
#'  tcga_signatures = cnv_extractSignatures(tcga_sample_component_matrix, nsig = 3, cores = 1)
#' }
cnv_extractSignatures <-
    function(sample_by_component,
             nsig,
             seed = 123456,
             nmfalg = "brunet",
             cores = 1)
    {
        message("Running NMF based on specified rank...")
        #suppressMessages(library(NMF))
        NMF::nmf(
            t(sample_by_component),
            nsig,
            seed = seed,
            nrun = 1000,
            method = nmfalg,
            .opt = paste0("p", cores)
        )
    }

#---------------------------------------------------------------------------
#' @title Quantify exposure for samples using Linear Combination Decomposition (LCD)
#'
#' @inheritParams cnv_chooseSigNumber
#' @param component_by_signature a componet by signature matrix,
#' default is `NULL`, it will use pre-compiled data from CNV signature paper
#' https://www.nature.com/articles/s41588-018-0179-8.
#' @author Geoffrey Macintyre, Shixiang Wang
#' @return a `list` contains absolute/relative exposure.
#' @export
#' @inherit cnv_extractSignatures seealso
#' @family CNV analysis functions
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = cnv_derivefeatures(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' ## fit mixture model  (this will take some time)
#' tcga_components = cnv_fitMixModels(CN_features = tcga_features, cores = 1)
#' ## generate a sample-by-component matrix
#' tcga_sample_component_matrix = cnv_generateSbCMatrix(tcga_features, tcga_components, cores = 1)
#' ## optimal rank survey
#'  tcga_sig_choose = cnv_chooseSigNumber(tcga_sample_component_matrix,
#'  nrun = 10, cores = 1, plot = FALSE)
#'  tcga_signatures = cnv_extractSignatures(tcga_sample_component_matrix, nsig = 3, cores = 1)
#'  w = NMF::basis(tcga_signatures) # signature matrix
#'  tcga_exposure = cnv_quantifySigExposure(sample_by_component =
#'  tcga_sample_component_matrix, component_by_signature = w)
#' }
cnv_quantifySigExposure <-
    function(sample_by_component,
             component_by_signature = NULL)
    {
        if (!requireNamespace("YAPSA", quietly = TRUE)) {
            stop("Package \"YAPSA\" needed for this function to work. Please install it.",
                 call. = FALSE)
        }
        message("Quantifying exposure of signatures by LCD analysis...")
        if (is.null(component_by_signature))
        {
            message("Using reference component_by_signature data from Nat.Gen paper.")
            component_by_signature <-
                readRDS(system.file("extdata", "feat_sig_mat.rds", package = "VSHunter"))
        }
        signature_by_sample <- YAPSA::LCD(
            t(sample_by_component),
            YAPSA::normalize_df_per_dim(component_by_signature, 2)
        )
        absolute_exposure = as.matrix(signature_by_sample)
        relative_exposure = normaliseMatrix(signature_by_sample)
        return(
            list(
                absolute_exposure = absolute_exposure,
                relative_exposure = relative_exposure
            )
        )
    }

#-------------------------------------------------------------------------------------
#' Auto-capture signature and coresponding exposure
#' @description  this is a wrapper of \code{cnv_chooseSigNumber}, \code{cnv_extractSignatures}
#' and \code{cnv_quantifySigExposure} these three functions.
#'
#' @inheritParams cnv_chooseSigNumber
#' @author Geoffrey Macintyre, Shixiang Wang
#' @return a `list` contains results of NMF best rank survey, run, signature matrix, exposure list etc..
#' @importFrom NMF basis coef
#' @export
#' @inherit cnv_extractSignatures seealso
#' @family CNV analysis functions
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = cnv_derivefeatures(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' ## fit mixture model  (this will take some time)
#' tcga_components = cnv_fitMixModels(CN_features = tcga_features, cores = 1)
#' ## generate a sample-by-component matrix
#' tcga_sample_component_matrix = cnv_generateSbCMatrix(tcga_features, tcga_components, cores = 1)
#' ## optimal rank survey13
#' tcga_results = cnv_autoCaptureSignatures(tcga_sample_component_matrix, nrun=10, cores = 1)
#' }
cnv_autoCaptureSignatures = function(sample_by_component,
                                     nTry = 12,
                                     nrun = 10,
                                     cores = 1,
                                     seed = 123456,
                                     plot = TRUE,
                                     consensusmap_name = "nmf_consensus",
                                     testRandom = TRUE,
                                     nmfalg = "brunet") {
    choose_res = cnv_chooseSigNumber(
        sample_by_component,
        nTry,
        nrun,
        cores,
        seed,
        plot = plot,
        consensusmap_name = consensusmap_name,
        testRandom = testRandom,
        nmfalg = nmfalg
    )
    NMF_res = cnv_extractSignatures(sample_by_component,
                                    nsig = choose_res$bestRank,
                                    cores = cores, nmfalg = nmfalg)
    #-- Signatures
    w = NMF::basis(NMF_res)
    w = apply(w, 2, function(x)
        x / sum(x)) # Scale signatures (basis)
    colnames(w) = paste("Signature", 1:ncol(w), sep = "_")

    #-- Contributions
    h = NMF::coef(NMF_res)
    h = apply(h, 2, function(x)
        x / sum(x)) # Scale contributions (coefs)

    rownames(h) = paste("Signature", 1:ncol(w), sep = "_")
    message("Done.")

    return(
        list(
            NMF = NMF_res,
            signatures = w,
            contributions = h,
            nmfEstimate = choose_res$nmfEstimate,
            nmfEstimate.random = choose_res$nmfEstimate.random,
            bestRank = choose_res$bestRank,
            survey = choose_res$survey,
            survey_plot = choose_res$survey_plot,
            seed = choose_res$seed
        )
    )
}


#' Calling CNV signature pipeline
#' @description  this pipeline integrate multiple independent steps in `VSHunter`.
#' @inheritParams cnv_derivefeatures
#' @inheritParams cnv_fitMixModels
#' @inheritParams cnv_chooseSigNumber
#' @inheritParams cnv_extractSignatures
#' @param ranks a integer vector, manually specify `ranks` to capture instead of using auto-capture feature.
#' @param tmp whether create a tmp directory to store temp result or not, default is `FALSE`.
#' @param plot_survey logical. If `TRUE`, plot rank survey.
#' @param de_novo default is `TRUE`. If set to `FALSE`, it will use reference components to
#' generate sample-by-component matrix and then extract signatures.
#' @param reference_components the object result from [cnv_fitMixModels],
#' default is `NULL`. When `de_novo` is `FALSE` and this
#' argument is `NULL`, it will use reference components from Nature Genetics paper.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a `list` contains results of NMF best rank survey, run, signature matrix, exposure list etc..
#' @export
#' @inherit cnv_extractSignatures seealso
#' @family CNV analysis functions
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## run cnv signature pipeline
#' result = cnv_pipe(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' }
cnv_pipe = function(CN_data,
                    cores = 1,
                    genome_build = c("hg19", "hg38"),
                    de_novo = TRUE,
                    reference_components = NULL,
                    min_comp = 2,
                    max_comp = 10,
                    min_prior = 0.001,
                    model_selection = "BIC",
                    nrep = 1,
                    niter = 1000,
                    nTry = 12,
                    ranks = NULL,
                    nrun = 10,
                    seed = 123456,
                    plot_survey = TRUE,
                    consensusmap_name = "nmf_consensus",
                    testRandom = FALSE,
                    nmfalg = "brunet",
                    tmp = FALSE) {
    stopifnot(exists("CN_data"),
              nTry < 100,
              is.logical(plot_survey),
              is.logical(testRandom))

    genome_build = match.arg(genome_build)
    cat("===============================================\n")
    cat("===== VSHunter: CNV Signature Pipeline ========\n")
    cat("===============================================\n")

    cat("Part 1 - Reading arguments\n")
    cat("==========\n")
    cat("Thread number                          :", cores, "\n")
    cat("Genome build                           :", genome_build, "\n")
    cat("De novo signature analysis?            :", de_novo, "\n")
    cat("Use reference components from paper?   :",
        !is.null(reference_components),
        "\n")
    cat("Minimal number of components           :", min_comp, "\n")
    cat("Maximal number of components           :", max_comp, "\n")
    cat("Minimal prior value                    :", min_prior, "\n")
    cat("Model selection strategy               :", model_selection, "\n")
    cat("Number of run for each component       :", nrep, "\n")
    cat("Maximal number of iteration            :", niter, "\n")
    cat("Maximum NMF rank in survey             :", if(is.null(ranks)) nTry else "Disable", "\n")
    cat("NMF ranks (set by hand)                :", if(!is.null(ranks)) ranks else "Disable", "\n")
    cat("Number of run for each rank            :", nrun, "\n")
    cat("Seed number                            :", seed, "\n")
    cat("Plot survey                            :", plot_survey, "\n")
    cat("Use random data in survey              :", testRandom, "\n")
    cat("Algorithm of NMF                       :", nmfalg, "\n")
    cat("Store temp result data?                :", tmp, "\n")
    cat("==========\n")

    cat("Part 2 - Derive feature distributions\n")
    cat("==========\n")
    features = cnv_derivefeatures(CN_data = CN_data,
                                  cores = cores,
                                  genome_build = genome_build)

    if (tmp) {
        tmp_dir = file.path(getwd(), "tmp")
        dir.create(tmp_dir,
                   showWarnings = FALSE,
                   recursive = TRUE)
    }

    if (tmp)
        save(features, file = file.path(tmp_dir, "VSHunter_CNV_features.RData"))

    cat("Part 3 - Fit model components (this may take some time)\n")
    cat("==========\n")
    if (!de_novo) {
        cat("Detect de_novo argument is FALSE, skip this step...\n")
    } else {
        components = cnv_fitMixModels(
            CN_features = features,
            seed = seed,
            min_comp = min_comp,
            max_comp = max_comp,
            min_prior = min_prior,
            model_selection = model_selection,
            nrep = nrep,
            niter = niter,
            cores = cores
        )

        if (tmp)
            save(components,
                 file = file.path(tmp_dir, "VSHunter_CNV_components.RData"))
    }


    cat("Part 4 - Generate a sample-by-component matrix\n")
    cat("==========\n")

    if (!de_novo) {
        cat("Using reference components...\n")
        sample_component_matrix = cnv_generateSbCMatrix(features, reference_components, cores = cores)
    } else {
        sample_component_matrix = cnv_generateSbCMatrix(features, components, cores = cores)
    }

    if (tmp)
        save(sample_component_matrix,
             file = file.path(tmp_dir, "VSHunter_CNV_SbCMatrix.RData"))

    cat("Part 5 - Capture signatures of copy number profile (this may take much time)\n")
    cat("==========\n")

    if (!is.null(ranks)) {
        cat("Specified ranks: ", ranks, "\n")
        results = cnv_extractSignatures(sample_component_matrix, nsig = ranks,
                                        seed = seed, nmfalg = nmfalg,
                                        cores = cores)
    } else {
        results = cnv_autoCaptureSignatures(
            sample_component_matrix,
            nTry = nTry,
            nrun = nrun,
            cores = cores,
            seed = seed,
            plot = plot_survey,
            consensusmap_name = consensusmap_name,
            testRandom = testRandom,
            nmfalg = nmfalg
        )
    }

    cat("Pipeline done.\n")

    if (!de_novo) {
        if (is.null(reference_components))
            message(
                "Of note, the components can obtain from https://github.com/ShixiangWang/absoluteCNVdata"
            )
        return(c(
            list(
                features = features,
                components = reference_components,
                sample_component_matrix = sample_component_matrix
            ),
            results
        ))
    } else {
        return(c(
            list(
                features = features,
                components = components,
                sample_component_matrix = sample_component_matrix
            ),
            results
        ))
    }
}



utils::globalVariables(
    c(
        "centromeres.hg19",
        "centromeres.hg38",
        "chromsize.hg19",
        "chromsize.hg38",
        "feat",
        "i"
    )
)
