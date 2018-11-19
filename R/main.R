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
#' @param input a \code{data.frame} or a file or a directory contains copy number profile
#' @param is_dir is \code{input} a directory?
#' @param pattern a optional regular expression used to select part of files if input is a directory, more detail please see
#' \code{list.files} function.
#' @param ignore_case logical. Should pattern-matching be case-insensitive?
#' @param sep the field separator character. Values on each line of the file are separated by this character.
#' @param cols four characters used to specify chromosome, start position,
#'  end position and copy number value, respectively. Default use names from ABSOLUTE calling result.
#' @param have_sampleCol Whether input have sample column or not.
#' This argument must be \code{TRUE} and \code{sample_col} also
#' properly be assigned when input is a file or a \code{data.frame}.
#' @param sample_col a character used to specify the sample column name.
#' @author  Shixiang Wang <w_shixiang@163.com>
#' @return a \code{list} contains absolute copy-number profile for multiple samples.
#' @importFrom utils read.csv
#' @export
read_copynumbers = function(input, is_dir = FALSE, pattern = NULL, ignore_case = FALSE, sep = "\t",
                            cols = c("Chromosome", "Start.bp", "End.bp", "modal_cn"),
                            have_sampleCol = TRUE, sample_col = "sample") {
    stopifnot(is.logical(is_dir), is.logical(have_sampleCol),
              is.character(sample_col), length(sample_col) == 1)
    if (is_dir) {
        message("Treat input as a directory...")
        if (length(input) != 1) {
            stop("Only can take one directory as input!")
        }
        # get files and exclude directories
        all.files <- list.files(path = input, pattern = pattern,
                                all.files = FALSE, recursive = FALSE,
                                ignore.case = ignore_case)
        files = all.files[!file.info(file.path(input, all.files))$isdir]
        if (length(files) == 0) stop("No files exist, please check!")
        files_path = file.path(input, files)
        files_list = list()
        for (i in seq_along(files_path)) {
            temp = read.csv(file = files_path[i], sep = sep, comment.char = "#", stringsAsFactors = FALSE)
            if (!all(cols %in% colnames(temp))) stop("not all cols are in file, please check.")
            if (have_sampleCol) {
                tempName = unique(temp[, sample_col])
                if (length(tempName) > 1) {
                    stop("When input is a directory, a file contains only one sample.")
                }
                temp = temp[, cols]
                colnames(temp) = c("chromosome", "start", "end", "segVal")
                files_list[[tempName]] = temp
            } else {
                message("Select file names as sample names.")
                temp = temp[, cols]
                colnames(temp) = c("chromosome", "start", "end", "segVal")
                files_list[[files[i]]] == temp
            }
        }
        return(files_list)
    } else if (all(is.character(input))) {
        message("Treat input as a file...")
        if (length(input) > 1) {
            stop("Muliple files are not a valid input, please use directory as input.")
        }

        if (!file.exists(input)) stop("input file not exists")
        if (!have_sampleCol) stop("When input is a file, sample column must set.")
        input = read.csv(file = input, sep = sep, comment.char = "#", stringsAsFactors = FALSE)
    }

    if (!sample_col %in% colnames(input)) stop("sample column user set not exists in input file.")
    if (!all(cols %in% colnames(input))) stop("not all cols are in file, please check.")
    samples = unique(input[, sample_col])

    res_list = list()
    for (i in seq_along(samples)) {
        tempDF = input[input[, sample_col] == samples[i], ]
        tempDF = tempDF[, cols]
        colnames(tempDF) = c("chromosome", "start", "end", "segVal")
        res_list[[samples[i]]] = tempDF
    }

    return(res_list)
}

#-------------------------------------------

#' @title  Derive copy number feature distributions
#' @description This function summarise each copy-number profile using a number of different
#' feature distributions: sigment size, breakpoint number (per ten megabase), change-point copy-number,
#' segment copy-number, breakpoint number(per chromosome arm), length of segments with oscilating
#' copy-number
#' @param CN_data a \code{QDNAseqCopyNumbers} object or a list contains multiple \code{data.frame}s
#' each one \code{data.frame} stores copy-number profile for one sample with 'chromosome', 'start', 'end' and
#' 'segVal' these four necessary columns. Of note, 'segVal' column shoule be absolute copy number values.
#' @param cores number of compute cores to run this task.
#' @param genome_build genome build version, must be one of 'hg19' or 'hg38'.
#' @author Geoffrey Macintyre, Shixiang Wang
#' @return a \code{list} contains six copy number feature distributions.
#' @import foreach doMC QDNAseq Biobase
#' @export
#'
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = derive_features(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#'}
derive_features = function(CN_data,
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
    chrlen = chrlen[chrlen$chrom %in% centromeres$chrom, ]
    if (cores > 1) {
        #require(foreach)
        requireNamespace("foreach", quietly = TRUE)
        doMC::registerDoMC(cores)

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
#' of Gaussian or mixtures of Poison distributions using the \code{flexmix} package.
#'
#' @param CN_features a \code{list} generate from \code{derive_features} function.
#' @param seed seed number.
#' @param min_comp minimal number of components to fit, default is 2.
#' @param max_comp maximal number of components to fit, default is 10.
#' @param min_prior minimal prior value, default is 0.001. Details about custom setting please
#' refer to \code{flexmix} package.
#' @param model_selection model selection strategy, default is 'BIC'.Details about custom setting please
#' refer to \code{flexmix} package.
#' @param nrep number of run times fro each value of component, keep only the solution with maximum likelihood.
#' @param niter maximal number of iteration to achive converge.
#' @param cores number of compute cores to run this task.
#' @param featsToFit integer vector used for task assignment in parallel computation. Do not change it.
#' @author Geoffrey Macintyre, Shixiang Wang
#' @return a \code{list} contain \code{flexmix} object of copy-number features.
#' @import flexmix
#' @export
#'
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = derive_features(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' ## fit mixture model  (this will take some time)
#' tcga_components = fit_mixModels(CN_features = tcga_features, cores = 1)
#' }
fit_mixModels = function(CN_features,
                         seed = 77777,
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
        doMC::registerDoMC(cores)

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
#' @param CN_features a \code{list} contains six copy number feature distributions, obtain this from
#' \code{derive_features} function.
#' @param all_components a \code{list} contain \code{flexmix} object of copy-number features, obtain this
#' from \code{fit_mixModels} function or use pre-compiled components data which come from CNV signature paper
#' https://www.nature.com/articles/s41588-018-0179-8 (set this argument as \code{NULL}).
#' @param cores number of compute cores to run this task.
#' @param rowIter step size of iteration for rows of ech CNV feature \code{data.frame}.
#' @author Geoffrey Macintyre, Shixiang Wang
#' @import doMC
#' @return a numeric sample-by-component \code{matrix}
#' @export
#'
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = derive_features(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' ## fit mixture model  (this will take some time)
#' tcga_components = fit_mixModels(CN_features = tcga_features, cores = 1)
#' ## generate a sample-by-component matrix
#' tcga_sample_component_matrix = generate_sbcMatrix(tcga_features, tcga_components, cores = 1)
#' }
generate_sbcMatrix = function(CN_features,
                              all_components = NULL,
                              cores = 1,
                              rowIter = 1000)
{
    if (is.null(all_components))
    {
        # all_components <-
        #     readRDS(paste(this_path, "data/component_parameters.rds", sep = "/"))
        message("argument all_components is not set, will download reference components.")
        message("more detail please see https://github.com/ShixiangWang/absoluteCNVdata")
        download.file(url = "https://github.com/ShixiangWang/absoluteCNVdata/raw/master/component_parameters.rds",
                      destfile = "Nat_Gen_component_parameters.rds")
        all_components <-
            readRDS("Nat_Gen_component_parameters.rds")
    }

    if (cores > 1) {
        #require(foreach)
        requireNamespace("foreach", quietly = TRUE)
        feats = c("segsize",
                  "bp10MB",
                  "osCN",
                  "changepoint",
                  "copynumber",
                  "bpchrarm")
        doMC::registerDoMC(cores)

        full_mat = foreach(feat = feats, .combine = cbind) %dopar% {
            calculateSumOfPosteriors(
                CN_features[[feat]],
                all_components[[feat]],
                feat,
                rowIter = rowIter,
                cores = cores
            )
        }
    } else {
        full_mat <- cbind(
            calculateSumOfPosteriors(CN_features[["segsize"]], all_components[["segsize"]], "segsize"),
            calculateSumOfPosteriors(CN_features[["bp10MB"]], all_components[["bp10MB"]], "bp10MB"),
            calculateSumOfPosteriors(CN_features[["osCN"]], all_components[["osCN"]], "osCN"),
            calculateSumOfPosteriors(CN_features[["changepoint"]], all_components[["changepoint"]], "changepoint"),
            calculateSumOfPosteriors(CN_features[["copynumber"]], all_components[["copynumber"]], "copynumber"),
            calculateSumOfPosteriors(CN_features[["bpchrarm"]], all_components[["bpchrarm"]], "bpchrarm")
        )
    }

    rownames(full_mat) <- unique(CN_features[["segsize"]][, 1])
    full_mat[is.na(full_mat)] <- 0
    full_mat
}


#------------------------------------
#' @title Choose optimal number of signatures
#' @description This function use \code{NMF} package to evaluate the optimal number of signatures. The most
# common approach is to choose the smallest rank for which cophenetic correlation coefficient
# starts decreasing (Used by this function). Another approach is to choose the rank for which the plot of the residual
# sum of squares (RSS) between the input matrix and its estimate shows an inflection point.
#'
#' @param sample_by_component a sample-by-component \code{matrix}, generate from \code{generate_sbcMatrix} function.
#' @param nTry the maximal tried number of signatures, default is 12. Of note, this value should far less than number
#' of features or samples.
#' @param nrun a numeric giving the number of run to perform for each value in range of 2 to \code{nTry}, default is 50.
#' According to \code{NMF} package documentation, nrun set to 50 is enough to achieve robust result.
#' @param cores number of compute cores to run this task.
#' @param seed seed number.
#' @param plot \code{logical}. If \code{TRUE}, plot best rank survey.
#' @author Geoffrey Macintyre, Shixiang Wang
#' @import NMF
#' @import grDevices
#'
#' @return a \code{list} contains information of NMF run and rank survey.
#' @export
#'
#' @examples
#' \dontrun{
#' #' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = derive_features(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' ## fit mixture model  (this will take some time)
#' tcga_components = fit_mixModels(CN_features = tcga_features, cores = 1)
#' ## generate a sample-by-component matrix
#' tcga_sample_component_matrix = generate_sbcMatrix(tcga_features, tcga_components, cores = 1)
#' ## optimal rank survey
#'  tcga_sig_choose = choose_nSignatures(tcga_sample_component_matrix,
#'  nrun = 10, cores = 1, plot = FALSE)
#' }
#'
choose_nSignatures <-
    function(sample_by_component,
             nTry = 12,
             nrun = 50,
             cores = 1,
             seed = 77777,
             plot = TRUE)
    {
        message('Estimating best rank..')
        nmfalg <- "brunet"

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


        #--- copy from maftools and modified ---#
        nmf.sum = summary(estim.r) # Get summary of estimates
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

        if (plot) {
            # message('Creating nmf consensusmap plot...')
            # pdf('nmf_consensus.pdf', bg = 'white', pointsize = 9, width = 12, height = 12, paper = "special")
            # NMF::consensusmap(estim.r)
            # dev.off()


            message('Creating nmf rank survey plot...')
            p <- NMF::plot(
                estim.r,
                estim.r.random,
                what = c(
                    "cophenetic",
                    "dispersion",
                    "sparseness",
                    "silhouette",
                    "residuals",
                    "rss"
                ),
                xname = "Observed",
                yname = "Randomised",
                main = "NMF Rank Survey"
            )

            print(p)

            pdf(
                'nmf_rank_survey.pdf',
                bg = 'white',
                pointsize = 9,
                width = 8,
                height = 6,
                paper = "special"
            )
            print(p)
            dev.off()

            return(
                list(
                    nmfEstimate = estim.r,
                    bestRank = n,
                    survey = nmf.sum,
                    survey_plot = p,
                    seed = seed
                )
            )

        }

        return(list(
            nmfEstimate = estim.r,
            bestRank = n,
            survey = nmf.sum,
            seed = seed
        ))
    }

#--------------------------
# extract signatures
#' @title Extract signature based on specified rank value
#' @inheritParams choose_nSignatures
#' @param nsig specification of the factorization rank.
#' @param nmfalg specification of the NMF algorithm.
#' @author Geoffrey Macintyre, Shixiang Wang
#' @import NMF
#'
#' @return a object of \code{NMF} run.
#' @export
#'
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = derive_features(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' ## fit mixture model  (this will take some time)
#' tcga_components = fit_mixModels(CN_features = tcga_features, cores = 1)
#' ## generate a sample-by-component matrix
#' tcga_sample_component_matrix = generate_sbcMatrix(tcga_features, tcga_components, cores = 1)
#' ## optimal rank survey
#'  tcga_sig_choose = choose_nSignatures(tcga_sample_component_matrix, nrun = 10,
#'  cores = 1, plot = FALSE)
#'  tcga_signatures = extract_Signatures(tcga_sample_component_matrix, nsig = 3, cores = 1)
#' }
extract_Signatures <-
    function(sample_by_component,
             nsig,
             seed = 77777,
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
            .opt = paste0("vp", cores)
        )
    }

#---------------------------------------------------------------------------
#' @title quantify exposure for samples using Linear Combination Decomposition (LCD)
#'
#' @inheritParams choose_nSignatures
#' @param component_by_signature a componet by signature matrix, default is \code{NULL},
#' it will use pre-compiled data from CNV signature paper
#' https://www.nature.com/articles/s41588-018-0179-8
#' @author Geoffrey Macintyre, Shixiang Wang
#'
#' @return a \code{list} contains absolute/relative exposure.
#' @import YAPSA
#' @export
#'
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = derive_features(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' ## fit mixture model  (this will take some time)
#' tcga_components = fit_mixModels(CN_features = tcga_features, cores = 1)
#' ## generate a sample-by-component matrix
#' tcga_sample_component_matrix = generate_sbcMatrix(tcga_features, tcga_components, cores = 1)
#' ## optimal rank survey
#'  tcga_sig_choose = choose_nSignatures(tcga_sample_component_matrix,
#'  nrun = 10, cores = 1, plot = FALSE)
#'  tcga_signatures = extract_Signatures(tcga_sample_component_matrix, nsig = 3, cores = 1)
#'  w = NMF::basis(tcga_signatures) # signature matrix
#'  tcga_exposure = quantify_Signatures(sample_by_component =
#'  tcga_sample_component_matrix, component_by_signature = w)
#' }
quantify_Signatures <-
    function(sample_by_component,
             component_by_signature = NULL)
    {
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
#' @description  this is a wrapper of \code{choose_nSignatures}, \code{extract_Signatures}
#' and \code{quantify_Signatures} these three functions.
#'
#' @inheritParams choose_nSignatures
#' @author Geoffrey Macintyre, Shixiang Wang
#' @return a \code{list} contains results of NMF best rank survey, run, signature matrix, exposure list etc..
#' @import doMC NMF YAPSA
#' @export
#'
#' @examples
#' \dontrun{
#' ## load example copy-number data from tcga
#' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' ## generate copy-number features
#' tcga_features = derive_features(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' ## fit mixture model  (this will take some time)
#' tcga_components = fit_mixModels(CN_features = tcga_features, cores = 1)
#' ## generate a sample-by-component matrix
#' tcga_sample_component_matrix = generate_sbcMatrix(tcga_features, tcga_components, cores = 1)
#' ## optimal rank survey13
#' tcga_results = autoCapture_Signatures(tcga_sample_component_matrix, nrun=10, cores = 1)
#' }
autoCapture_Signatures = function(sample_by_component,
                                  nTry = 12,
                                  nrun = 50,
                                  cores = 1,
                                  seed = 77777,
                                  plot = TRUE) {
    choose_res = choose_nSignatures(sample_by_component, nTry, nrun, cores, seed, plot = plot)
    NMF_res = extract_Signatures(sample_by_component,
                                 nsig = choose_res$bestRank,
                                 cores = cores)
    w = NMF::basis(NMF_res)
    #h = NMF::coef(NMF_res)
    exposure = quantify_Signatures(sample_by_component = sample_by_component,
                                   component_by_signature = w)
    message("Done.")
    return(
        list(
            NMF = NMF_res,
            signature = w,
            exposure = exposure,
            bestRank = choose_res$n,
            survey = choose_res$survey,
            survey_plot = choose_res$survey_plot,
            seed = choose_res$seed
        )
    )
}

utils::globalVariables(c("centromeres.hg19", "centromeres.hg38", "chromsize.hg19", "chromsize.hg38",
                         "feat", "i"))
