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


# Code Part ---------------------------------------------------------------

# locationOfThisScript = function()
#     # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
# {
#     this.file = NULL
#     # This file may be 'sourced'
#     for (i in-(1:sys.nframe())) {
#         if (identical(sys.function(i), base::source))
#             this.file = (normalizePath(sys.frame(i)$ofile))
#     }
#
#     if (!is.null(this.file))
#         return(dirname(this.file))
#
#     # But it may also be called from the command line
#     cmd.args = commandArgs(trailingOnly = FALSE)
#     cmd.args.trailing = commandArgs(trailingOnly = TRUE)
#     cmd.args = cmd.args[seq.int(from = 1,
#                                 length.out = length(cmd.args) - length(cmd.args.trailing))]
#     res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
#
#     # If multiple --file arguments are given, R uses the last one
#     res = tail(res[res != ""], 1)
#     if (0 < length(res))
#         return(dirname(res))
#
#     # Both are not the case. Maybe we are in an R GUI?
#     return(NULL)
# }
# this_path <- locationOfThisScript()
# source(paste(this_path, "helper_functions.R", sep = "/"))

quantifySignatures <-
    function(sample_by_component,
             component_by_signature = NULL)
    {
        if (is.null(component_by_signature))
        {
            component_by_signature <-
                readRDS(paste(this_path, "data/feat_sig_mat.rds", sep = "/"))
        }
        signature_by_sample <- YAPSA::LCD(
            t(sample_by_component),
            YAPSA:::normalize_df_per_dim(component_by_signature, 2)
        )
        signature_by_sample <- normaliseMatrix(signature_by_sample)
        signature_by_sample
    }

generateSignatures <-
    function(sample_by_component,
             nsig,
             seed = 77777,
             nmfalg = "brunet",
             cores = 1)
    {
        NMF::nmf(
            t(sample_by_component),
            nsig,
            seed = seed,
            nrun = 1000,
            method = nmfalg,
            .opt = paste0("p", cores)
        )
    }

chooseNumberSignatures <-
    function(sample_by_component,
             outfile = "numSigs.pdf",
             min_sig = 3,
             max_sig = 12,
             iter = 100,
             cores = 1)
    {
        nmfalg <- "brunet"
        seed <- 77777

        estim.r <-
            NMF::nmfEstimateRank(
                t(sample_by_component),
                min_sig:max_sig,
                seed = seed,
                nrun = iter,
                verbose = FALSE,
                method = nmfalg,
                .opt = list(shared.memory = FALSE, paste0("p", cores))
            )

        V.random <- NMF::randomize(t(sample_by_component))
        estim.r.random <-
            NMF::nmfEstimateRank(
                V.random,
                min_sig:max_sig,
                seed = seed,
                nrun = iter,
                verbose = FALSE,
                method = nmfalg,
                .opt = list(shared.memory = FALSE, paste0("p", cores))
            )

        p <- NMF::plot(
            estim.r,
            estim.r.random,
            what = c("cophenetic", "dispersion", "sparseness", "silhouette"),
            xname = "Observed",
            yname = "Randomised",
            main = ""
        )
        pdf(file = outfile,
            width = 10,
            height = 10)
        p
        dev.off()

        return(p)

    }

extractCopynumberFeatures <- function(CN_data, cores = 1, genome_build = c("hg19", "hg38")){

    genome_build = match.arg(genome_build)
    # get chromosome lengths and centromere locations
    if (genome_build == "hg19") {
        data("chromsize.hg19", package = "cnPattern")
        data("centromeres.hg19", package = "cnPattern")
        chrlen = chromsize.hg19
        centromeres = centromeres.hg19
    } else {
        data("chromsize.hg38", package = "cnPattern")
        data("centromeres.hg38", package = "cnPattern")
        chrlen = chromsize.hg38
        centromeres = centromeres.hg38
    }

    # chrlen <-
    #     read.table(
    #         paste(this_path, "data/hg19.chrom.sizes.txt", sep = "/"),
    #         sep = "\t",
    #         stringsAsFactors = F
    #     )[1:24, ]
    #
    # gaps <-
    #     read.table(
    #         paste(this_path, "data/gap_hg19.txt", sep = "/"),
    #         sep = "\t",
    #         header = F,
    #         stringsAsFactors = F
    #     )
    # centromeres <- gaps[gaps[, 8] == "centromere", ]

    if (cores > 1) {
        require(foreach)
        doMC::registerDoMC(cores)

        temp_list = foreach::foreach(i = 1:6) %dopar% {
            if (i == 1) {
                list(segsize = getSegsize(CN_data))
            } else if (i == 2) {
                list(bp10MB = getBPnum(CN_data, chrlen))
            } else if (i == 3) {
                list(osCN = getOscilation(CN_data, chrlen))
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

fitMixtureModels <-
    function(CN_features,
             seed = 77777,
             min_comp = 2,
             max_comp = 10,
             min_prior = 0.001,
             model_selection = "BIC",
             nrep = 1,
             niter = 1000,
             cores = 1,
             featsToFit = seq(1, 6))
    {
        if (cores > 1) {
            require(foreach)

            doMC::registerDoMC(cores)

            temp_list = foreach(i = 1:6) %dopar% {
                if (i == 1 & i %in% featsToFit) {
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

generateSampleByComponentMatrix <-
    function(CN_features,
             all_components = NULL,
             cores = 1,
             rowIter = 1000,
             subcores = 2)
    {
        if (is.null(all_components))
        {
            all_components <-
                readRDS(paste(this_path, "data/component_parameters.rds", sep = "/"))
        }

        if (cores > 1) {
            require(foreach)

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
                    cores = subcores
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
