
#------------------------------------------------------------------
# read copy number as a list of data.frame from data.frame or files
read_copynumbers = function(input){

}

#-------------------------------------------
# derive copy number feature distributions
derive_features = function(CN_data, cores = 1, genome_build = c("hg19", "hg38")) {
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

    if (cores > 1) {
        require(foreach)
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
# fit optimal number of mixture model components
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
#------------------------------------
# generate sample by component matrix
generate_sbcMatrix = function(CN_features,
             all_components = NULL,
             cores = 1,
             rowIter = 1000)
    {
        if (is.null(all_components))
        {
            # all_components <-
            #     readRDS(paste(this_path, "data/component_parameters.rds", sep = "/"))
            all_components <-
                readRDS(system.file("extdata", "component_parameters.rds", package = "cnPattern"))
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
# choose optimal number of signatures

# The most
# common approach is to choose the smallest rank for which cophenetic correlation coefficient
# starts decreasing. Another approach is to choose the rank for which the plot of the residual
# sum of squares (RSS) between the input matrix and its estimate shows an inflection point.
choose_nSignatures <-
    function(sample_by_component,
             nTry = 12,
             iter = 100,
             cores = 1)
    {
        message('Estimating best rank..')
        nmfalg <- "brunet"
        seed <- 77777

        estim.r <-
            NMF::nmfEstimateRank(
                t(sample_by_component),
                seq(2,nTry),
                seed = seed,
                nrun = iter,
                verbose = FALSE,
                method = nmfalg,
                .opt = list(shared.memory = FALSE, paste0("p", cores))
            )

        #--- copy from maftools and modified ---#
        nmf.sum = summary(estim.r) # Get summary of estimates
        # data.table::setDT(nmf.sum)
        print(nmf.sum)
        nmf.sum$diff = c(0, diff(nmf.sum$cophenetic))
        bestFit = nmf.sum$rank[which(nmf.sum$diff < 0)][1]
        #bestFit = nmf.sum[diff < 0, rank][1] #First point where cophenetic correlation coefficient starts decreasing

        plot(nmf.sum$rank, nmf.sum$cophenetic, axes = FALSE, pch = 16, col = "#D8B365", cex = 1.2, xlab = NA, ylab = NA)
        axis(side = 1, at = nmf.sum$rank, labels = nmf.sum$rank, lwd = 3, font = 2, cex.axis = 1.2)
        lines(x = nmf.sum$rank, y = round(nmf.sum$cophenetic, digits = 4), lwd = 3)
        points(nmf.sum$rank, nmf.sum$cophenetic, pch = 16, col = "#D8B365", cex = 1.6)
        axis(side = 2, at = round(nmf.sum$cophenetic, digits = 4), lwd = 3, font = 2, las = 2, cex = 1.4, cex.axis = 1.2)
        segments(x0 = bestFit, y0 = 0, x1 = bestFit, y1 = nmf.sum[rank == bestFit, cophenetic], lwd= 3, lty = 2, col = "maroon")
        title(main = "cophenetic metric", adj = 0, font.main = 4)


        #bestFit = nmf.sum[which(nmf.sum$cophenetic == max(nmf.sum$)),'rank'] #Get the best rank based on highest cophenetic correlation coefficient
        message(paste('Using ',bestFit, ' as a best-fit rank based on decreasing cophenetic correlation coefficient.', sep=''))
        n = as.numeric(bestFit)

        V.random <- NMF::randomize(t(sample_by_component))
        estim.r.random <-
            NMF::nmfEstimateRank(
                V.random,
                seq(2,nTry),
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
            main = "NMF Rank Survey"
        )
        return(p)

    }

#--------------------------
# extract signatures
extract_Signatures <-
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


# quantify exposure for samples
quantify_Signatures <-
    function(sample_by_component,
             component_by_signature = NULL)
    {
        if (is.null(component_by_signature))
        {
            component_by_signature <-
                readRDS(system.file("extdata", "feat_sig_mat.rds", package = "cnPattern"))
        }
        signature_by_sample <- YAPSA::LCD(
            t(sample_by_component),
            YAPSA:::normalize_df_per_dim(component_by_signature, 2)
        )
        signature_by_sample <- normaliseMatrix(signature_by_sample)
        signature_by_sample
    }
