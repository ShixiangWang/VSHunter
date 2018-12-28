#' @inherit maftools::read.maf
#' @importFrom  maftools read.maf
#' @family SNV analysis functions
#' @examples
#' \dontrun{
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- snv_readprofile(maf = laml.maf)
#' }
snv_readprofile = function(
    maf, clinicalData = NULL, removeDuplicatedVariants = TRUE,
    useAll = TRUE, gisticAllLesionsFile = NULL,
    gisticAmpGenesFile = NULL, gisticDelGenesFile = NULL,
    gisticScoresFile = NULL, cnLevel = "all", cnTable = NULL,
    isTCGA = FALSE, vc_nonSyn = NULL, verbose = TRUE
) {
    maftools::read.maf(
        maf, clinicalData = clinicalData,
        removeDuplicatedVariants = removeDuplicatedVariants,
        useAll = useAll,
        gisticAllLesionsFile = gisticAllLesionsFile,
        gisticAmpGenesFile = gisticAmpGenesFile,
        gisticDelGenesFile = gisticDelGenesFile,
        gisticScoresFile = gisticScoresFile,
        cnLevel = cnLevel, cnTable = cnTable,
        isTCGA = isTCGA, vc_nonSyn = vc_nonSyn,
        verbose = verbose
    )
}


#' @inherit maftools::trinucleotideMatrix
#' @importFrom  maftools trinucleotideMatrix
#' @family SNV analysis functions
#' @examples
#' \dontrun{
#' laml.tnm <- snv_trinucleotideMatrix(maf = laml, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
#'     prefix = 'chr', add = TRUE, useSyn = TRUE)
#' }
snv_trinucleotideMatrix = function(
    maf, ref_genome = NULL, prefix = NULL,
    add = TRUE, ignoreChr = NULL, useSyn = TRUE, fn = NULL
){
    maftools::trinucleotideMatrix(
        maf, ref_genome = ref_genome, prefix = prefix,
        add = add, ignoreChr = ignoreChr, useSyn = useSyn, fn = fn
    )
}

#' @inherit maftools::extractSignatures
#' @importFrom  maftools extractSignatures
#' @param nrun number of runs to perform on user specified rank.
#' @param nrun_Try number of runs to perform when obtain a best-fit rank for NMF.
#' @family SNV analysis functions
#' @examples
#' \dontrun{
#' laml.tnm <- snv_trinucleotideMatrix(maf = laml, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
#'     prefix = 'chr', add = TRUE, useSyn = TRUE)
#' laml.sign <- snv_extractSignatures(mat = laml.tnm, plotBestFitRes = FALSE)
#' }
snv_extractSignatures = function(mat,
                                 n = NULL,
                                 nrun = 1,
                                 nTry = 6,
                                 nrun_Try = 10,
                                 plotBestFitRes = FALSE,
                                 parallel = NULL,
                                 pConstant = NULL) {

    #------ Directly call extractSignature from maftools
    # maftools::extractSignatures(
    #     mat,
    #     n = n,
    #     nTry = nTry,
    #     plotBestFitRes = plotBestFitRes,
    #     parallel = parallel,
    #     pConstant = pConstant
    # )
    #

    #suppressPackageStartupMessages(require(NMF, quietly = TRUE))
    #transpose matrix
    mat = t(mat$nmf_matrix)

    #Validation
    zeroMutClass = names(which(rowSums(mat) == 0))

    if (length(zeroMutClass)) {
        stop(paste(
            'Warning : Found zero mutations for conversions ',
            zeroMutClass,
            sep = ''
        ))
        #Add small value to avoid zero counts (maybe not appropriate). This happens when sample size is low or in cancers with low mutation rate.
        #mat[which(rowSums(mat) == 0),] = 0.1
    }

    #To avoid error due to non-conformable arrays
    if (!is.null(pConstant)) {
        if (pConstant < 0 | pConstant == 0) {
            stop("pConstant must be > 0")
        }
        mat = mat + pConstant
    }

    #Notes:
    #Available methods for nmf decompositions are 'brunet', 'lee', 'ls-nmf', 'nsNMF', 'offset'.
    #But based 21 breast cancer signatures data, defualt brunet seems to be working close to the results.
    #Sticking with default for now.

    if (is.null(n)) {
        message('Estimating best rank..')
        if (!is.null(parallel)) {
            nmfTry = NMF::nmfEstimateRank(
                mat,
                seq(2, nTry),
                method = 'brunet',
                nrun = nrun_Try,
                seed = 123456,
                .opt = parallel
            ) #try nmf for a range of values
        } else{
            nmfTry = NMF::nmfEstimateRank(
                mat,
                seq(2, nTry),
                method = 'brunet',
                nrun = nrun_Try,
                seed = 123456
            ) #try nmf for a range of values
        }

        if (plotBestFitRes) {
            pdf(
                'nmf_consensus.pdf',
                bg = 'white',
                pointsize = 9,
                width = 12,
                height = 12,
                paper = "special"
            )
            NMF::consensusmap(nmfTry)
            dev.off()
            message('created nmf_consensus.pdf')
            #print(NMF::plot(nmfTry, 'cophenetic'))
        }

        nmf.sum = summary(nmfTry) # Get summary of estimates
        data.table::setDT(nmf.sum)
        print(nmf.sum)
        nmf.sum$diff = c(0, diff(nmf.sum$cophenetic))
        bestFit = nmf.sum[diff < 0, rank][1] #First point where cophenetic correlation coefficient starts decreasing

        plot(
            nmf.sum$rank,
            nmf.sum$cophenetic,
            axes = FALSE,
            pch = 16,
            col = "#D8B365",
            cex = 1.2,
            xlab = NA,
            ylab = NA
        )
        axis(
            side = 1,
            at = nmf.sum$rank,
            labels = nmf.sum$rank,
            lwd = 3,
            font = 2,
            cex.axis = 1.2
        )
        lines(
            x = nmf.sum$rank,
            y = round(nmf.sum$cophenetic, digits = 4),
            lwd = 3
        )
        points(
            nmf.sum$rank,
            nmf.sum$cophenetic,
            pch = 16,
            col = "#D8B365",
            cex = 1.6
        )
        axis(
            side = 2,
            at = round(nmf.sum$cophenetic, digits = 4),
            lwd = 3,
            font = 2,
            las = 2,
            cex = 1.4,
            cex.axis = 1.2
        )
        segments(
            x0 = bestFit,
            y0 = 0,
            x1 = bestFit,
            y1 = nmf.sum[rank == bestFit, cophenetic],
            lwd = 3,
            lty = 2,
            col = "maroon"
        )
        title(main = "cophenetic metric",
              adj = 0,
              font.main = 4)


        #bestFit = nmf.sum[which(nmf.sum$cophenetic == max(nmf.sum$)),'rank'] #Get the best rank based on highest cophenetic correlation coefficient
        message(
            paste(
                'Using ',
                bestFit,
                ' as a best-fit rank based on decreasing cophenetic correlation coefficient.',
                sep = ''
            )
        )
        n = as.numeric(bestFit)
    }

    if (!is.null(parallel)) {
        conv.mat.nmf = NMF::nmf(
            x = mat,
            rank = n,
            nrun = nrun,
            .opt = parallel,
            seed = 123456
        )
    } else{
        conv.mat.nmf = NMF::nmf(x = mat,
                                rank = n,
                                nrun = nrun,
                                seed = 123456)
    }

    #Signatures
    w = NMF::basis(conv.mat.nmf)
    w = apply(w, 2, function(x)
        x / sum(x)) #Scale the signatures (basis)
    colnames(w) = paste('Signature', 1:ncol(w), sep = '_')

    #Contribution
    h = NMF::coef(conv.mat.nmf)
    colnames(h) = colnames(mat) #correct colnames (seems to be mssing with low mutation load)
    #For single signature, contribution will be 100% per sample
    if (n == 1) {
        h = h / h
        rownames(h) = paste('Signature', '1', sep = '_')
    } else{
        h = apply(h, 2, function(x)
            x / sum(x)) #Scale contributions (coefs)
        rownames(h) = paste('Signature', 1:nrow(h), sep = '_')
    }


    #conv.mat.nmf.signatures.melted = melt(conv.mat.nmf.signatures)
    #levels(conv.mat.nmf.signatures.melted$X1) = colOrder

    sigs = data.table::fread(
        input = system.file('extdata', 'signatures.txt', package = 'maftools'),
        stringsAsFactors = FALSE,
        data.table = FALSE
    )
    colnames(sigs) = gsub(pattern = ' ',
                          replacement = '_',
                          x = colnames(sigs))
    rownames(sigs) = sigs$Somatic_Mutation_Type
    sigs = sigs[, -c(1:3)]
    #sigs = sigs[,1:22] #use only first 21 validated sigantures
    sigs = sigs[rownames(w), ]

    aetiology = structure(
        list(
            aetiology = c(
                "spontaneous deamination of 5-methylcytosine",
                "APOBEC Cytidine Deaminase (C>T)",
                "defects in DNA-DSB repair by HR",
                "exposure to tobacco (smoking) mutagens",
                "Unknown",
                "defective DNA mismatch repair",
                "UV exposure",
                "Unknown",
                "defects in polymerase-eta",
                "defects in polymerase POLE",
                "exposure to alkylating agents",
                "Unknown",
                "APOBEC Cytidine Deaminase (C>G)",
                "Unknown",
                "defective DNA mismatch repair",
                "Unknown",
                "Unknown",
                "Unknown",
                "Unknown",
                "defective DNA mismatch repair",
                "unknown",
                "exposure to aristolochic acid",
                "Unknown",
                "exposures to aflatoxin",
                "Unknown",
                "defective DNA mismatch repair",
                "Unknown",
                "Unknown",
                "exposure to tobacco (chewing) mutagens",
                "Unknown"
            )
        ),
        .Names = "aetiology",
        row.names = c(
            "Signature_1",
            "Signature_2",
            "Signature_3",
            "Signature_4",
            "Signature_5",
            "Signature_6",
            "Signature_7",
            "Signature_8",
            "Signature_9",
            "Signature_10",
            "Signature_11",
            "Signature_12",
            "Signature_13",
            "Signature_14",
            "Signature_15",
            "Signature_16",
            "Signature_17",
            "Signature_18",
            "Signature_19",
            "Signature_20",
            "Signature_21",
            "Signature_22",
            "Signature_23",
            "Signature_24",
            "Signature_25",
            "Signature_26",
            "Signature_27",
            "Signature_28",
            "Signature_29",
            "Signature_30"
        ),
        class = "data.frame"
    )

    message(
        'Comparing against experimentally validated 30 signatures.. (See http://cancer.sanger.ac.uk/cosmic/signatures for details.)'
    )
    #corMat = c()
    coSineMat = c()
    for (i in 1:ncol(w)) {
        sig = w[, i]
        coSineMat = rbind(coSineMat, apply(sigs, 2, function(x) {
            round(crossprod(sig, x) / sqrt(crossprod(x) * crossprod(sig)),
                  digits = 3) #Estimate cosine similarity against all 30 signatures
        }))
        #corMat = rbind(corMat, apply(sigs, 2, function(x) cor.test(x, sig)$estimate[[1]])) #Calulate correlation coeff.
    }
    #rownames(corMat) = colnames(w)
    rownames(coSineMat) = colnames(w)

    for (i in 1:nrow(coSineMat)) {
        ae = aetiology[names(which(coSineMat[i, ] == max(coSineMat[i, ]))), ]
        ae = paste0("Aetiology: ",
                    ae,
                    " [cosine-similarity: ",
                    max(coSineMat[i, ]),
                    "]")
        message(
            'Found ',
            rownames(coSineMat)[i],
            ' most similar to validated ',
            names(which(coSineMat[i, ] == max(coSineMat[i, ]))),
            '. ',
            ae,
            sep = ' '
        )
    }

    return(
        list(
            signatures = w,
            contributions = h,
            coSineSimMat = coSineMat,
            nmfObj = conv.mat.nmf
        )
    )
}


#' @inherit maftools::plotSignatures
#' @importFrom  maftools plotSignatures
#' @family SNV analysis functions
#' @return a base plot
#' @examples
#' \dontrun{
#' laml.tnm <- snv_trinucleotideMatrix(maf = laml, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
#'     prefix = 'chr', add = TRUE, useSyn = TRUE)
#' laml.sign <- snv_extractSignatures(mat = laml.tnm, plotBestFitRes = FALSE)
#' }
snv_plotSignatures = function(
    nmfRes = NULL, contributions = FALSE, color = NULL,
    patient_order = NULL, font_size = 1.2, show_title = TRUE,
    axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE,
    yaxisLim = 0.3, ...
){
    maftools::plotSignatures(
        nmfRes = nmfRes, contributions = contributions,
        color = color,
        patient_order = patient_order,
        font_size = font_size, show_title = show_title,
        axis_lwd = axis_lwd, title_size = title_size,
        show_barcodes = show_barcodes,
        yaxisLim = yaxisLim, ...
    )
}


#' @inherit maftools::signatureEnrichment
#' @importFrom  maftools signatureEnrichment
#' @family SNV analysis functions
#' @examples
#' \dontrun{
#' laml.tnm <- snv_trinucleotideMatrix(maf = laml, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
#'     prefix = 'chr', add = TRUE, useSyn = TRUE)
#' laml.sign <- snv_extractSignatures(mat = laml.tnm, plotBestFitRes = FALSE)
#' }
snv_signatureEnrichment = function(
    maf, sig_res, minMut = 5, useCNV = FALSE,
    fn = NULL
){
    maftools::signatureEnrichment(
        maf, sig_res, minMut = minMut, useCNV = useCNV,
        fn = fn
    )
}


#' Extract groups (sample clustering) from NMF run results
#'
#' X = W x H
#' W is the feature matrix, H is the sample matrix
#' After NMF run, use this function to select import features and assign groups for the two matrix
#'
#' This code is modified from [source](https://github.com/naikai/sake/blob/master/R/nmf_utils.R).
#' @param nmfObj NMF run results
#' @keywords NMF groups
#' @import NMF
#' @export
#' @examples
#' \dontrun{
#' nmf_extract_group(nmf_res)
#' }
nmf_extract_group <- function(nmfObj, type="consensus", matchConseOrder=F){
    data <- NULL
    if(type=="consensus"){
        predict.consensus <- predict(nmfObj, what="consensus")
        silhouette.consensus <- silhouette(nmfObj, what="consensus")
        # It turns out the factor levels is the NMF_assigned_groups from consensus matrix
        # that matches the original sampleNames(nmfObj) order
        # The attributes(a.predict.consensus)$iOrd is the idx order for it to match the
        # order of the samples in consensusmap(nmfObj). It is just for displaying
        # Therefore, the merged data frame sampleNames(nmfObj) + a.predict.consensus is the final
        # consensus results.
        data <- data.frame(Sample_ID=sampleNames(nmfObj),
                           nmf_subtypes = predict.consensus,
                           sil_width = signif(silhouette.consensus[, "sil_width"], 3))
        # If we want to display as we see in consensusmap, we just need to reoder everything.
        # Now re-order data to match consensusmap sample order
        if(matchConseOrder){
            sample.order <- attributes(predict.consensus)$iOrd
            data <- data[sample.order, ]
        }
    }else if(type=="samples"){
        predict.samples <- predict(nmfObj, what="samples", prob=T)
        silhouette.samples <- silhouette(nmfObj, what="samples")
        data <- data.frame(Sample_ID=names(predict.samples$predict),
                           nmf_subtypes = predict.samples$predict,
                           sil_width = signif(silhouette.samples[, "sil_width"], 3),
                           prob = signif(predict.samples$prob, 3))
    }else{
        stop(paste("Wrong type:", type, "Possible options are: 'consensus', 'samples' "))
    }
    return(data)
}
