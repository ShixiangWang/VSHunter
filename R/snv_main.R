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
#' @family SNV analysis functions
#' @examples
#' \dontrun{
#' laml.tnm <- snv_trinucleotideMatrix(maf = laml, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
#'     prefix = 'chr', add = TRUE, useSyn = TRUE)
#' laml.sign <- snv_extractSignatures(mat = laml.tnm, plotBestFitRes = FALSE)
#' }
snv_extractSignatures = function(
    mat, n = NULL, nTry = 6, plotBestFitRes = FALSE,
    parallel = NULL, pConstant = NULL
){
    maftools::extractSignatures(
        mat, n = n, nTry = nTry, plotBestFitRes = plotBestFitRes,
        parallel = parallel, pConstant = pConstant
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

