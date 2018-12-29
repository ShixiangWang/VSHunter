library("VSHunter")

context("General test for common functions")

test_that("example data exists", {
    skip_on_cran()
    example_data = system.file("extdata", "example_cn_list.RData", package = "VSHunter")
    expect_equal(file.exists(example_data), TRUE)
})

test_that("cnv_* (analysis) separated functions work normal", {
    skip_on_cran()
    example_data = system.file("extdata", "example_cn_list.RData", package = "VSHunter")
    load(example_data)
    features = cnv_derivefeatures(tcga_segTabs, cores = 1, genome_build = "hg19")
    expect_is(features, "list")
    components = cnv_fitMixModels(features)
    expect_is(components, "list")
    sbc_matrix = cnv_generateSbCMatrix(features, components)
    expect_is(sbc_matrix, "matrix")

    skip_on_travis()
    res = cnv_autoCaptureSignatures(sbc_matrix, nTry = 3)
    expect_is(res, "list")

    # nmf_extract_groups work normal?
    nmf_groups = nmf_extract_group(res$NMF)
    expect_is(nmf_groups, "data.frame")

})

