tcga_test_res = cnv_pipe(tcga_segTabs, cores = 4, genome_build = "hg19")

load("inst/test/abs_copynumber_prad.RData")

# de novo signature capture
tcga_prad_snp = cnv_pipe(CN_data = abs_copynumber_prad, genome_build = "hg38", cores = 4, nrun = 50)

# refernce signature capture
tcga_prad_snp_ref = cnv_pipe(CN_data = abs_copynumber_prad, genome_build = "hg38",
                             de_novo = FALSE,
                             cores = 4, nrun = 50)

save(tcga_prad_snp, file = "~/Documents/tcga_prad_cnv_signature_res.RData")


load("~/Documents/tcga_prad_cnv_signature_res.RData")

# prad_features = derive_features(CN_data = abs_copynumber_prad, cores = 1, genome_build = "hg38")
# prad_components = fit_mixModels(CN_features = prad_features, cores = 1, min_comp = 2, max_comp = 10)
# prad_sample_component_matrix = generate_sbcMatrix(prad_features, prad_components, cores = 8)
# prad_sig_choose = choose_nSignatures(prad_sample_component_matrix, nrun = 10, cores = 8)
# prad_signatures = extract_Signatures(prad_sample_component_matrix, nsig = 3, cores = 8)
#
#
# lapply(prad_features2, function(x) any(is.na(x)))

cnv_plotSignatures(tcga_signatures)
cnv_plotSignatures(tcga_signatures, contributions = TRUE)
cnv_plotSignatures(tcga_signatures, contributions = TRUE, show_barcodes = T)

cnv_plotSignatures(tcga_results)
cnv_plotSignatures(tcga_results, contributions = TRUE)

cnv_plotSignatures(tcga_prad_snp$NMF)
cnv_plotFeatureDistribution(tcga_features)
cnv_plotMixComponents(tcga_features, tcga_components)

#----------
cnv_plotSignatures(tcga_prad_snp)
cnv_plotSignatures(tcga_prad_snp, contributions = TRUE)
cnv_plotMixComponents(tcga_prad_snp)

prad_frac =  cnv_getLengthFraction(abs_copynumber_prad, genome_build = "hg38")
