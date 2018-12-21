tcga_test_res = cnv_pipe(tcga_segTabs, cores = 4, genome_build = "hg19")

load("~/Documents/GitHub/CNAResolver/report_Rmarkdown/data/abs_copynumber_prad.RData")

# de novo signature capture
tcga_prad_snp = cnv_pipe(CN_data = abs_copynumber_prad, genome_build = "hg38", cores = 4, nrun = 50)

# refernce signature capture
tcga_prad_snp_ref = cnv_pipe(CN_data = abs_copynumber_prad, genome_build = "hg38",
                             de_novo = FALSE,
                             cores = 4, nrun = 50)

save(tcga_prad_snp, file = "~/Documents/tcga_prad_cnv_signature_res.RData")


prad_res = readRDS("~/Documents/GitHub/CNAResolver/report_Rmarkdown/data/tcga_prad_cnv_signature_res.rds")

cnv_plotSignatures(prad_res)
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
cnv_plotFeatureDistribution(prad_res$features)
cnv_plotMixComponents(prad_res$features, prad_res$components)

prad_frac =  cnv_getLengthFraction(abs_copynumber_prad, genome_build = "hg38")

cnv_plotDistributionProfile(prad_frac, rm_normal = TRUE, mode = "cd", genome_build = "hg38")



q = ggplot_build(p)
q$data[[1]]
ggplot(q$data[[1]], aes(x = factor(x, levels = as.character(1:22)), y = count,
                        fill = factor(fill, levels = c("#F8766D","#00BA38", "#619CFF")))) + geom_bar(stat = "identity")

all(q$data[[1]]$y ==q$data[[1]]$ymax)

q$data[[1]]["y"] = q$data[[1]]["y"] / 1000
q$data[[1]]["count"] = q$data[[1]]["count"] / 1000
q$data[[1]]["ymin"] = q$data[[1]]["ymin"] / 1000
q$data[[1]]["ymax"] = q$data[[1]]["ymax"] / 1000

q <- ggplot_gtable(q)
ggdraw(q) + scale_y_continuous(limits = c(0, 10))
##q
##plot(q) + theme_cowplot()
q1= ggplot_build(ggdraw(q))
