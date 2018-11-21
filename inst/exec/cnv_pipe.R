#!/usr/bin/env Rscript
# Copyright 2018 Shixiang Wang <w_shixiang@163.com>

if (!require("optparse")) {
    message("package optparse is not installed, try installing it")
    install.package("optparse", dependencies = TRUE)
    suppressPackageStartupMessages(library("optparse"))
}

if (!require("devtools")) {
    message("package devtools is not installed, try installing it")
    install.package("devtools", dependencies = TRUE)
    suppressPackageStartupMessages(library("devtools"))
}

if (!require("VSHunter")) {
    message("package VSHunter is not installed, try installing it")
    devtools::install_github("ShixiangWang/VSHunter", dependencies = TRUE)
    suppressPackageStartupMessages(library("VSHunter"))
}

# specify desired options in a list
option_list = list(
    make_option(c("-i", "--input"), type = "character", default = NULL, help = "Path to input file or directory"),
    make_option("--RDS", action = "store_true", default = FALSE, help = "Input is a RDS file which is a list contain cnv of multiple samples"),
    make_option(c("-o", "--output"), type = "character", default = "cns.RData", help = "Output RData file name, [default %default], directory can be specified as [directory/filename]"),
    make_option(c("-d", "--directory"), action = "store_true", default = FALSE, help = "Input is a directory"),
    make_option("--pattern", type = "character", default = NULL, help = "An optional regular expression used to select part of files if input is a directory [default %default]"),
    make_option(c("-f", "--sep"), type = "character", default = "\t", help = "The field separator character of file [default is \\t]"),
    make_option("--colnames", type = "character", default = "Chromosome,Start.bp,End.bp,modal_cn", help = "Comma [,] seperated and ordered colnames to chromosome, start position, end position and copy number value, respectively. [default %default]"),
    make_option("--noSampleCol", action = "store_true", default = FALSE, help = "Input file has no sample column [default %default]"),
    make_option("--sampleCol", type = "character", default = "sample", help = "Sample column name [default %default]"),
    make_option("--refversion", type = "character", default = "hg19", help = "Genome build version, must be one of hg19 or hg38, [default %default]"),
    make_option("--mincomp", type = "integer", default = 2, help = "Minimal number of components to fit [default %default]"),
    make_option("--maxcomp", type = "integer", default = 10, help = "Maximal number of components to fit [default %default]"),
    make_option("--minprior", type = "double", default = 0.001, help = "Minimal prior value [default %default]"),
    make_option("--selectMethod", type = "character", default = "BIC", help = "Model selection strategy, [default %default]"),
    make_option("--nrep", type = "integer", default = 1, help = "Number of run times fro each value of component [default %default]"),
    make_option("--niter", type = "integer", default = 1000, help = "Maximal number of iteration to achive converge [default %default]"),
    make_option("--nTry", type = "integer", default = 12, help = "Maximal tried number of signatures, [default %default]. Of note, this value should far less than number of features or samples"),
    make_option("--nrun", type = "integer", default = 10, help = "The number of run to perform for each value [default %default]"),
    make_option("--seed", type = "integer", default = 123456, help = "Seed number [default %default]"),
    make_option("--plot", action = "store_true", default = FALSE, help = "Plot best rank survey"),
    make_option("--random", action = "store_true", default = FALSE, help = "Generate random data from input to test measurements"),
    make_option("--nmfalg", type = "character", default = "brunet", help = "Specification of the NMF algorithm [default %default]"),
    make_option(c("-T", "--thread"), type = "integer", default = 1, help = "Number of cores to run this program [default %default]")
)

# get command line options
opt <- parse_args(OptionParser(option_list=option_list))

# opt = parse_args(OptionParser(option_list=option_list), args = c("-d", "--plot", "--random", "--input=wsx", "--colnames=Chromosome, Start, End, SegVal"))

if (is.null(opt$input)) stop("input must be specified with -i or --input option")

if (opt$RDS) {
    if (!file.exists(opt$input)) stop("input file does not exist")
    cnp = readRDS(opt$input)
} else {
    # check some option
    cols = trimws(unlist(strsplit(opt$colnames, split = ",")))
    if (length(cols) != 4) stop("--colnames option should follow by 4 character sperated by comma.")
    sample_col = trimws(opt$sampleCol)

    # read input
    cnp = cnv_readprofile(input = opt$input, is_dir = opt$directory, pattern = opt$pattern,
                          sep = opt$sep, cols = cols, have_sampleCol = !opt$noSampleCol, sample_col = sample_col)
}


# run pipeline
res = cnv_pipe(CN_data = cnp, cores = opt$thread, genome_build = opt$refversion,
               min_comp = opt$mincomp, max_comp = opt$maxcomp, min_prior = opt$minprior,
               model_selection = opt$selectMethod, nrep = opt$nrep, niter = opt$niter,
               nTry = opt$nTry, nrun = opt$nrun, seed = opt$seed, plot_survey = opt$plot,
               testRandom = opt$random, nmfalg = opt$nmfalg)

cat("Result RData file is saved to ", opt$output, "\n")
save(res, file = opt$output, compress = TRUE)
