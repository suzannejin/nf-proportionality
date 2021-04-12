#!/usr/local/bin/Rscript

# ========= #
# ARGUMENTS #
# ========= #

library(optparse)

option_list = list(
  make_option("--data", 
              type="character", 
              default=NULL, 
              help="Input count data file (in .rds format)", 
              metavar="character"),
  make_option("--output",
              type="character",
              default=NULL,
              help="Processed data (in .rds format)",
              metavar="character"),
  make_option("--norm", 
              type="character", 
              default="NA", 
              help="Data normalization. Default=NA. Options=[rpkm,tpm,tmm,NA]", 
              metavar="character"),
  make_option("--ivar", 
              type="character", 
              default="clr", 
              help="Data transformation for propr. Default=clr(centered log-ratio trasformation). Options=[clr,log2,NA]", 
              metavar="character"),
  make_option("--tissue",
              type="character",
              default=NULL,
              help="Tissue name. Use samples from this tissue only",
              metavar="character"),
  make_option("--donor",
              type="character",
              default=NULL,
              help="Donor ID. Use samples from this donor only",
              metavar="character"),
  make_option("--test",
              action="store_true",
              default=FALSE,
              help="Run a test on a small subset"),
  make_option("--ntest",
              type="integer",
              default=100,
              help="n: defines the number of gene features of the test"),
  make_option("--seed",
              type="integer",
              default=0,
              help="Random seed for the sampling of test gene features")   
); 
 
opt_parser = OptionParser(option_list=option_list, usage="Normalize and/or log-transform data for proportionality computation")
opt = parse_args(opt_parser)

# ========= #
# LIBRARIES #
# ========= #

library(recount)
library(zCompositions)
library(propr)
library(edgeR)

# # TODO use this when the docker is updated
# library(this.path)
# scriptdir = this.dir()
# source(file.path(scriptdir, "print_msg.R"))
print_msg <- function(...){

    vector = c(...)
    msg = paste(vector, collapse="")
    msg = paste(Sys.time(), msg, "\n", collapse="\t")
    cat(msg)
}

# ================================ #
# READ AND PROCESS EXPRESSION DATA #
# ================================ #

print_msg("Reading input files")

# read gtex data
rse_gene = readRDS(file.path(opt$data))
# get count matrix
m = assay(rse_gene)


# run test on randomly selected n gene features
if (opt$test){
    print_msg("Warning: only ", opt$ntest, " randomly selected genes will be used for testing purposes")
    features = seq(1:nrow(m))
    set.seed(opt$seed); pos = sample(features, opt$ntest)
    rse_gene = rse_gene[pos,]
    m = m[pos,]
}

print_msg("Starting processing count data of dim[", nrow(m), ",", ncol(m), "]")

# normalize the counts if required
# https://www.reneshbedre.com/blog/expression_units.html
# https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
if (opt$norm != "NA"){
    print_msg("Normalizing gene expression data[", opt$norm, "]")
}
if (opt$norm == "rpkm"){
  m = getRPKM(rse_gene)
}else if (opt$norm == "tpm"){
  m = getTPM(rse_gene)
}else if (opt$norm == "tmm"){
  # Note that normalization in edgeR is model-based, and the original read counts are not themselves transformed. This means that users should not transform the read counts in any way
  # before inputing them to edgeR. For example, users should not enter RPKM or FPKM values to edgeR in place of read counts. Such quantities will prevent edgeR from correctly
  # estimating the mean-variance relationship in the data, which is a crucial to the statistical
  # strategies underlying edgeR. Similarly, users should not add artificial values to the counts
  # before inputing them to edgeR

  # For further information check: https://www.biostars.org/p/84087/ and https://www.biostars.org/p/317701/
  y = DGEList(counts=m)
  y = calcNormFactors(y)
  m = cpm(y)
}

# log transformation if required
# https://www.biostars.org/p/100926/
# https://www.biostars.org/p/242573/
# https://blog.qbaseplus.com/seven-tips-for-bio-statistical-analysis-of-gene-expression-data
if (opt$ivar != "NA"){
    print_msg("Log transforming data [", opt$ivar, "]")
}
if (opt$ivar == "clr"){
  M = t(m)
  M = cmultRepl(M,method="CZM",label=0,output="p-counts")  # replace zeros
  M = propr:::proprCLR(M)  # get clr
  m = t(M)
}else if (opt$ivar == 'log2'){
  m = log2(m+1)
}

# 

# subset data, if required
subsetname = NA
if (!is.null(opt$tissue)){
    subsetname = opt$tissue
}else if (!is.null(opt$donor)){
    subsetname = opt$donor
}
if (!is.na(subsetname)){
    print_msg("Subset: only use samples from ", subsetname)

    # extract samples from the given tissue
    if (!is.null(opt$tissue)){
        s=which(colData(rse_gene)[,"smtsd"]==subsetname)

    # extract samples from the given donor
    }else if(!is.null(opt$donor)){
        donors = colData(rse_gene)[,"sampid"]
        donors = stringr::str_extract(donors, "[A-Z]+-[0-9A-Z]+")
        colData(rse_gene)[,"sampid"] = donors

        s=which(colData(rse_gene)[,"sampid"]==subsetname)
    }

    # subset matrix
    m = m[,s]
}

# write output matrix : this matrix will be used as input for propr
print_msg("Saving processed expression count data")
M = t(m)
saveRDS(M, file=opt$output)
