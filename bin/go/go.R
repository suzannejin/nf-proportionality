#!/usr/local/bin/Rscript

# ========= #
# ARGUMENTS #
# ========= #

library(optparse)

option_list = list(
    make_option("--pro", 
                type="character", 
                default=NULL, 
                help="Input propr object (.rds file)", 
                metavar="character"),
    make_option("--outdir",
                type="character",
                default=getwd(),
                help="Output directory",
                metavar="character"),
    make_option("--minK",
              type="integer",
              default=2,
              help="Minimun GO term size (min number of target genes annotated with the same GO)"),
    make_option("--maxK",
              type="integer",
              default=NA,
              help="Maximum GO term size (max number of target genes annotated with the same GO)")
); 
 
opt_parser = OptionParser(option_list=option_list, usage="Coexpression network vs GO: ROC and PR curves")
opt = parse_args(opt_parser)

# ========= #
# LIBRARIES #
# ========= #

library(precrec)

# load custom functions
library(this.path)
scriptdir = this.dir()
bindir = normalizePath(file.path( this.dir(), ".." ))
source(file.path(bindir, "utils/print_msg.R"))
source(file.path(bindir, "go/go-functions.R"))

# ====================== #
# LOAD AND ORGANIZE DATA
# ====================== #

# read input propr object 
print_msg("Read and process propr object")
pro = readRDS(opt$pro)
matrix = pro@matrix
mygenes = rownames(matrix)
rm(pro)

# reference GO graph
print_msg("Extract reference GO graph")
K = getK2(mygenes, 
            keytype="ENSEMBL", 
            ont=c("BP"),
            evi=c("EXP", "IDA", "IPI", "IMP", "IGI", "TAS", "IC"),
            minK=opt$minK,
            maxK=opt$maxK)
K = K2adjmatrix(K,mygenes)

# ============ #
# ROC ANALYSIS #
# ============ #

# scores and labels
print_msg("Organize scores and labels")
scores = matrix[lower.tri(matrix)]
labels = K[lower.tri(K)]

print_msg("Compute ROC and PR curves for ", length(labels), " items")

# create data
dat = mmdata(scores, labels)
print(dat)

# Calculate ROC and Precision-Recall curves for multiple models
mscurves <- evalmod(dat)

# Plot and save ROC and Precision-Recall curves 
jpeg(file=file.path(opt$outdir, "curve.jpg"))
ggplot2::autoplot(mscurves)
dev.off()

# Retrieve a dataframe of AUCs
aucs.df <- auc(mscurves)
roc = aucs.df[which(aucs.df$curvetypes=="ROC"),"aucs"]
pr = aucs.df[which(aucs.df$curvetypes=="PRC"),"aucs"]
writeLines(as.character(roc), file.path(opt$outdir, "roc.txt"))
writeLines(as.character(pr), file.path(opt$outdir, "pr.txt"))
