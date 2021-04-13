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
    make_option("--kegg", 
                type="character", 
                default=NULL, 
                help="KEGG pathways (.rds file) with ENSEMBL gene names", 
                metavar="character"),
    make_option("--outdir",
                type="character",
                default=getwd(),
                help="Output directory",
                metavar="character"),
    make_option("--clique",
              action="store_true",
              default=FALSE,
              help="KEGG graph as clique")
); 
 
opt_parser = OptionParser(option_list=option_list, usage="Coexpression network vs KEGG: ROC and PR curves")
opt = parse_args(opt_parser)

# ========= #
# LIBRARIES #
# ========= #

library(biomaRt)
library(graphite)
library(propr)
library(precrec)

# load custom functions
library(this.path)
scriptdir = this.dir()
bindir = normalizePath(file.path( this.dir(), ".." ))
source(file.path(bindir, "utils/print_msg.R"))
source(file.path(bindir, "kegg/kegg-functions.R"))

# ====================== #
# LOAD AND ORGANIZE DATA
# ====================== #

# read input propr object 
print_msg("Read and process propr object")
pro = readRDS(opt$pro)
matrix = pro@matrix
rm(pro)

# load kegg
print_msg("Load KEGG pathways")
kegg = readRDS(opt$kegg)

# match genes in kegg and our data
mygenes = rownames(matrix)
mygenes = kegg2nodes(kegg, mygenes)
matrix = matrix[mygenes, mygenes]

# ============ #
# ROC ANALYSIS #
# ============ #

print_msg("Organize scores and labels")

# get reference adjacency matrix
# and convert it to label list
K = kegg2matrix(kegg, mygenes, clique=opt$clique)
labels = K[lower.tri(K)]

# get prediction scores
scores = matrix[lower.tri(matrix)]

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
