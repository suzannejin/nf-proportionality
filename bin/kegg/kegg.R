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
              help="KEGG graph as clique"),
    make_option("--minK",
              type="integer",
              default=2,
              help="Minimun pathway size (min number of target genes in the pathway)"),
    make_option("--maxK",
              type="integer",
              default=NA,
              help="Maximum pathway size (max number of target genes in the pathway)")
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
source(file.path(bindir, "utils/compute_roc.R"))
source(file.path(bindir, "kegg/kegg-functions.R"))

# ====================== #
# LOAD AND ORGANIZE DATA
# ====================== #

# read input association matrix, and genes
print_msg("Read and process propr object")
pro = readRDS(opt$pro)
matrix = pro@matrix
mygenes = rownames(matrix)
rm(pro)

# load kegg
print_msg("Load KEGG pathways")
kegg = readRDS(opt$kegg)
# get pathways containing a reasonable number of target genes
kegg = keggWithMygenes(kegg, mygenes, min=opt$minK, max=opt$maxK)

# # match genes in kegg and our data
# mygenes = kegg2nodes(kegg, mygenes)
# matrix = matrix[mygenes, mygenes]

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
compute_roc(scores, labels, opt$outdir)