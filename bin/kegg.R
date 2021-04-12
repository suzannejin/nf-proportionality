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
    make_option("--output",
                type="character",
                default=NULL,
                help="Output file with aucROC and aucPR values",
                metavar="character"),
    make_option("--image",
                type="character",
                default=NULL,
                help="Output file ROC and PR curves",
                metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list, usage="Normalize and/or log-transform data for proportionality computation")
opt = parse_args(opt_parser)

# ========= #
# LIBRARIES #
# ========= #

library(biomaRt)
library(graphite)
library(propr)
library(precrec)


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

print_msg("Organize scores and labels")

# get reference adjacency matrix
# and convert it to label list
edges = kegg2edges(kegg, mygenes)
K = edges2matrix(edges, genes_df$SYMBOL)
labels = K[lower.tri(K)]

# get prediction scores
scores = matrix[lower.tri(matrix)]

# ============ #
# ROC ANALYSIS #
# ============ #

print_msg("Compute ROC and PR curves for ", length(labels), " items")


# create data
dat = mmdata(l_scores, l_labels, modnames = names(l_scores))
print(dat)

# Calculate ROC and Precision-Recall curves for multiple models
mscurves <- evalmod(dat)

# Plot and save ROC and Precision-Recall curves 
jpeg(file=opt$output)
autoplot(mscurves)
dev.off()

# Retrieve a dataframe of AUCs
aucs.df <- auc(mscurves)

# Use knitr::kable to display the result in a table format
knitr::kable(aucs.df)