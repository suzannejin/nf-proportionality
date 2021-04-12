#!/usr/local/bin/Rscript

#
# This script was used to
#       1. Get KEGG pathways from KEGG database using graphite
#       2. Match my genes, and only keep the pathways that have >=2 target genes
#   the output was stored into [cluster] nf-proportionality/data/kegg/kegg_ensembl_enzyme.rds
#                          and [github] nf-proportionality/example/data/kegg_ensembl_enzyme.rds
#
#


# arguments
args = commandArgs(trailingOnly=TRUE)
input = args[1]   # input propr object (.rds file)
output = args[2]  # output kegg pathways into .rds file


# ======================= #
# LIBRARIES AND FUNCTIONS #
# ======================= #

# load libraries
library(graphite)

# load custom functions
library(this.path)
scriptdir = this.dir()
bindir = normalizePath(file.path( this.dir(), ".." ))
source(file.path(bindir, "utils/print_msg.R"))
source(file.path(bindir, "kegg/kegg-functions.R"))

# get the kegg pathways that have at least 2 genes from my set
keggWithMygenes <- function(kegg, mygenes){

    p = nnode_per_pathway(kegg, mygenes)
    p = p[p$n_nodes >= 2,]
    kegg = kegg[p$pathway]

    return(kegg)
}


# ============ #
# PROCESS DATA #
# ============ #

# read input propr object 
print_msg("Read and process propr object")
pro = readRDS(input)
matrix = pro@matrix
rm(input)
rm(pro)

# load kegg
print_msg("Load KEGG pathways")
kegg = pathways("hsapiens", "kegg")
kegg = convertIdentifiers(kegg, "ENSEMBL")

# match genes in kegg and our data
mygenes = rownames(matrix)
mygenes = kegg2nodes(kegg, mygenes)

# keep the pathways with >= 2 mygenes
kegg = keggWithMygenes(kegg, mygenes)

# store KEGG data
print_msg("Save output with ", as.character(length(kegg)), " KEGG pathways")
saveRDS(kegg, output)
