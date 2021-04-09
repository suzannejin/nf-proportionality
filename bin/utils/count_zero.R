#
# This script reads a count data
# and counts the percentage of zeros for each tissue & donor
#

library(recount)
library(matrixStats)

# arguments
args = commandArgs(trailingOnly=TRUE)
input = args[1]
output = args[2] 

# read count data
load(input)   # rse_gene
if (exists("rse_gene")){
    m = assay(rse_gene)
}else if(exists("rse_gene2")){
    rse_gene = rse_gene2
    m = assay(rse_gene)
}else{
    stop("ERROR: make sure that you give the correct input")
}

# get the percentage of genes that have 0 expression in at least 100*q % of samples
get_perc_gene <- function(counts, q=.8){
    rowq = rowQuantiles(counts, probs=q)
    pos = which(rowq == 0)
    perc = length(pos) / length(rowq) * 100
    return(perc)
}

# get the percentage of samples that have 0 expression in at least 100*q % of genes
get_perc_sample <- function(counts, q=.8){
    colq = colQuantiles(counts, probs=q)
    pos = which(colq == 0)
    perc = length(pos) / length(colq) * 100
    return(perc)
}

l = list()
l$pgene50 = get_perc_gene(m, q=.5)
l$pgene80 = get_perc_gene(m, q=.8)
l$psample50 = get_perc_sample(m, q=.5)
l$psample80 = get_perc_sample(m, q=.8)
df = as.data.frame(l)
write.table(df, output, sep=",", row.names=FALSE, quote=FALSE)

