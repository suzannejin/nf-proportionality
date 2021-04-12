library(recount)
library(biomaRt)

# 
# this script divides the original count data (rse_gene data downloaded from recount)
# into subsets of data according to the biotype (protein_coding, lncRNA, processed_pseudogene, etc)
# 
### I used this script to divide the original data, then saved the output data in [cluster] nf-proportionality/data/counts/.

# arguments
args = commandArgs(trailingOnly=TRUE)
input = args[1]
outdir = args[2]
if (!dir.exists(outdir)) {dir.create(outdir, recursive=TRUE)}

# read count data
rse_gene = readRDS(input)   # rse_gene
m = assay(rse_gene)

# filter low-expression genes
# leave out the lower 25% of genes for a reasonable expression
av=apply(m,1,mean)
mygenes_av=rownames(m)[which(av>boxplot(av)$stats[2,1])] 
rse_gene = rse_gene[mygenes_av,] 
m = assay(rse_gene)

# rename gene names
mygenes = rownames(m)
mygenes = stringr::str_extract(mygenes, "[A-Z0-9]+")
rownames(m) = mygenes
rownames(rse_gene) = mygenes

# create biomart object
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
mart <- useMart(biomart="ensembl", 
                dataset="hsapiens_gene_ensembl", 
                host="http://aug2020.archive.ensembl.org")   # ensembl version 101

# filter according to ensemble gene id
ensembl = getBM(attributes = c( "ensembl_gene_id", "gene_biotype"), 
                 filters=c("ensembl_gene_id"),
                 values=list(mygenes), 
                 mart=mart)   # 42807 genes
# > dim(biotypes)
# [1] 42807     2
# > dim(m)
# [1] 43527  9662
# > length(setdiff(mygenes, biotypes$ensembl_gene_id))
# [1] 720
# NOTE that there are genes in recount2-GTEx that are no longer in the recent Ensembl version
biotypes = as.data.frame(table(ensembl$gene_biotype))


# divide count matrix into subsets according to biotype
for (biotype in biotypes$Var1){

    # get the genes of one biotype
    pos = which(ensembl$gene_biotype == biotype)
    genes = ensembl[pos,]$ensembl_gene_id
    rse_gene2 = rse_gene[genes,]

    # save subset data
    filename = paste(biotype, ".rds", sep="")
    saveRDS(rse_gene2, file=file.path(outdir, filename))
}
