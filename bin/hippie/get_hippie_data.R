
library(AnnotationDbi)
library(org.Hs.eg.db)

# arguments
args = commandArgs(trailingOnly=TRUE)
output = args[1]  # output hippie data file 

# read hippie data
src = url("http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt")
hippie = read.table(src, sep="\t")
colnames(hippie) = c("symbol1", "entrez1", "symbol2", "entrez2", "score", "info")

# entrez id to ensembl id
hippie$ensembl1 = AnnotationDbi::mapIds(org.Hs.eg.db,
                          keys=as.character(hippie$entrez1), 
                          column="ENSEMBL",
                          keytype="ENTREZID")
hippie$ensembl2 = AnnotationDbi::mapIds(org.Hs.eg.db,
                          keys=as.character(hippie$entrez2), 
                          column="ENSEMBL",
                          keytype="ENTREZID")

# reorder columns
hippie = hippie[,c("symbol1","symbol2","entrez1","entrez2","ensembl1","ensembl2","score","info")]

# save output
write.table(hippie, file=output, quote=FALSE, sep="\t", row.names=FALSE)
