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
    make_option("--hippie", 
                type="character", 
                default=NULL, 
                help="Input hippie tsv file", 
                metavar="character"),
    make_option("--outdir",
                type="character",
                default=getwd(),
                help="Output directory",
                metavar="character"),
    make_option("--confidence",
                type="integer",
                default=0,
                help="HIPPIE score confidence level. 1=medium level (score 0.63 - second quantile of the full HIPPIE distribution). 2=high level (score 0.73 - third quantile). 3=filter out the lower quantile of the remaining HIPPIE pairs containing target genes")
); 
 
opt_parser = OptionParser(option_list=option_list, usage="Coexpression network vs HIPPIE protein-protein interaction: ROC and PR curves")
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
source(file.path(bindir, "utils/compute_roc.R"))
source(file.path(bindir, "kegg/kegg-functions.R"))


# ====================== #
# LOAD AND ORGANIZE DATA
# ====================== #

# read input propr object 
print_msg("Read and process propr object")
pro = readRDS(opt$pro)
matrix = pro@matrix
mygenes = rownames(matrix)
rm(pro)

# read hippie data
print_msg("Read and process hippie PPI data")
hippie = read.csv(opt$hippie, sep="\t")
# remove self-interacting pairs
hippie = hippie[which(hippie$ensembl1 != hippie$ensembl2),]
# match our genes with the database
hippie = hippie[which(hippie$ensembl1 %in% mygenes & hippie$ensembl2 %in% mygenes),]
# filter interactions based on the confidence level
if (opt$confidence >=1 ){
    if (opt$confidence == 1){
        filt = 0.63
    }else if (opt$confidence == 2){
        filt = 0.73
    }else if (opt$confidence == 3){
        filt = as.numeric(quantile(hippie$score, prob=.25))  # 0.49
    }
    print_msg("Filtering out the PPIs with a score lower than ", filt)
    hippie = hippie[which(hippie$score >= filt),]
}
print_msg("Obtained a final number of ", nrow(hippie), " protein-protein interactions")

# reference hippie graph
print_msg("Extract the reference hippie graph")
K = edges2matrix(hippie, mygenes, col1="ensembl1", col2="ensembl2")


# ============ #
# ROC ANALYSIS #
# ============ #

# scores and labels
print_msg("Organize scores and labels")
scores = matrix[lower.tri(matrix)]
labels = K[lower.tri(K)]

print_msg("Compute ROC and PR curves for ", length(labels), " items")
compute_roc(scores, labels, opt$outdir)
