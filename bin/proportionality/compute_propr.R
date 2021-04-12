#!/usr/local/bin/Rscript

# ========= #
# ARGUMENTS #
# ========= #

library(optparse)

option_list = list(
  make_option("--data", 
              type="character", 
              default=NULL, 
              help="Preprocessed expression count matrix: rows=samples, columns=genes in .rds format", 
              metavar="character"),
  make_option("--method",
              type="character",
              default="rho",
              help="Method for the computation of association coefficients. Default = %default"),
  make_option("--outdir",
              type="character",
              default=getwd(),
              help="Output directory",
              metavar="character"),
  make_option("--cutoff_interval",
              type="double",
              default=0.005,
              help="Interval for FDR computation (updateCutoffs step). Default = %default"),
  make_option("--interval_min",
              type="double",
              default=0.1,
              help="Minimum cutoff for FDR computation. Default = %default"),
  make_option("--interval_max",
              type="double",
              default=0.995,
              help="Maximum cutoff for FDR computation. Default = %default"),
  make_option("--permutation",
              type="integer",
              default=20,
              help="Permutation number for FDR computation. Default = %default",
              metavar="number")
); 
 
opt_parser = OptionParser(option_list=option_list, usage="Preprocess data for proportionality computation")
opt = parse_args(opt_parser)

# ========= #
# LIBRARIES #
# ========= #

library(propr)

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

# ===================== #
# COMPUTE PROPR/PEARSON #
# ===================== #


# load count matrix
print_msg("Load count matrix")
M = readRDS(file.path(opt$data))  

# calculate association coefficients
print_msg("Calculate coefficients [", opt$method, "]")
pro=propr(M, metric=opt$method, ivar=NA, p=opt$permutation)  # use p=20 for faster computation. Moreover, p is only useful for the calculation of FDR, and now we don't care that much on the FDR.
pro=updateCutoffs(pro,seq(opt$interval_min, opt$interval_max, opt$cutoff_interval))

# ============ #
# WRITE OUTPUT #
# ============ #

print_msg("Saving output files")

# create output directory
if (!dir.exists(opt$outdir)) {dir.create(opt$outdir, recursive=TRUE)}

# write Rdata
rdata = file.path(opt$outdir, "propr_results.rds")
saveRDS(pro, file=rdata)

# # write matrix
# matout = file.path(opt$outdir, "propr_results.mat")
# write.table(pro@matrix, matout)

# write cutoff
cutoff = min(pro@fdr[which(pro@fdr[,"FDR"]<0.05), "cutoff"])
cutout = file.path(opt$outdir, "propr_results_cutoff.txt")
writeLines(as.character(cutoff), cutout)

# write FDR
fdrout = file.path(opt$outdir, "propr_results_fdr.txt")
write.table(pro@fdr, fdrout, row.names = FALSE, quote=FALSE, sep=",", dec=".")
