#!/usr/local/bin/Rscript

# ========= #
# ARGUMENTS #
# ========= #

library(optparse)

option_list = list(
  make_option("--data", 
              type="character", 
              default=NULL, 
              help="Processed expression count matrix: rows=samples, columns=genes in .rds format", 
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
              metavar="number"),
  make_option("--ncores",
              type="integer",
              default=1,
              help="Number of cores",
              metavar="number"),
  make_option("--chunk_size",
              type="integer",
              default=1000,
              help="Divide the data into chunks of n gene features and compute them separately to decrease the RAM usage",
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

# This function divides the counts data into chunks of size up to 2n gene features
#   Here n is the size of one gene group
#   The chunk is formed by two gene groups, so that the genes of each pair of groups will be compared pairwisely
#   NOTE: since one gene group might be used many times (group1 to group2, group1 to group3, etc),
#         the computation takes longer time than the original code which uses all gene features at once.
# Eg. if we have 100 gene features and we want chunks of size 30, 
# we will get 4 groups of 30 gene features each, except the last group that will have 10 gene features  
blocks2combs <- function(counts, n){

  # define blocks
  nblocks = ncol(counts) %/% n
  ngroup <- rep(1:nblocks, each = n)
  leftover <- ncol(counts) - length(ngroup)
  if(leftover > 0) ngroup <- c(ngroup, rep(nblocks + 1, leftover))

  # check size
  # here I decided to stop computation if no more than 2 groups are generated (so only 1 chunk)
  # so if this happens you better change n
  # alternatively you can change this code to run normal propr function without dividing the data into chunks
  if (length(unique(ngroup)) <= 2){
    stop(paste("ERROR: chunk size ", n, " is too big for data frame of size [", nrow(counts), "][", ncol(counts), "]", sep=""))
  }

  # split groups
  # each row in combs define a chunk
  split <- split(1:ncol(counts), ngroup)
  combs <- expand.grid(1:length(split), 1:length(split))
  combs <- t(apply(combs, 1, sort))
  combs <- unique(combs) 
  combs <- combs[combs[,1] != combs[,2],]

  print_msg("Running propr for ", nrow(combs), " chunks of size ", n)

  return(list(combs, split))
}

chunk2propr <- function(i, chunk, metric="rho", ivar=NA, symmetrize=FALSE,
                        alpha, p=100, interval=seq(0.4,0.95,0.01), fdr=0.05 ){

  # compute propr for chunk
  rho.i <- propr(chunk, metric = metric, ivar = ivar, alpha = alpha, p=p)
  rho.i <- updateCutoffs(rho.i,interval)

  # get cutoff | fdr
  df <- data.frame('cutoff'=rho.i@fdr[,'cutoff'], 'FDR'=rho.i@fdr[,'FDR'])

  return(list(rho.i@matrix, df))
}

chunk2full <- function(counts, RES, split, combs, fdr){
  
  # define variables
  QUILT <- matrix(0, ncol(counts), ncol(counts))
  d <- data.frame('cutoff'=RES[[1]][[2]][,'cutoff'])

  # collect chunks
  for(i in 1:nrow(combs)){

    # add fdr.i
    d[paste('FDR', i, sep="")] <- RES[[i]][[2]][,'FDR']

    # Fill final matrix with each chunk
    batch1 <- split[[combs[i,1]]]
    batch2 <- split[[combs[i,2]]]
    patch.i <- c(batch1, batch2)
    QUILT[patch.i, patch.i] <- RES[[i]][[1]]
  }

  # average fdr
  df <- data.frame('cutoff'=RES[[1]][[2]][,'cutoff'], 'FDR'=round(rowSums(d[,2:ncol(d)])/i,4))
  # cutoff
  cutoff <- min(df[df[,"FDR"]<fdr,"cutoff"])

  # rename columns & rows
  matrix <- QUILT
  rownames(matrix) <- colnames(counts)
  colnames(matrix) <- colnames(counts)

  return(list(matrix, cutoff, df))
}

file2res <- function(RES, combs, dir){

  print_msg("Collecting ", nrow(combs), " chunk files previously computed")

  for(i in 1:nrow(combs)){
    file2 <- paste0(dir, "/job-", combs[i,1], "+", combs[i,2], ".csv")
    csv <- read.csv(file2, row.names = 1, header= TRUE)
    matrix <- data.matrix(csv)
    RES[[i]][[1]] <- matrix
  }

  # remove tmp files
  unlink(dir, recursive=TRUE)

  return(RES)
}

propr.chunk <- function(counts, metric = c("rho", "phi", "phs", "cor", "vlr"),
                        ivar = NA, symmetrize = FALSE, alpha=NA, p=100, fdr=0.05,
                        n=100, ncores = 1, interval=seq(0.3,0.95,0.01),
                        dir = NA){
  
  # divide data into chunks
  # I recommend using large chunk size n,
  # otherwise dividing the data into too many chunks will slow down the computation as well
  l <- blocks2combs(counts, n)
  combs <- l[[1]]
  split <- l[[2]]
  
  # parallelize propr computation for i chunks
  doMC::registerDoMC(cores = ncores)
  `%dopar%` <- foreach::`%dopar%`
  
  RES <- foreach::foreach(i = 1:nrow(combs)) %dopar% {

    # get chunk
    # the chunk will have size [nsample, ngenes] : 
    #   nsample = number of samples in counts 
    #   ngenes = number of gene features defined in the given chunk (maximun 2n gene features)
    batch1 <- split[[combs[i,1]]]
    batch2 <- split[[combs[i,2]]]
    chunk = subset(counts, select = c(batch1, batch2))

    # compute propr 
    l_propr <- chunk2propr(i, chunk, metric=metric, ivar=ivar, symmetrize=symmetrize, alpha=alpha, p=p, interval=interval, fdr=fdr)

    if(is.na(dir)){
      # return propr matrix, cutoff, and FDR
      l_propr
    }else{
      # save data if required
      file2 <- paste0(dir, "/job-", combs[i,1], "+", combs[i,2], ".csv")
      print_msg("Saving tmp file[", i, "] to ", file2)
      write.csv(l_propr[[1]], file=file2)

      # return cutoff and FDR
      list(NULL, l_propr[[2]])
    }
  }

  # collect files
  if(!is.na(dir)){
    RES <- file2res(RES, combs, dir)
  }

  # merge chunks
  l_propr <- chunk2full(counts, RES, split, combs, fdr)

  return(l_propr)
}


# ===================== #
# COMPUTE PROPR/PEARSON #
# ===================== #

# load count matrix
print_msg("Load count matrix")
M = readRDS(file.path(opt$data)) 

# output folder
if (!dir.exists(opt$outdir)) {dir.create(opt$outdir, recursive=TRUE)}
dir2 = paste(opt$outdir, "tmp", sep="/")
if (!dir.exists(dir2)) {dir.create(dir2, recursive=TRUE)}

# calculate association coefficients
print_msg("Calculate coefficients [", opt$method, "]")
ch = propr.chunk(M, metric=opt$method, ivar=NA, alpha=NA, p=opt$permutation, 
                 interval=seq(opt$interval_min, opt$interval_max, opt$cutoff_interval), 
                 ncores=opt$ncores,  n=opt$chunk_size, dir=dir2)
matrix = ch[[1]]
cutoff = ch[[2]]
fdr = ch[[3]]

# ============ #
# WRITE OUTPUT #
# ============ #

print_msg("Saving output files")

# write matrix
matout = file.path(opt$outdir, "propr_results.mat")
write.table(matrix, matout)

# write cutoff
cutout = file.path(opt$outdir, "propr_results_cutoff.txt")
writeLines(as.character(cutoff), cutout)

# write FDR
fdrout = file.path(opt$outdir, "propr_results_fdr.txt")
write.table(fdr, fdrout, row.names = FALSE, quote=FALSE, sep=",", dec=".")
