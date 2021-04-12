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
    make_option("--output",
                type="character",
                default=NULL,
                help="Output kegg pathways into .rds file",
                metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list, usage="Organize KEGG pathways. This script keeps the pathways that have at least 2 genes of our interest")
opt = parse_args(opt_parser)

# ======================= #
# LIBRARIES AND FUNCTIONS #
# ======================= #

library(graphite)

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

# get the kegg pathways that have at least 2 genes from my set
keggWithMygenes <- function(kegg, mygenes){

    p = nnode_per_pathway(kegg, mygenes)
    p = p[p$n_nodes >= 2,]
    kegg = kegg[p$pathway]

    return(kegg)
}

# get nodes from one kegg pathway
kegg2node <- function(pathway, genes=NULL){
    nod = graphite::nodes(pathway)
    nod = sapply(strsplit(nod, ":"), "[[", 2)  # convert SYMBOL:ID  to ID 
    if (length(nod)>=1 & !is.null(genes)){
        nod = nod[which(nod %in% genes)]
    }
    return(nod)
}

# get total genes in PathwayList
kegg2nodes <- function(kegg, genes=NULL){
    nodes = c()
    for (k in 1:length(kegg)){
        pathway = kegg[[k]]
        nod = kegg2node(pathway, genes)
        if (length(nod)>=1){
            nodes = c(nodes, nod)
        }
    }
    nodes = unique(nodes)
    return(nodes)
}

# get number of nodes and edges per pathway
nnode_per_pathway <- function(kegg, genes=NULL){
    nodes = c(); edges = c()
    for (k in 1:length(kegg)){
        pathway = kegg[[k]]
        
        # get nodes and edges
        nod = kegg2node(pathway, genes)
        ed = kegg2edge(pathway, genes)
        
        # update counts
        nodes = c(nodes, length(nod))
        edges = c(edges, nrow(ed))
    }
    df = data.frame(pathway=names(kegg), n_nodes=nodes, n_edges=edges)
    return(df)
}

# get edges from one kegg pathway
kegg2edge <- function(pathway, genes=NULL, full=FALSE, type=NULL, onlydirected=FALSE, avoidcyclic=FALSE){
    ed = graphite::edges(pathway)
    if (nrow(ed)>=1 & !is.null(genes)){
        ed = ed[which(ed$src %in% genes & ed$dest %in% genes),]
    }
    if (nrow(ed)>=1 & !is.null(type)){
        ed = ed[which(ed$type %in% type),]
    }
    if (nrow(ed)>=1 & onlydirected){
        ed = ed[which(ed$direction=="directed"),]
    }
    if (full==FALSE){
        ed = ed[,c("src","dest")]
    }
    if (avoidcyclic){
        rows = rownames(ed)
        for (row in rows){
            e = ed[rownames(ed) == row,]
            if (nrow(e)<1){print(row)}
            # if duplicated in different direction
            if (e$src %in% ed$dest && e$dest %in% ed$src){
                ed = ed[rownames(ed) != row,]
            }
        }
    }
    return(ed)
}


# ============ #
# PROCESS DATA #
# ============ #

# read input propr object 
print_msg("Read and process propr object")
pro = readRDS(opt$pro)
matrix = pro@matrix
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
saveRDS(kegg, opt$output)
