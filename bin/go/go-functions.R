
# =============== #
# REFERENCE GRAPH #
# =============== #

getK2 <- function (keys, columns = "GO", keytype = "ENSEMBL", minK = 10, maxK=NA, ont=NA, evi=NA) 
{
    if (length(columns) > 1) 
        stop("Please provide a single column.")
    if (!is.na(ont)){
        columns = c(columns, "ONTOLOGY")
    }
    if (!is.na(evi)){
        columns = c(columns, "EVIDENCE")
    }
    # check packages
    library("AnnotationDbi")
    library("org.Hs.eg.db")
    db <- org.Hs.eg.db::org.Hs.eg.db
    
    # get pathways-gene data frame
    godf <- AnnotationDbi::select(db, keys = keys, columns = columns, keytype = keytype)
    
    # filter if required
    if (!is.na(ont)){
        godf = godf[which(godf$ONTOLOGY %in% ont),]
    }
    if (!is.na(evi)){
        godf = godf[which(godf$EVIDENCE %in% evi),]
    }
    # reorganize data frame
    gotab <- table(godf[, 1:2])
    gotab[gotab > 1] <- 1
#     gotab <- gotab[keys, ]
    if (length(rownames(gotab)) > length(which(rownames(gotab) %in% keys))) {
        stop("Uh oh! Unexpected mapping.")
    }
    gotab = gotab[, colSums(gotab) >= minK]   # only keep the pathways with a minimun number of nodes 
    if (!is.na(maxK)){
        gotab = gotab[, colSums(gotab) < maxK]
    }
    return(gotab)
}

# convert gene-pathway matrix into an adjacency matrix
K2adjmatrix <- function(K, genes) {
    
    # initialize matrix
    ngenes = length(genes)
    mat = matrix(0, ncol=ngenes, nrow=ngenes)
    rownames(mat) = genes
    colnames(mat) = genes
    
    # for kegg pathway
    for (k in 1:ncol(K)){
        nodes = names(which(K[,k]==1))
        if (length(nodes)>=2){
            for (i in nodes){
                for (j in nodes){
                    if (i != j){  # avoid diagonal
                        mat[i,j] = 1
                    }
                }
            }
        }
    }
    return(mat)
}