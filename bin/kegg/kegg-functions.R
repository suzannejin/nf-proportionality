
# =============== #
# NODES AND EDGES #
# =============== #

# NODES

# get nodes from one kegg pathway
kegg2node <- function(pathway, genes=NULL){
    nod = graphite::nodes(pathway)
    nod = sapply(strsplit(nod, ":"), "[[", 2)  # convert SYMBOL:ID  to ID,  ENSEMBL:ID to ID, etc
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

# EDGES

# get edges from one kegg pathway
kegg2edge <- function(pathway, genes=NULL, full=FALSE, type=NULL, onlydirected=FALSE, avoidcyclic=FALSE){
    ed = graphite::edges(pathway)
    if (nrow(ed)>=1 & !is.null(genes)){  # if only interested in a specific set of genes
        ed = ed[which(ed$src %in% genes & ed$dest %in% genes),]
    }
    if (nrow(ed)>=1 & !is.null(type)){   # keep the interested interaction types
        ed = ed[which(ed$type %in% type),]
    }
    if (nrow(ed)>=1 & onlydirected){    # keep only the directed edges
        ed = ed[which(ed$direction=="directed"),]
    }
    if (full==FALSE){   # keep all columns if full, otherwise keep only src and dest columns
        ed = ed[,c("src","dest")]
    }
    if (avoidcyclic){   # avoid cyclic pairs
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

# get total edges in PathwayList
kegg2edges <- function(kegg, genes=NULL, full=FALSE, type=NULL, onlydirected=FALSE){
    edges = data.frame()
    for (k in 1:length(kegg)){
        pathway = kegg[[k]]
        ed = kegg2edge(pathway, genes, full=full, type=type, onlydirected=onlydirected)
        if (nrow(ed)>=1){
            edges = rbind(edges, ed)
        }
    }
    edges = unique(edges)
    return(edges)
}

# PATHWAY INFO

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

# get the kegg pathways that have a reasonable number of target genes
keggWithMygenes <- function(kegg, mygenes, min=2, max=NA){

    p = nnode_per_pathway(kegg, mygenes)
    p = p[p$n_nodes >= min,]
    if (!is.na(max)){
        p = p[p$n_nodes < max,]
    }
    kegg = kegg[p$pathway]

    return(kegg)
}


# ===================== #
# KEGG ADJACENCY MATRIX #
# ===================== #

# convert KEGG PathwayList object into adjacency matrix
kegg2matrix <- function(kegg, genes, clique=TRUE) {
    
    # initialize matrix
    ngenes = length(genes)
    mat = matrix(0, ncol=ngenes, nrow=ngenes)
    rownames(mat) = genes
    colnames(mat) = genes
    
    # for kegg pathway
    for (k in 1:length(kegg)){
        
        # update matrix values
        # if clique, then assign an edge between each pair of nodes in the same pathway
        if (clique){
            nodes = kegg2node(kegg[[k]], genes)  # get the nodes in pathway 
            if (length(nodes)>=2){
                for (i in nodes){
                    for (j in nodes){
                        if (i != j){  # avoid diagonal
                            mat[i,j] = 1
                        }
                    }
                }
            }
        }else{
            # if not clique, then only use the edges annotated in kegg database
            edges = kegg2edge(kegg[[k]], genes, type=NULL, onlydirected=FALSE, avoidcyclic=FALSE)   # get annotated edges
            if (nrow(edges)>=2){
                for(row in 1:nrow(edges)){
                    edge = edges[row,]
                    i = edge$src
                    j = edge$dest
                    mat[i,j] = 1
                    mat[j,i] = 1
                }
            }
        }
    }
    return(mat)
}

# convert KEGG edges data frame into adjacency matrix
edges2matrix <- function(edges, mygenes){
    
    # initialize matrix
    ngenes = length(mygenes)
    mat = matrix(0, ncol=ngenes, nrow=ngenes)
    colnames(mat) = mygenes
    rownames(mat) = mygenes
    
    # update value for edge
    for (e in 1:nrow(edges)){
        i = edges[e,"src"]
        j = edges[e, "dest"]
        mat[i,j] = 1
        mat[j,i] = 1
    }
    
    return(mat)
}
