rm.comp <- function(s) { gsub("(\\[[a-z]\\]|\\([a-z]\\))$", "" , s) }


unique.network <- function(net) {
    net$edges <- unique(net$edges)
    net$meta <- net$meta[match(unique(net$meta$v), net$meta$v),]
    return(net)
}

rm.comp.network <- function(net) {
    rm.comp1 <- function(u) {
        u.withComp <- net$meta[match(u, net$meta$v), "nodeType"] == "met"
        u.withComp[is.na(u.withComp)] = F
        u[u.withComp] <- sapply(u[u.withComp], rm.comp)
        return(u)        
    }
    
    net$edges$u <- rm.comp1(net$edges$u)
    net$edges$v <- rm.comp1(net$edges$v)    
    net$meta$v <- rm.comp1(net$meta$v)
    net$meta$nodeComp <- NULL
    return(unique.network(net))
}

rm.nodes <- function(net, to_remove) {
    to_keep_u <- is.na(match(net$edges$u, to_remove))
    to_keep_v <- is.na(match(net$edges$v, to_remove))
    net$edges <- net$edges[to_keep_u & to_keep_v,]
    to_keep_meta <- is.na(match(net$meta$v, to_remove))
    net$meta <- net$meta[to_keep_meta,]
    return(net)
}

nodes.network <- function(net) {
    return(unique(c(net$edges$u, net$edges$v)))
}

rm.ex.network <- function(net) {
    nodes <- nodes.network(net)
    to_remove <- nodes[grep("Ex$", nodes)]
    return(rm.nodes(net, to_remove))
}

convert.ids.network <- function(net, from, to) {
    convert.ids1 <- function(u) {
        u.match <- match(u, from)        
        u[!is.na(u.match)] <- to[u.match[!is.na(u.match)]]
        return(u)
    }
    net$edges$u <- convert.ids1(net$edges$u)
    net$edges$v <- convert.ids1(net$edges$v)
    net$meta$v <- convert.ids1(net$meta$v)
    return(unique.network(net))
}

add.square.network <- function(net) {    
    edges <- net$edges
    
    edges_rev <- data.frame(u=edges$v, e=edges$e, v=edges$u, stringsAsFactors=F)
    edges <- unique(rbind(edges, edges_rev))
    
    edges2 <- edges
    colnames(edges2) <- c("v", "e", "w")
    edges_sq <- merge(edges, edges2, by="v")
    edges_sq <- data.frame(u=edges_sq$u, e="act", v=edges_sq$w, stringsAsFactors=F)
    edges_sq <- unique(rbind(edges, edges_sq))
    net$edges <- edges_sq
    return(net)
}

escape.names <- function(net) {
    from=nodes.network(net)
    to=gsub("'", "_prime_", from, fixed=T)
    return(convert.ids.network(net, from, to))
}

edgelist <- function(network) {
    edges <- edges(network)
    vs <- unlist(edges)
    us <- rep(names(edges), sapply(edges, length))
    res <- data.frame(u=us, v=vs, stringsAsFactors=F)
    return(res)
}

convert.node.names <- function(network, from, to) {
    edges <- edgelist(network)
    m1 <- match(edges$u, from)
    edges$u[!is.na(m1)] <- to[m1[!is.na(m1)]]
    m2 <- match(edges$v, from)
    edges$v[!is.na(m2)] <- to[m2[!is.na(m2)]]
    edges <- unique(edges)
    network2 <- graph.edgelist(as.matrix(edges), directed=F)
    network2 <- simplify(network2, remove.multiple=T)
    network2 <- igraph.to.graphNEL(network2)
    
    nodeDataDefaults(network2) <- nodeDataDefaults(network)
    
    old_nodes <- intersect(nodes(network), nodes(network2))
    
    for (attr in names(nodeDataDefaults(network))) {
        nodeData(network2, old_nodes, attr) <- nodeData(network, old_nodes, attr)
    }
    
    network2   
}
