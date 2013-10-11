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

#' Get data.frame with network edges
#' @param network Network
#' @return Data.frame with edges
#' @export
edgelist <- function(network) {
    edges <- edges(network)
    vs <- unlist(edges)
    us <- rep(names(edges), sapply(edges, length))
    res <- data.frame(u=us, v=vs, stringsAsFactors=F)
    return(res)
}

#' @export
splitMappingByConnectivity <- function(connections, from, to) {
    connections <- data.frame(connections)
    colors <- seq_along(from)
    names(colors) <- from
    names(to) <- from
    connections <- connections[connections[,1] != connections[,2],]
    c1 <- connections[to[connections[,1]] == to[connections[,2]],]
    while (T) {
        c2 <- data.frame(x=colors[c1[,1]], y=colors[c1[,2]])
        
        c2[c2$x > c2$y, ] <- c2[c2$x > c2$y, c("y", "x")]
        c1 <- c1[c2$x != c2$y, ]
        c2 <- c2[c2$x != c2$y, ]
        if (nrow(c2) == 0) {
            break;
        }
        c2 <- aggregate(x ~ y, c2, min)
        
        matched <- colors %in% c2$y
        colors[matched] <- c2$x[match(colors[matched], c2$y)]
    }
#     while (T) {
#         c2 <- data.table(x=colors[c1[,1]], y=colors[c1[,2]])
#         
#         c2[x > y, ] <- c2[x > y, list(y,x)]
#         c2 <- c2[x != y, ]
#         if (nrow(c2) == 0) {
#             break;
#         }
#         c2 <- c2[, list(x=min(x)), by=y]
#         colors[match(c2$y, colors)] <- c2$x
#     }
    
    
    res <- paste(to, colors, sep=".")
    names(res) <- from
    res
}

# splitMappingByConnectivity <- function(connections, from, to) {
#     graph <- graphNEL.from.tables(edge.table=connections, directed=F)
#     names(to) <- from
#     
#     es <- edges(graph)
#     to.new <- to
#     to.new[] <- NA
#     
#     comp_names <- unique(to)
#     comp_counts <- rep(0, length(comp_names))
#     names(comp_counts) <- comp_names
#     
#     comp_sizes <- rep(0, length(comp_names))
#     names(comp_sizes) <- comp_names
#     
#     dfs <- function(node, mark.allow, new.name) {    
#         if (is.na(to[node]) || (to[node] != mark.allow)) {
#             return(0)
#         }
#         
#         if (!is.na(to.new[node])) {
#             return(0)
#         }
#         to.new[node] <<- new.name
#         
#         res <- 1
#         for (v in es[[node]]) {
#             res <- res + dfs(v, mark.allow, new.name)
#         }
#         res
#     }
#     
#     for (node in from) {
#         if (!is.na(to.new[node])) {
#             next
#         }
#         comp_n <- comp_counts[to[node]] + 1
#         comp_counts[to[node]] <- comp_n
#         comp_name <- paste(to[node], comp_n, sep=".")
#         comp_size <- dfs(node, to[node], comp_name)    
#         comp_sizes[comp_name] <- comp_size
#     }
#     
#     to.new
# }


#" @param graph igraph object
addNodeAttributes <- function(graph, node.table=list(), node.col=1, name.as.label=T) {
    if (is.list(node.table)) {
        for (name in names(node.table)) {
            node.table[[name]]$nodeType <- name
            node.table[[name]] <- moveColumnsToFront(node.table[[name]], node.col)
        }
        node.table <- do.call(rbind.fill, node.table)
    } else if (!is.null(node.table)) {
        node.table <- moveColumnsToFront(node.table, node.col)
    }
    
    if ("name" %in% names(node.table) && name.as.label) {
        node.table <- rename(node.table, c("name"="label"))
    }
    
    node.table <- node.table[node.table[,1] %in% V(graph)$name, ]
    
    res <- graph
    for (name in names(node.table)[-1]) {
        res <- set.vertex.attribute(res, name, node.table[,1], node.table[, name])
    }
    
    return(res)
}

moveColumnsToFront <- function(d, cols) {
    if (is.character(cols)) {
        cols <- match(cols, colnames(d))
    }
    d[, c(cols, seq_along(d)[-cols])]
}

#' @importFrom igraph graph.edgelist igraph.to.graphNEL simplify graph.data.frame
#' @importFrom plyr rbind.fill
#' @export
graphNEL.from.tables <- function(node.table=NULL, edge.table,
                                 node.col=1, edge.cols=c(1,2),
                                 directed=T,
                                 name.as.label=T) {    
    edge.table <- moveColumnsToFront(edge.table, edge.cols)
    
    if ("name" %in% names(edge.table) && name.as.label) {
        edge.table <- rename(edge.table, c("name"="label"))
    }
    
    net1 <- graph.data.frame(edge.table, directed=directed)
    net1 <- delete.vertices(net1, V(net1)[degree(net1) == 0])
    net1 <- addNodeAttributes(net1, node.table, node.col, name.as.label)
    net1 <- igraph.to.graphNEL(net1)
    return(net1)
}
