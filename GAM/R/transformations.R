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


#' Add attributes for nodes from data.frame
#' @param graph igraph object
#' @param node.table Table with attribute values
#' @param node.col Column (number of character) with vertex IDs
#' @param name.as.label If TRUE rename "name" attribute to "label"
#' @import igraph
addNodeAttributes <- function(graph, node.table, node.col=1, name.as.label=TRUE) {
    if (class(node.table) == "list") {
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

#' @import igraph
#' @importFrom plyr rbind.fill
graph.from.tables <- function(node.table=NULL, edge.table,
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
    return(net1)
}
