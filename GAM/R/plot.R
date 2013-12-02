.node.color <- function(network, colors)
{ 
    colors <- colors[V(network)$name]
    colors2 <- colors
    # set red colors
    if(max(abs(colors))<5)
    {
        colors <- colors*5
    }
    if(any(colors>0))
    {
        max.red <- max(ceiling(abs(colors[which(colors>0)])))
        reds <- colorRampPalette(colors=c("white", "red"))
        red.vec <- reds(max.red+1)
        colors2[which(colors>0)] <- red.vec[ceiling(abs(colors[which(colors>0)]))+1]
    }
    # set green colors
    if(any(colors<0))
    {
        max.green <- max(ceiling(abs(colors[which(colors<0)])))
        greens <- colorRampPalette(colors=c("white", "green"))
        green.vec <- greens(max.green+1)
        colors2[which(colors<0)] <- green.vec[ceiling(abs(colors[which(colors<0)]))+1]
    }
    return(colors2)
}


#' Plot network 
#' @param module Module to plot
#' @param scale Scale factor for vertex sizes and label fonts
#' @param attr.label Attribute to use as a label
#' @param attr.shape Attribute to use for shape
#' @param layout Layout to use
#' @param ... Arguments for plot
#' @import igraph 
#' @examples
#' data(examplesGAM)
#' \dontrun{
#' plotNetwork(module.re)
#' }
#' @export
plotNetwork <- function(module, scale=1, attr.label="label", attr.shape="nodeType", layout=layout.kamada.kawai, ...) {
    network <- module
    
    if (is.null(V(network)$name)) {
        V(network)$name <- as.character(V(network))
    }
    
    if ("log2FC.norm" %in% list.vertex.attributes(network)) {
        de <- V(network)$log2FC.norm
        de <- de * 10
    } else {
        de <- V(network)$log2FC
    }
    names(de) <- V(network)$name
    de[is.na(de)] <- 0
    
    
    # Hack for coloring
    de[de < 0] <- de[de < 0] - 1
    de[de > 0] <- de[de > 0] + 1
    
    attr.labels = get.vertex.attribute(network, attr.label)
    names(attr.labels) <- V(network)$name
    attr.labels <- na.omit(attr.labels)
    
    node.labels = V(network)$name
    names(node.labels) <-node.labels
    all.labels = c(attr.labels, node.labels[!names(node.labels) %in% names(attr.labels)])
    
    labels = NULL
    scores = NULL
    main = NULL
    vertex.size = NULL
                         
    diff.expr <- de
    labels <- all.labels[node.labels]
    
    shapes <- rep("circle", length(V(network)))
    if (attr.shape %in% list.vertex.attributes(network)) {
        shapes.possible <- c("circle", "csquare")
        shapes <- shapes.possible[as.factor(get.vertex.attribute(network, attr.shape))]
    }
    
    names(shapes) <- V(network)$name
    
    if (!is.null(diff.expr) && !is.null(names(diff.expr))) {
        coloring <- .node.color(network, diff.expr)
    }
    else {
        coloring <- "SkyBlue2"
    }
    if (is.null(diff.expr) && "diff.expr" %in% list.vertex.attributes(network)) {
        diff.exprs = V(network)$diff.expr
        names(diff.exprs) <- V(network)$name
        coloring <- .node.color(network, diff.exprs)
    }
    max.labels <- max(nchar(labels))
    vertex.size2 <- 3    
    cex = 0.6
    network.size = length(V(network))
    layout <- layout(network)
    layout <- layout.norm(layout, -1, 1, -1, 1)
    layout.coordinates <- layout
    xs <- layout.coordinates[,1]
    ys <- layout.coordinates[,2]
    
    
    
    es <- get.edges(network, E(network))
    #es <- es + 1 # :ToDo: remove when moving from igraph0
    
    #message(paste("par:", par("pin")))
    scale <- scale * min(par("pin")) / 10
    
    if (nrow(es) > 0) {
        dist.sum <- 0
        for (i in 1:nrow(es)) {
            u <- es[i, 1]
            v <- es[i, 2]
            d <- sqrt((xs[u] - xs[v])^2 + (ys[u] - ys[v])^2)
            dist.sum <- dist.sum + d
        }
        
        
        dist.avg <- dist.sum / nrow(es)
        #message(paste("average distance:", dist.avg))
        scale <- scale * dist.avg * 10
    }
    
    #message(cex)
    
    vertex.size2 <- vertex.size2 * scale
    cex <- cex * scale

    random.state <- .Random.seed
    set.seed(42)
    plot(network, layout = layout.coordinates, vertex.size = vertex.size2, 
         vertex.label = labels, vertex.label.cex = cex, 
         vertex.label.dist = 0 , vertex.label.degree = -pi/2,
         vertex.color = coloring,
         vertex.label.family = "sans", 
         vertex.shape = shapes, 
         edge.label.cex = cex,
         main = main, ...)
    .Random.seed  <- random.state
}

