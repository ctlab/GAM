#' @import BioNet 
NULL

#' Load data only if its absent from global environment
lazyData <- function(name, ...) {
    if (!name %in% ls(envir=.GlobalEnv)) {
        print(paste0("No ", name, ", loading"))
        data(name, ...)
    }
}

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


#' Plots network module
#' @param module Module to plot
#' @param scale Scale factor for vertex sizes and label fonts
#' @param attr.label Attribute to use as a label
#' @param attr.shape Attribute to use for shape
#' @param layout Layout to use
#' @param ... Arguments for plot
#' @importFrom igraph E V igraph.from.graphNEL list.vertex.attributes layout.kamada.kawai layout.norm get.edges plot.igraph get.vertex.attribute
#' @export
plotNetwork <- function(module, scale=1, attr.label="label", attr.shape="nodeType", layout=layout.kamada.kawai, ...) {
    if (!(attr.label %in% names(nodeDataDefaults(module)))) {
        return()        
    }
    # Hack for coloring
    
    if ("logFC.norm" %in% names(nodeDataDefaults(module))) {
        de <- unlist(nodeData(module, attr="logFC.norm"))
        de <- de * 10
    } else {
        de <- unlist(nodeData(module, attr="logFC"))    
    }
    de[is.na(de)] <- 0
    
    de[de < 0] <- de[de < 0] - 1
    de[de > 0] <- de[de > 0] + 1
    
    attr.labels = na.omit(unlist(nodeData(module, nodes(module), attr=attr.label)))
    node.labels = nodes(module)
    names(node.labels) <-node.labels
    all.labels = c(attr.labels, node.labels[!names(node.labels) %in% names(attr.labels)])
    
    labels = NULL
    scores = NULL
    main = NULL
    vertex.size = NULL
                         
    diff.expr <- de
    labels <- all.labels[node.labels]
    network <- module
    
    if (is(network, "graphNEL")) {
        network <- igraph.from.graphNEL(network)
    }
    if (is.null(V(network)$name)) {
        V(network)$name <- as.character(V(network))
    }
    if (is.null(labels)) {
        if ("geneSymbol" %in% list.vertex.attributes(network)) {
            labels <- V(network)$geneSymbol
        }
        else {
            labels <- V(network)$name
        }
    }
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
    
    #print(paste("par:", par("pin")))
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
        #print(paste("average distance:", dist.avg))
        scale <- scale * dist.avg * 10
    }
    
    #print(cex)
    
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

#' Add node attribute with normalized log-foldchange inside node groups
#' @param module Module to add attribute to
#' @param logFC.attr Attribute with node log-foldchange
#' @param logFC.norm.attr Name of a new attribute
#' @param group.by Attribute by which node grouping shoud happen
#' @return Modified module with normalized log-foldchange node attribute
#' @export
addNormLogFC <- function(module, logFC.attr="logFC", logFC.norm.attr="logFC.norm", group.by="nodeType") {
    nodeDataDefaults(module, logFC.norm.attr) <- NA
    module.nodes <- nodes(module)
    node.types <- unique(na.omit(unlist(nodeData(module, module.nodes, group.by))))
        
    for (node.type in node.types) {
        module.type.nodes <- module.nodes[which(unlist(nodeData(module, module.nodes, group.by)) == node.type)]
        module.type.de <- unlist(nodeData(module, module.type.nodes, logFC.attr))        
        module.type.norm.de <- module.type.de / max(abs(na.omit(module.type.de)))
        nodeData(module, module.type.nodes, logFC.norm.attr) <- unname(module.type.norm.de)
    }
    
    module
}

#' Scores p-value by fitted BUM-model
#' This is a helper function based on BioNet::scoreFunction()
scoreValue <- function (fb, pval, fdr = 0.01) 
{
    return((fb$a - 1) * (log(pval) - log(fdrThreshold(fdr, fb))))
}


#' Preprocess experiment set's differential expression data for genes and metabolites
#' Convert metabolite and gene IDs to be the same as in network. Computes
#' reaction differential expression data.
#' @param es Experiment set with DE data
#' @param met.ids Type of IDs used in metabolite DE data (@see met.id.map for possible values)
#' @param gene.ids Type of IDs used in gene DE data (@see gene.id.map for possible values)
# :ToDo: this function is ugly, refactor
preprocessPvalAndMetDE <- function(es, met.ids, gene.ids, plot=T) {
    if (!is.null(es$gene.de)) {
        print("Processing gene p-values...")
        es$gene.de$logFC <- fixInf(es$gene.de$logFC)
        if (!is.null(gene.ids)) {
            lazyData("gene.id.map")
            es$gene.de <- convertPval(es$gene.de, 
                                       from=gene.id.map[,gene.ids], 
                                       to=gene.id.map[,es$network$gene.ids])
        }    
        es$gene.de$origin <- NULL
        print("Converting gene p-values to reactions...")
        es$rxn.de <- convertPval(es$gene.de, from=es$network$rxn2gene$gene, to=es$network$rxn2gene$rxn)
        
    }
    
    if (!is.null(es$met.de)) {
        print("Processing metabolite p-values...")
        es$met.de$logFC <- fixInf(es$met.de$logFC)        
        if (!is.null(met.ids)) {
            lazyData("met.id.map")
            es$met.de <- convertPval(es$met.de, 
                                      from=met.id.map[,met.ids], 
                                      to=met.id.map[,es$network$met.ids])
        }
        
        es$met.pval <- es$met.de$pval
        names(es$met.pval) <- es$met.de$ID
        es$fb.met <- fitBumModel(es$met.pval, plot=F)
        if (plot) {
            hist(es$fb.met, main="Histogram of metabolite p-values")
            plot(es$fb.met, main="QQ-Plot for metabolite BUM-model")
        }
    } else {
        es$met.pval <- NULL
    }
    
    return(es)
}

#' Make preprocessing necessary for findModule function
#' Converts gene and metabolite ids, makes an actual to work with.
#' @param network Network object
#' @param met.de Differential expression data for metabolites
#' @param gene.de Differential expression data for genes
#' @param rxn.de Differential expression data for reactions
#' @param met.ids Type of IDs used in met.de
#' @param gene.ids Type of IDs used in gene.de
#' @param reactions.as.edges If TRUE represent reaction as edges betwen metabolites,
#'                          otherwise represent them as nodes with connections to
#'                          compounds
#' @param collapse.reactions If TRUE collapse reaction nodes if they share enzyme
#'                           and at least one metabolite
#' @param use.rpairs If TRUE only rpairs will be used as reaction edges
#' @param plot If TRUE plot BUM-models
#' @export
#' @importFrom plyr rename
#' @import data.table
makeExperimentSet <- function(network, 
                              met.de=NULL, gene.de=NULL, rxn.de=NULL,
                              met.ids=NULL, gene.ids=NULL,
                              reactions.as.edges=T,
                              collapse.reactions=T,
                              use.rpairs=T,
                              plot=T) {
    es <- newEmptyObject()
    es$network <- network
    es$met.de <- met.de
    es$gene.de <- gene.de
    es$rxn.de <- rxn.de 
    
    
    es$graph.raw <- es$network$graph.raw
    es$reactions.as.edges <- reactions.as.edges
    
    
    es <- preprocessPvalAndMetDE(es, met.ids=met.ids, gene.ids=gene.ids, plot=plot)
    
    if (!is.null(es$rxn.de)) {
        print("Processing reaction p-values...")
        es$rxn.de$logFC <- fixInf(es$rxn.de$logFC)                        
        if ("origin" %in% colnames(es$rxn.de)) {
            
            es$rxn.de <- es$rxn.de[es$rxn.de$ID %in% es$graph.raw$rxn, ]            
            
            lazyData("gene.id.map")
            unknown.rxn <- setdiff(es$rxn.de$ID, es$network$rxn2name$rxn)
            unknown.rxn2name <- do.call(cbind,  c(
                list(unknown.rxn), 
                rep(list(rep("", length(unknown.rxn))), ncol(es$network$rxn2name) - 1)))
            unknown.rxn2name <- as.data.frame(unknown.rxn2name, stringsAsFactors=F)
            colnames(unknown.rxn2name) <- colnames(es$network$rxn2name)
            
            es$network$rxn2name <- rbind(es$network$rxn2name, unknown.rxn2name)
            
            es$network$rxn2name$name[match(es$rxn.de$ID, es$network$rxn2name$rxn)] <-
                gene.id.map$name[match(es$rxn.de$origin, gene.id.map[,es$network$gene.ids])]            
            
            
        }
        
        es$rxn.pval <- es$rxn.de$pval
        names(es$rxn.pval) <- es$rxn.de$ID
        
        es$fb.rxn <- fitBumModel(es$rxn.pval, plot=F)
        if (plot) {
            hist(es$fb.rxn, main="Histogram of gene p-values")
            plot(es$fb.rxn, main="QQ-Plot for gene BUM-model")
        }
    } else {
        es$rxn.pval <- NULL
    }
        
    print("Building network")    
    
    met.de.ext <- data.frame(ID=unique(c(es$graph.raw$met.x, es$graph.raw$met.y)), stringsAsFactors=F)
    met.de.ext <- merge(met.de.ext, es$met.de, all.x=T) # all.x=T — we keep mets if there is no MS data
    met.de.ext <- merge(met.de.ext, es$network$met2name, by.x="ID", by.y="met", all.x=T)
    met.de.ext$logPval <- log(met.de.ext$pval)
    es$met.de.ext <- met.de.ext
    
    
    rxn.de.ext <- data.frame(ID=unique(es$graph.raw$rxn), stringsAsFactors=F)
    rxn.de.ext <- merge(rxn.de.ext, es$rxn.de, all.x=F) # we drop reaction if it's not expressed
    rxn.de.ext <- merge(rxn.de.ext, es$network$rxn2name, by.x="ID", by.y="rxn", all.x=T)
    rxn.de.ext$logPval <- log(rxn.de.ext$pval)
    es$rxn.de.ext <- rxn.de.ext
    
    es$graph.raw <- es$graph.raw[es$graph.raw$rxn %in% rxn.de.ext$ID,]
    es$graph.raw <- es$graph.raw[es$graph.raw$met.x %in% met.de.ext$ID,]
    es$graph.raw <- es$graph.raw[es$graph.raw$met.y %in% met.de.ext$ID,]
    
    
    if (reactions.as.edges) {    
        net.edges.ext <- merge(es$graph.raw, rxn.de.ext, by.x="rxn", by.y="ID")
        edges2rev <- net.edges.ext$met.x > net.edges.ext$met.y
        net.edges.ext[edges2rev, c("met.x", "met.y")] <- 
            net.edges.ext[edges2rev, c("met.y", "met.x")]        
        
        
        net.edges.pval <- (aggregate(pval ~ met.x + met.y, data=net.edges.ext, min))
        net.edges.ext <- merge(net.edges.pval, net.edges.ext)
        net.edges.ext <- net.edges.ext[!duplicated(net.edges.ext[,c("met.x", "met.y")]),]        
        net.edges.ext <- net.edges.ext[net.edges.ext$met.x != net.edges.ext$met.y,]    
        
        
        
        if (use.rpairs) {
            edges.keep <- net.edges.ext$rptype %in% "main"
        } else { 
            edges.keep <- TRUE
        }
        
        net.edges.ext.all <- net.edges.ext
        net.edges.ext <- net.edges.ext[edges.keep,]
        
        es$net.edges.ext <- net.edges.ext
        es$net.edges.ext.all <- net.edges.ext.all
        
        net1 <- graphNEL.from.tables(node.table=es$met.de.ext, edge.table=es$net.edges.ext,
                                     node.col="ID", edge.cols=c("met.x", "met.y"),
                                     directed=F)
        
        es$subnet <- net1
    } else {    
        
        edges <- rbind(rename(es$graph.raw[, c("met.x", "rxn")], c("met.x" = "met")),
                       rename(es$graph.raw[, c("met.y", "rxn")], c("met.y" = "met")))
        
        if (collapse.reactions) {
            print("Collapsing reactions by common most significant enzymes")
            
            edges <- data.table(edges, key="met")
            rxn.net <- merge(rename(edges, c("rxn"="rxn.x")), rename(edges, c("rxn" = "rxn.y")), by="met", allow.cartesian=T)
            rxn.net <- unique(rxn.net[, c("rxn.x", "rxn.y")])
            
            #rxn.de.origin.split <- es$rxn.de$origin
            print("before split")
            es$rxn.de.origin.split <- splitMappingByConnectivity(rxn.net, rxn.de.ext$ID, rxn.de.ext$origin)
            print("after split")
            
            t <- data.frame(from=es$rxn.de.ext$ID, to=es$rxn.de.origin.split, stringsAsFactors=F)                
            #from.table <- aggregate(from ~ to, t, function(x) paste(x, collapse="+"))
            
            es$rxn.de.ext$rxns <- es$rxn.de.ext$ID
            es$rxn.de.ext$ID <- es$rxn.de.origin.split            
            
            uniqueJoin <- function(x) {
                u  <-  unique(unlist(strsplit(as.character(x), split="\\+")))
                return(paste(u, collapse="+"))
            }
            
            print("before fixing")
            cols.to.fix <- setdiff(c(colnames(es$network$rxn2name), "rxns"), c("rxn", "name"))
            f <- as.formula(paste0("cbind(", paste0(cols.to.fix, collapse=","), ") ~ ID"))
            cols.fixed <- aggregate(f , es$rxn.de.ext, uniqueJoin)
            print("after fixing")
            
            es$rxn.de.ext <- unique(es$rxn.de.ext[, !names(es$rxn.de.ext) %in% cols.to.fix])
            es$rxn.de.ext <- merge(es$rxn.de.ext, cols.fixed)
            
            edges$rxn <- es$rxn.de.origin.split[edges$rxn]
            edges <- unique(edges)
            print("Done collapsing")
        }
        
        es$edges <- edges
        net1 <- graphNEL.from.tables(node.table=list(met=es$met.de.ext, rxn=es$rxn.de.ext), edge.table=edges,
                                     node.col="ID",
                                     directed=F)
        
        es$subnet <- net1
    }
    
    
    
    es$all.pval <- c(es$rxn.pval, es$met.pval)
    
    es$fb.all <- fitBumModel(es$all.pval, plot=F)
    if (plot) {
        hist(es$fb.rxn, main="Histogram of all p-values")
        plot(es$fb.rxn, main="QQ-Plot for joint BUM-model")
    }
    return(es)
}

appendModule <- function(res, module.graph) {        
    res[[length(res)+1]] <- module.graph
    res
}

#' Find significant module in the network
#' @param es Experiment set object
#' @param fdr FDR for both metabolites and genes/reactions, if set met.fdr and gene.fdr aren't used
#' @param met.fdr FDR for metabolites only
#' @param gene.fdr FDR for genes/reactions only
#' @param absent.met.score Score for metabolites absent from data
#' @param score.separately Score metabolites and reactions separately
#' @param heinz.py Path to heinz.py executable
#' @param heinz.nModules Number of modules to search for
#' @param heinz.tolerance tolerance parameter for heinz
#' @param heinz.subopt_diff subopt_diff parameter for heinz
#' @param simplify If TRUE and only one module was found return just the module, not a list.
#' @return List of most significant modules
#' @export 
findModule <- function(es,                         
                         fdr=NULL, met.fdr=NULL, gene.fdr=NULL,
                         absent.met.score=NULL,
                         score.separately=T,
                         heinz.py, heinz.nModules=1,
                         heinz.tolerance=10,
                         heinz.subopt_diff=100,
                         simplify=T) {
    
    if (!is.null(fdr)) {
        met.fdr <- fdr
        gene.fdr <- fdr
    }
    
            
    met.scores <- NULL    
    if (!is.null(es$met.de) && !is.null(met.fdr)) {            
        fb <- if (score.separately) es$fb.met else es$fb.all
        met.scores <- scoreValue(fb, na.omit(es$met.de.ext$pval), met.fdr)
        names(met.scores) <- es$met.de.ext$ID[!is.na(es$met.de.ext$pval)]
    }
    
    
    rxn.scores <- NULL
    if (!is.null(es$rxn.de) && !is.null(gene.fdr)) {
        fb <- if (score.separately) es$fb.rxn else es$fb.all
        rxn.scores <- scoreValue(fb, na.omit(es$rxn.de.ext$pval), gene.fdr)
        names(rxn.scores) <- es$rxn.de.ext$ID[!is.na(es$rxn.de.ext$pval)]
    }
    
    
        
    if (!is.null(met.scores)) {
        if (is.null(absent.met.score)) {
            absent.met.score <- mean(met.scores[met.scores < 0])
            print(paste0("absent.met.score <- ", absenet.met.score))
        }
        
        nodeDataDefaults(es$subnet, "score") <- absent.met.score
        met.scores <- met.scores[names(met.scores) %in% nodes(es$subnet)]
        nodeData(es$subnet, names(met.scores), "score") <- met.scores
    }
    
    if (!is.null(rxn.scores)) {
        if (es$reactions.as.edges) {
            edgeDataDefaults(es$subnet, "score") <- NA
            edgeData(es$subnet, from=es$net.edges.ext$met.x, to=es$net.edges.ext$met.y, attr="score") <- 
                rxn.scores[es$net.edges.ext$rxn]
        } else {
            rxn.scores <- rxn.scores[names(rxn.scores) %in% nodes(es$subnet)]
            nodeData(es$subnet, names(rxn.scores), "score") <- rxn.scores
        }
    }
    
    
    score.edges <- "score" %in% names(edgeDataDefaults(es$subnet))
    score.nodes <- "score" %in% names(nodeDataDefaults(es$subnet))
    
    if (!is.null(heinz.py)) {
        res <- runHeinz(
            subnet=es$subnet, 
            heinz.py=heinz.py, 
            score.edges=score.edges,
            score.nodes=score.nodes,
            heinz.nModules=heinz.nModules, 
            heinz.tolerance=heinz.tolerance,
            heinz.subopt_diff=heinz.subopt_diff)        
    } else if (score.edges) {
        stop("Can't use heuristic search with scored edges")
    } else {
        res <- list(runFastHeinz(es$subnet, unlist(nodeData(es$subnet, attr="score"))))
    }
        
    if (simplify && length(res) == 1) {
        return(res[[1]])
    }
    return(res)
}



runHeinz <- function(subnet,
                      heinz.py, 
                      score.edges=F,
                      score.nodes=T,                      
                      heinz.nModules=1, 
                      heinz.tolerance=10,
                      heinz.subopt_diff=100) {
    graph.dir <- tempfile("graph")
    dir.create(graph.dir)
    edges_file <- paste(graph.dir, "edges.txt", sep="/")
    nodes_file <- paste(graph.dir, "nodes.txt", sep="/")
    
    writeHeinzEdges(subnet, file=edges_file, use.score=score.edges)
    
    if (!score.nodes) {
        # Hack to make writeHeinzeNodes working
        nodeDataDefaults(subnet, "score")  <- 0
        nodeData(subnet, attr="score") <- 0
    }
    writeHeinzNodes(subnet, file=nodes_file, use.score=T)
    
    
    wd.bak <- getwd()
    heinz.dir <- dirname(heinz.py)
    setwd(heinz.dir)
    heinz.tmpdir <- tempfile("heinztmp")
    system2(paste0("./", basename(heinz.py)),
            c("-n", nodes_file,
              "-e", edges_file,
              "-N", if (score.nodes) "True" else "False",
              "-E", if (score.edges) "True" else "False",              
              "-s", heinz.nModules,
              "--tolerance", heinz.tolerance,
              "--subopt_diff", heinz.subopt_diff,
              "--heinztmp", heinz.tmpdir,
              "-v"),
            env="ILOG_LICENSE_FILE=./access.ilm")
    setwd(wd.bak)
    
    
    res <- list()
    for (i in 0:(heinz.nModules-1)) {
        module.graph = readHeinzGraph(node.file = paste(nodes_file, i, "hnz", sep="."), 
                                      network = subnet)
        res <- appendModule(res, module.graph)
        
    }
    return(res)
}

#' Replace simple reaction nodes with edges
#' For every reaction node that have just two connections
#' with metabolides on different side of the reaction
#' node for the reaction is replaced with direct edge
#' between these two metabolites
#' @param module Module to modify
#' @param es Experiment set object
#' @return Modified module
#' @export
removeSimpleReactions <- function(module, es) {
    stopifnot(!es$reactions.as.edges)
    rxn.nodes <- nodes(module)[unlist(nodeData(module, attr="nodeType")) == "rxn"]
    rxn.edges <- edges(module, rxn.nodes)
    bm.rxns <- names(rxn.edges)[sapply(rxn.edges, length) == 2]
    original.rxns <- unlist(nodeData(module, bm.rxns, attr="rxns"))
    
    res <- module
    for (attr in names(nodeDataDefaults(res))) {
        if (!attr %in% edgeDataDefaults(res)) {
            edgeDataDefaults(res, attr) <- nodeDataDefaults(res, attr)
        }
    }
    
    for (new.rxn in bm.rxns) {
        is.reaction <- F
        c1 <- rxn.edges[[new.rxn]][[1]]
        c2 <- rxn.edges[[new.rxn]][[2]]
        
        for (old.rxn in names(es$rxn.de.origin.split[es$rxn.de.origin.split == new.rxn])) {
            is.reaction <- 
                is.reaction || 
                any((es$network$graph.raw$met.x == c1) & 
                    (es$network$graph.raw$met.y == c2) & 
                    (es$network$graph.raw$rxn == old.rxn))
            
            is.reaction <- 
                is.reaction || 
                any((es$network$graph.raw$met.x == c2) & 
                    (es$network$graph.raw$met.y == c1) & 
                    (es$network$graph.raw$rxn == old.rxn))
        }
        if (is.reaction) {
            if (!c1 %in% unlist(edges(res, c2)) && !c2 %in% unlist(edges(res, c1))) {
                res <- addEdge(from=c1, to=c2, res)
            }
            # :ToDo: what if there already was an edge
            if (is.na(edgeData(res, from=c1, to=c2, attr="pval")) || 
                    unlist(edgeData(res, from=c1, to=c2, attr="pval")) > unlist(nodeData(res, new.rxn, "pval"))) {
                
                for (attr in names(nodeDataDefaults(res))) {
                    edgeData(res, from=c1, to=c2, attr=attr) <- nodeData(res, new.rxn, attr)
                }
            }
            
            res <- removeNode(new.rxn, res)
        }
        
    }
    return(res)
}

#' Add all metabolites connected with reactions
#' For every reaction node that have connections only with
#' metabolites on one side of the reaction connections are
#' added to metabolites on the other side of the reaction.
#' @param module Module to modify
#' @param es Experiment set object
#' @return Modified module
#' @export
addMetabolitesForReactions<- function(module, es) {
    stopifnot(!es$reactions.as.edges)
    rxn.nodes <- nodes(module)[unlist(nodeData(module, attr="nodeType")) == "rxn"]
    rxn.edges <- edges(module, rxn.nodes)
    leaf.rxns <- names(rxn.edges)[sapply(rxn.edges, length) == 1]
    original.rxns <- unlist(nodeData(module, leaf.rxns, attr="rxns"))
    
    res <- module
    
    for (new.rxn in rxn.nodes) {
        
        cs <- rxn.edges[[new.rxn]]
        
        old.rxns <- names(es$rxn.de.origin.split[es$rxn.de.origin.split == new.rxn])
        
        candidate.rxns <- c()
        for (old.rxn in old.rxns) {
            rxn.net <- es$network$graph.raw[es$network$graph.raw$rxn == old.rxn,]
            net.cs <- unique(c(rxn.net$met.x, rxn.net$met.y))
            
            if (length(setdiff(cs, net.cs)) == 0) {
                #candidate.rxns <- c(candidate.rxns, old.rxn)
            }
            candidate.rxns <- c(candidate.rxns, old.rxn)

            if (any((rxn.net$met.x %in% cs) & (rxn.net$met.y %in% cs))) {
                next
            }
        }
        
        if (length(candidate.rxns) > 1) {
            #next
        }
        
        for (old.rxn in candidate.rxns) {
            rxn.net <- es$network$graph.raw[es$network$graph.raw$rxn == old.rxn,]
            rxn.net <- rbind(rxn.net, rename(rxn.net, c("met.x"="met.y", "met.y"="met.x")))
            
            c2s <- unique(c(rxn.net$met.x, rxn.net$met.y))
            
            for (c2 in c2s) {
                if (!c2 %in% nodes(res)) {
                    res <- addNode(c2, res)
                    nodeData(res, c2, attr="nodeType") <- "met"
                }
                if (!new.rxn %in% unlist(edges(res, c2))) {
                    res <- addEdge(new.rxn, c2, res)
                }
            }
        }
        
    }
    res <- addNodeAttributes(res, list(met=es$met.de.ext))
    return(res)
}

#' Add reactions that connect metabolites in module
#' Metabolites are connected by reaction if they are on different sides
#' @param module Module to work with
#' @param es Experiment set object
#' @return Module with interconnecting reactions
#' @export
addInterconnections <- function(module, es) {
    met.nodes <- nodes(module)[unlist(nodeData(module, attr="nodeType")) == "met"]
    interconnects <- es$graph.raw$rxn[es$graph.raw$met.x %in% met.nodes & es$graph.raw$met.y %in% met.nodes]
    interconnects <- unique(es$rxn.de.origin.split[interconnects])
    
    new.nodes <- unique(c(nodes(module), interconnects))
    res <- subNetwork(new.nodes, es$subnet)
    return(res)
}

#' Copy attributes from reaction node to adjacent edges
#' @param module Module do work with
#' @return Module with attributes copied
#' @export
expandReactionNodeAttributesToEdges <- function(module) {
    rxn.nodes <- nodes(module)[unlist(nodeData(module, attr="nodeType")) == "rxn"]
    
    res <- module
    for (attr in names(nodeDataDefaults(res))) {
        if (!attr %in% edgeDataDefaults(res)) {
            edgeDataDefaults(res, attr) <- nodeDataDefaults(res, attr)
        }
    }
    
    for (rxn.node in rxn.nodes) {
        for (met.node in unlist(edges(res, rxn.node))) {
            for (attr in names(nodeDataDefaults(res))) {
                edgeData(res, from=rxn.node, to=met.node, attr=attr) <-
                    nodeData(res, rxn.node, attr)
            }
        }
    }
    return(res)
}

#' Remove hanging nodes without data
#' @param module Module to work with
#' @return Module with hanging nodes removed
#' @export
removeHangingNodes <- function(module) {
    absent.data.nodes <- nodes(module)[is.na(unlist(nodeData(module, attr="pval")))]
    edges <- edgelist(module)
    absent.data.edges <- edges[edges$u %in% absent.data.nodes,]
    absent.nodes.degree <- table(absent.data.edges$u)
    hanging.nodes <- names(absent.nodes.degree)[absent.nodes.degree == 1]
    res <- subNetwork(setdiff(nodes(module), hanging.nodes), module)
    
    return(res)
    
}

#' Add reaction trans- edges connecting metabolites in module
#' @param module Module to work with
#' @param es Experiment set object
#' @return Modified module
#' @export
addTransEdges <- function(module, es) {
    edges.keep <- es$net.edges.ext.all$rptype %in% c("main", "trans")
    net.with.trans <- graphNEL.from.tables(node.table=es$met.de.ext, edge.table=es$net.edges.ext.all[edges.keep,],
                                 node.col="ID", edge.cols=c("met.x", "met.y"),
                                 directed=F)
    return(subNetwork(nodes(module), net.with.trans))
}
