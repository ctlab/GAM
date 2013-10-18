#' @import BioNet 
NULL

setdiff.data.frame <-
    function(A,B) A[ !duplicated( rbind(B,A) )[ -seq_len(nrow(B))] , ]

#' Load data only if its absent from global environment
#' @param name Name of a dataset
#' @param ... Additional arguments for data()
lazyData <- function(name, ...) {
    if (!name %in% ls(envir=.GlobalEnv)) {
        print(paste0("No ", name, ", loading"))
        do.call(data, list(name, ...))
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
#' @importFrom igraph E V list.vertex.attributes layout.kamada.kawai layout.norm get.edges plot.igraph get.vertex.attribute
#' @export
plotNetwork <- function(module, scale=1, attr.label="label", attr.shape="nodeType", layout=layout.kamada.kawai, ...) {
    network <- module
    
    if (is.null(V(network)$name)) {
        V(network)$name <- as.character(V(network))
    }
    # Hack for coloring
    
    if ("logFC.norm" %in% list.vertex.attributes(network)) {
        de <- V(network)$logFC.norm
        de <- de * 10
    } else {
        de <- V(network)$logFC
    }
    names(de) <- V(network)$name
    de[is.na(de)] <- 0
    
    
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
    node.types <- unique(na.omit(get.vertex.attribute(module, group.by)))
        
    for (node.type in node.types) {
        module.type.nodes <- V(module)[get.vertex.attribute(module, group.by) == node.type]
        module.type.de <- get.vertex.attribute(module, logFC.attr, module.type.nodes)
        module.type.norm.de <- module.type.de / max(abs(na.omit(module.type.de)))
        module <- set.vertex.attribute(module, logFC.norm.attr, module.type.nodes, module.type.norm.de)
    }
    
    module
}

# Scores p-value by fitted BUM-model
# This is a helper function based on BioNet::scoreFunction()
scoreValue <- function (fb, pval, fdr = 0.01) 
{
    return((fb$a - 1) * (log(pval) - log(fdrThreshold(fdr, fb))))
}


#' Preprocess experiment set's differential expression data for genes and metabolites
#' 
#' Convert metabolite and gene IDs to be the same as in network. Computes
#' reaction differential expression data.
#' @param es Experiment set with DE data
#' @param met.ids Type of IDs used in metabolite DE data (see met.id.map for possible values)
#' @param gene.ids Type of IDs used in gene DE data (see gene.id.map for possible values)
#' @param plot If TRUE plot BUM fit
preprocessPvalAndMetDE <- function(es, met.ids, gene.ids, plot=T) {
    if (!is.null(es$gene.de)) {
        print("Processing gene p-values...")
        es$gene.de$logFC <- fixInf(es$gene.de$logFC)
        if (!is.null(gene.ids)) {
            lazyData("gene.id.map")
            gene.id.map <- get("gene.id.map")
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
            met.id.map <- get("met.id.map")
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
            gene.id.map <- get("gene.id.map")
            
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
    
    if (is.null(es$met.de)) {
        es$met.de <- data.frame(
            ID=unique(c(es$graph.raw$met.x, es$graph.raw$met.y)), 
            pval=NA,
            logFC=NA,
            stringsAsFactors=F)
    }
    
    met.de.ext <- data.frame(ID=unique(c(es$graph.raw$met.x, es$graph.raw$met.y)), stringsAsFactors=F)
    met.de.ext <- merge(met.de.ext, es$met.de, all.x=T) # all.x=T — we keep mets if there is no MS data
    met.de.ext <- merge(met.de.ext, es$network$met2name, by.x="ID", by.y="met", all.x=T)
    met.de.ext$logPval <- log(met.de.ext$pval)
    es$met.de.ext <- met.de.ext
    
    
    if (is.null(es$rxn.de)) {
        es$rxn.de <- data.frame(ID=es$graph.raw$rxn, pval=NA, logFC=NA, stringsAsFactors=F)
    }
    
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
        
        
        net.edges.pval <- (aggregate(pval ~ met.x + met.y, data=net.edges.ext, min, na.action=na.pass))
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
        
        net1 <- graph.from.tables(node.table=list(met=es$met.de.ext), edge.table=es$net.edges.ext,
                                     node.col="ID", edge.cols=c("met.x", "met.y"),
                                     directed=F)
        
        es$subnet <- net1
    } else {    
        
        edges <- rbind(rename(es$graph.raw[, c("met.x", "rxn")], c("met.x" = "met")),
                       rename(es$graph.raw[, c("met.y", "rxn")], c("met.y" = "met")))
        
        if (collapse.reactions && "origin" %in% names(rxn.de.ext)) {
            print("Collapsing reactions by common most significant enzymes")
            
            edges.dt <- data.table(edges, key="met")
            rxn.net <- merge(rename(edges.dt, c("rxn"="rxn.x")), rename(edges.dt, c("rxn" = "rxn.y")), by="met", allow.cartesian=T)
            rxn.net <- unique(rxn.net[, list(rxn.x, rxn.y)])
            
            #rxn.de.origin.split <- es$rxn.de$origin
            print("before split")
            es$rxn.de.origin.split <- GAM:::splitMappingByConnectivity(rxn.net, rxn.de.ext$ID, rxn.de.ext$origin)
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
        net1 <- graph.from.tables(node.table=list(met=es$met.de.ext, rxn=es$rxn.de.ext), edge.table=edges,
                                     node.col="ID",
                                     directed=F)
        
        es$subnet <- net1
    }
    
    
    
    es$all.pval <- c(es$rxn.pval, es$met.pval)
    
    es$fb.all <- fitBumModel(es$all.pval, plot=F)
    if (plot) {
        hist(es$fb.all, main="Histogram of all p-values")
        plot(es$fb.all, main="QQ-Plot for joint BUM-model")
    }
    return(es)
}

appendModule <- function(res, module.graph) {        
    res[[length(res)+1]] <- module.graph
    res
}

runHeinz <- function(subnet,
                      heinz.py, 
                      score.edges=F,
                      score.nodes=T,                      
                      nModules=1, 
                      tolerance=10,
                      subopt_diff=100,
                      cplexTimeLimit=1e+75) {
    
    graph.dir <- tempfile("graph")
    dir.create(graph.dir)
    edges_file <- paste(graph.dir, "edges.txt", sep="/")
    nodes_file <- paste(graph.dir, "nodes.txt", sep="/")
    
    writeHeinzEdges(subnet, file=edges_file, use.score=score.edges)
    
    if (!score.nodes) {
        # Hack to make writeHeinzeNodes working
        V(subnet)$score <- 0
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
              "-s", nModules,
              "--tolerance", tolerance,
              "--subopt_diff", subopt_diff,
              "--heinztmp", heinz.tmpdir,
              "--additional", paste0("'cplexTimeLimit ", cplexTimeLimit, "'"),
              "-v"),
            env="ILOG_LICENSE_FILE=./access.ilm")
    setwd(wd.bak)
    
    
    res <- list()
    for (i in 0:(nModules-1)) {
        module.graph <- readHeinzGraph(node.file = paste(nodes_file, i, "hnz", sep="."), 
                                      network = subnet, format="igraph")
        res <- appendModule(res, module.graph)
        
    }
    return(res)
}

#' Solves MWCS using heinz
#' @param heinz.py Path to heinz.py executable
#' @param nModules Number of modules to search for
#' @param tolerance tolerance parameter for heinz
#' @param subopt_diff subopt_diff parameter for heinz
#' @param timeLimit Time limit for execution
#' @return solver function
#' @export
heinz.solver <- function(heinz.py,
                         nModules=1,
                         tolerance=10,
                         subopt_diff=100,
                         timeLimit=1e+75
                         ) {
    
    
    function(network) {
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        
        res <- runHeinz(
            subnet=network, 
            heinz.py=heinz.py, 
            score.edges=score.edges,
            score.nodes=score.nodes,
            nModules=nModules, 
            tolerance=tolerance,
            subopt_diff=subopt_diff,
            cplexTimeLimit=timeLimit
            )        
    }
}

#' Solves MWCS with BioNet::runFastHeinz algorithm
#' @param network Netowrk to find module in
#' @return Module
#' @export
fastHeinz.solver <- function(network) {
    score.edges <- "score" %in% list.edge.attributes(network) && !all(E(network)$score == 0)
    if (score.edges) {
        stop("Can't run fast heinz on network with scored edges")
    }
    res <- list(runFastHeinz(network, V(net)$score))
}

#' Find significant module in the network
#' @param es Experiment set object
#' @param fdr FDR for both metabolites and genes/reactions, if set met.fdr and gene.fdr aren't used
#' @param met.fdr FDR for metabolites only
#' @param gene.fdr FDR for genes/reactions only
#' @param absent.met.score Score for metabolites absent from data
#' @param absent.rxn.score Score for reactions when there is no genomic data
#' @param score.separately Score metabolites and reactions separately
#' @param solver Solver function of MWCS problem to use, first argument should be a network,
#'               result should be a module or list of modules
#' @param simplify If TRUE and only one module was found return just the module, not a list.
#' @param ... Additional arguments for solver
#' @return List of most significant modules
#' @importFrom igraph V<- E<-
#' @export 
findModule <- function(es,                         
                       fdr=NULL, met.fdr=NULL, gene.fdr=NULL,
                       absent.met.score=NULL,
                       absent.rxn.score=0,
                       score.separately=T,
                       solver = fastHeinz.solver,
                       simplify=T,
                       ...) {
    
    net <- es$subnet
    
    if (!is.null(fdr)) {
        met.fdr <- fdr
        gene.fdr <- fdr
    }
    
            
    met.scores <- NULL    
    if (!is.null(es$fb.met) && !is.null(met.fdr)) {            
        fb <- if (score.separately) es$fb.met else es$fb.all
        met.scores <- scoreValue(fb, na.omit(es$met.de.ext$pval), met.fdr)
        names(met.scores) <- es$met.de.ext$ID[!is.na(es$met.de.ext$pval)]
    }
    
    
    rxn.scores <- NULL
    if (!is.null(es$fb.rxn) && !is.null(gene.fdr)) {
        fb <- if (score.separately) es$fb.rxn else es$fb.all
        rxn.scores <- scoreValue(fb, na.omit(es$rxn.de.ext$pval), gene.fdr)
        names(rxn.scores) <- es$rxn.de.ext$ID[!is.na(es$rxn.de.ext$pval)]
    }
    
    
    if (is.null(absent.met.score)) {
        absent.met.score <- mean(met.scores[met.scores < 0])
        print(paste0("absent.met.score <- ", absent.met.score))
    }
        
    absent.met.scores <- sapply(V(net)[nodeType == "met"]$name, function(x) absent.met.score)
    met.scores <- c(met.scores, absent.met.scores[!names(absent.met.scores) %in% names(met.scores)])
    
    met.scores <- met.scores[names(met.scores) %in% V(net)$name]
    V(net)[names(met.scores)]$score <- met.scores
    
    absent.rxn.scores <- sapply(es$rxn.de$ID, function(x) absent.rxn.score)
    rxn.scores <- c(rxn.scores, absent.rxn.scores[!names(absent.rxn.scores) %in% names(rxn.scores)])
    
    if (es$reactions.as.edges) {
        
        E(net)$score <- rxn.scores[E(net)$rxn]
    } else {
        rxn.scores <- rxn.scores[names(rxn.scores) %in% V(net)$name]
        V(net)[names(rxn.scores)]$score <- rxn.scores
    }
    
    
    res <- solver(net, ...)
    
    if (simplify && is(res, "list") && length(res) == 1) {
        return(res[[1]])
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
#' @importFrom igraph add.edges delete.edges delete.vertices degree
#' @export
removeSimpleReactions <- function(module, es) {
    stopifnot(!es$reactions.as.edges)
    res <- module
    
    rxn.nodes <- V(res)[nodeType == "rxn" & degree(res) == 2]$name
    rxn.edges <- get.edges(res, E(res)[adj(rxn.nodes)])
    rxn.edges <- matrix(V(res)[rxn.edges]$name, ncol=2)
    rxn.edges.types <- matrix(V(res)[rxn.edges]$nodeType, ncol=2)
    rxn.edges[rxn.edges.types[,2] == "rxn"] <- rxn.edges[rxn.edges.types[,2] == "rxn", c(2, 1)]
    
    simple.edges <-  data.frame(do.call(rbind, split(rxn.edges[,2], rxn.edges[,1])), stringsAsFactors=F)
    colnames(simple.edges) <- c("met.x", "met.y")
    simple.edges$rxn <- rownames(simple.edges)
    
    new2old <- data.frame(new=es$rxn.de.origin.split, old=names(es$rxn.de.origin.split))
    simple.edges.ext <- merge(simple.edges, new2old, by.x = "rxn", by.y="new")
    
    simple.edges.pasted.1 <- do.call(paste, c(simple.edges.ext[, c("met.x", "old", "met.y")], sep="\r"))
    simple.edges.pasted.2 <- do.call(paste, c(simple.edges.ext[, c("met.y", "old", "met.x")], sep="\r"))
    
    graph.raw.pasted <- do.call(paste, c(es$network$graph.raw[, c("met.x", "rxn", "met.y")], sep="\r"))
    
    
    is.simple <- simple.edges.pasted.1 %in% graph.raw.pasted | simple.edges.pasted.2 %in% graph.raw.pasted
    simple.edges.ext <- simple.edges.ext[is.simple,]
    
    simple.edges <- simple.edges[simple.edges$rxn %in% simple.edges.ext$rxn,]
    simple.edges <- simple.edges[order(V(res)[simple.edges$rxn]$pval),]
    
    res <- add.edges(
        res, 
        rbind(simple.edges$met.x, simple.edges$met.y),
        attr=get.vertex.attributes(res, simple.edges$rxn)
        )
    res <- delete.edges(res, E(res, P=rbind(simple.edges$rxn, simple.edges$met.x)))
    res <- delete.edges(res, E(res, P=rbind(simple.edges$rxn, simple.edges$met.y)))
    res <- delete.vertices(res, simple.edges$rxn)
    res <- simplify(res, edge.attr.comb="first")
    return(res)
}

#' Add all metabolites connected with reactions
#' @param module Module to modify
#' @param es Experiment set object
#' @return Modified module
#' @importFrom igraph add.vertices
#' @export
addMetabolitesForReactions<- function(module, es) {
    stopifnot(!es$reactions.as.edges)
    res <- module
    
    rxn.nodes <- V(res)[nodeType == "rxn"]$name
    rxn.edges <- get.edges(res, E(res)[adj(rxn.nodes)])
    rxn.edges <- matrix(V(res)[rxn.edges]$name, ncol=2)
    rxn.edges.types <- matrix(V(res)[rxn.edges]$nodeType, ncol=2)
    rxn.edges[rxn.edges.types[,2] == "rxn"] <- rxn.edges[rxn.edges.types[,2] == "rxn", c(2, 1)]
    rxn.edges <- data.frame(rxn.edges, stringsAsFactors=F)
    colnames(rxn.edges) <- c("rxn", "met")
    
    new2old <- data.frame(new=es$rxn.de.origin.split, old=names(es$rxn.de.origin.split), stringsAsFactors=F)
    new2old <- new2old[new2old$new %in% rxn.nodes,]
    
    all.rxn.edges <- merge(es$network$graph.raw, new2old, by.x="rxn", by.y="old")
    all.rxn.edges <- 
        rbind(
            rename(all.rxn.edges[, c("new", "met.x")], c("met.x"="met", "new"="rxn")),
            rename(all.rxn.edges[, c("new", "met.y")], c("met.y"="met", "new"="rxn"))
        )
    all.rxn.edges <- unique(all.rxn.edges)
    new.edges <- setdiff.data.frame(all.rxn.edges, rxn.edges)
    new.vertices <- setdiff(new.edges$met, V(res)$name)
    res <- add.vertices(res, length(new.vertices), attr=list(name=new.vertices))
    res <- add.edges(
        res, 
        rbind(new.edges$rxn, new.edges$met),
        attr=list(weight=1)
    )
        
    res <- addNodeAttributes(res, list(met=es$met.de.ext))
    return(res)
}

#' Add reactions that connect metabolites in module
#' Metabolites are connected by reaction if they are on different sides
#' @param module Module to work with
#' @param es Experiment set object
#' @return Module with interconnecting reactions
#' @importFrom igraph induced.subgraph
#' @export
addInterconnections <- function(module, es) {
    met.nodes <- V(module)[nodeType == "met"]$name
    interconnects <- es$graph.raw$rxn[es$graph.raw$met.x %in% met.nodes & es$graph.raw$met.y %in% met.nodes]
    interconnects <- unique(es$rxn.de.origin.split[interconnects])
    
    new.nodes <- unique(c(V(module)$name, interconnects))
    res <- induced.subgraph(es$subnet, new.nodes)
    return(res)
}

#' Copy attributes from reaction node to adjacent edges
#' @param module Module do work with
#' @return Module with attributes copied
#' @importFrom igraph set.edge.attribute
#' @export
expandReactionNodeAttributesToEdges <- function(module) {
    res <- module
    
    rxn.nodes <- V(res)[nodeType == "rxn"]
    
    edge.ids <- E(res)[adj(rxn.nodes)]
    edges <- get.edges(res, edge.ids)
    z <- edges[, 2] %in% rxn.nodes
    edges[z, ] <- edges[z, 2:1] 
    
    for (attr in list.vertex.attributes(res)) {
        res <- set.edge.attribute(res, attr, edge.ids, get.vertex.attribute(res, attr, edges[,1]))
    }
    return(res)
}

#' Remove hanging nodes without data
#' @param module Module to work with
#' @return Module with hanging nodes removed
#' @export
removeHangingNodes <- function(module) {
    res <- module
    res <- delete.vertices(res, V(res)[degree(res) == 1 & is.na(pval)])
    return(res)
    
}

#' Add reaction trans- edges connecting metabolites in module
#' @param module Module to work with
#' @param es Experiment set object
#' @return Modified module
#' @export
addTransEdges <- function(module, es) {
    edges.keep <- es$net.edges.ext.all$rptype %in% c("main", "trans")
    net.with.trans <- graph.from.tables(
        node.table=es$met.de.ext,
        edge.table=es$net.edges.ext.all[edges.keep,],
        node.col="ID",
        edge.cols=c("met.x", "met.y"),
        directed=F)
    return(induced.subgraph(net.with.trans, V(module)$name))
}
