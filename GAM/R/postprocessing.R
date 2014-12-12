#' Add node attribute with normalized log-foldchange inside node groups.
#' 
#' For each node assigns normalized log-foldchange values equal to 
#' log-foldchange divided by maximum of absolute log-foldchange among
#' the same node type.
#' 
#' @param module Module to add attribute to
#' @param logFC.attr Attribute with node log-foldchange
#' @param logFC.norm.attr Name of a new attribute
#' @param group.by Attribute by which node grouping shoud happen
#' @return Modified module with normalized log-foldchange node attribute
#' @examples
#' library("mouseMacrophages")
#' data("examplesGAM")
#' 
#' data("kegg.mouse.network")
#' 
#' es.rn <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M1.M2,
#'                            gene.de=gene.de.M1.M2,
#'                            reactions.as.edges=FALSE,
#'                            plot=FALSE)
#'     
#' \dontrun{
#' heinz.py <- "/usr/local/lib/heinz/heinz.py"
#' solver <- heinz.solver(heinz.py)
#' 
#' module.rn <- findModule(es.rn,
#'                         met.fdr=1e-6,
#'                         rxn.fdr=1e-6,
#'                         absent.met.score=-20,
#'                         solver=solver)
#' }
#' 
#' module.rn <- addMetabolitesForReactions(module.rn, es.rn)
#' module.rn <- addInterconnections(module.rn, es.rn)
#' module.rn <- addNormLogFC(module.rn)
#' module.rn <- removeHangingNodes(module.rn)
#' module.rn <- simplifyReactionNodes(module.rn, es.rn)
#' module.rn <- expandReactionNodeAttributesToEdges(module.rn)
#' @export
addNormLogFC <- function(module, logFC.attr="log2FC", logFC.norm.attr="log2FC.norm", group.by="nodeType") {
    node.types <- unique(na.omit(get.vertex.attribute(module, group.by)))
        
    for (node.type in node.types) {
        module.type.nodes <- V(module)[get.vertex.attribute(module, group.by) == node.type]
        module.type.de <- get.vertex.attribute(module, logFC.attr, module.type.nodes)
        module.type.de <- fixInf(module.type.de)
        module.type.norm.de <- module.type.de / max(abs(na.omit(module.type.de)))
        module <- set.vertex.attribute(module, logFC.norm.attr, module.type.nodes, module.type.norm.de)
    }
    
    module
}

#' Replace simple reaction nodes with edges
#' 
#' For every reaction node that have just two connections
#' with metabolides on different side of the reaction,
#' a node for the reaction is replaced with a direct edge
#' between these two metabolites.
#' 
#' @param module Module to modify
#' @param es Experiment set object
#' @return Modified module
#' @import igraph 
#' @examples
#' library("mouseMacrophages")
#' data("examplesGAM")
#' 
#' data("kegg.mouse.network")
#' 
#' es.rn <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M1.M2,
#'                            gene.de=gene.de.M1.M2,
#'                            reactions.as.edges=FALSE,
#'                            plot=FALSE)
#'     
#' \dontrun{
#' heinz.py <- "/usr/local/lib/heinz/heinz.py"
#' solver <- heinz.solver(heinz.py)
#' 
#' module.rn <- findModule(es.rn,
#'                         met.fdr=1e-6,
#'                         rxn.fdr=1e-6,
#'                         absent.met.score=-20,
#'                         solver=solver)
#' }
#' 
#' module.rn <- addMetabolitesForReactions(module.rn, es.rn)
#' module.rn <- addInterconnections(module.rn, es.rn)
#' module.rn <- addNormLogFC(module.rn)
#' module.rn <- removeHangingNodes(module.rn)
#' module.rn <- simplifyReactionNodes(module.rn, es.rn)
#' module.rn <- expandReactionNodeAttributesToEdges(module.rn)
#' @export
simplifyReactionNodes <- function(module, es) {
    stopifnot(!es$reactions.as.edges)
    res <- module
    
    rxn.nodes <- V(res)[nodeType == "rxn" & degree(res) == 2]$name
    rxn.edges <- get.edges(res, E(res)[adj(rxn.nodes)])
    rxn.edges <- matrix(V(res)[rxn.edges]$name, ncol=2)
    rxn.edges.types <- matrix(V(res)[rxn.edges]$nodeType, ncol=2)
    rxn.edges[rxn.edges.types[,2] == "rxn"] <- rxn.edges[rxn.edges.types[,2] == "rxn", c(2, 1)]
    
    simple.edges <-  data.frame(do.call(rbind, split(rxn.edges[,2], rxn.edges[,1])), stringsAsFactors=FALSE)
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
#' 
#' For each reaction node adds all metabolites that
#' are part of this reaction.
#' 
#' @param module Module to modify
#' @param es Experiment set object
#' @return Modified module
#' @import igraph 
#' @examples
#' library("mouseMacrophages")
#' data("examplesGAM")
#' 
#' data("kegg.mouse.network")
#' 
#' es.rn <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M1.M2,
#'                            gene.de=gene.de.M1.M2,
#'                            reactions.as.edges=FALSE,
#'                            plot=FALSE)
#'     
#' \dontrun{
#' heinz.py <- "/usr/local/lib/heinz/heinz.py"
#' solver <- heinz.solver(heinz.py)
#' 
#' module.rn <- findModule(es.rn,
#'                         met.fdr=1e-6,
#'                         rxn.fdr=1e-6,
#'                         absent.met.score=-20,
#'                         solver=solver)
#' }
#' 
#' module.rn <- addMetabolitesForReactions(module.rn, es.rn)
#' module.rn <- addInterconnections(module.rn, es.rn)
#' module.rn <- addNormLogFC(module.rn)
#' module.rn <- removeHangingNodes(module.rn)
#' module.rn <- simplifyReactionNodes(module.rn, es.rn)
#' module.rn <- expandReactionNodeAttributesToEdges(module.rn)
#' @export
addMetabolitesForReactions<- function(module, es) {
    stopifnot(!es$reactions.as.edges)
    res <- module
    
    rxn.nodes <- V(res)[nodeType == "rxn"]$name
    rxn.edges <- get.edges(res, E(res)[adj(rxn.nodes)])
    rxn.edges <- matrix(V(res)[rxn.edges]$name, ncol=2)
    rxn.edges.types <- matrix(V(res)[rxn.edges]$nodeType, ncol=2)
    rxn.edges[rxn.edges.types[,2] == "rxn"] <- rxn.edges[rxn.edges.types[,2] == "rxn", c(2, 1)]
    rxn.edges <- data.frame(rxn.edges, stringsAsFactors=FALSE)
    colnames(rxn.edges) <- c("rxn", "met")
    
    new2old <- data.frame(new=es$rxn.de.origin.split, old=names(es$rxn.de.origin.split), stringsAsFactors=FALSE)
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

#' Add reactions that connect metabolites in a module
#' 
#' Adds nodes for reactions that have at least one metabolite from the module
#' on both sides. Works only when reactions are nodes.
#' 
#' @param module Module to work with
#' @param es Experiment set object
#' @return Module with interconnecting reactions
#' @import igraph 
#' @examples
#' library("mouseMacrophages")
#' data("examplesGAM")
#' 
#' data("kegg.mouse.network")
#' 
#' es.rn <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M1.M2,
#'                            gene.de=gene.de.M1.M2,
#'                            reactions.as.edges=FALSE,
#'                            plot=FALSE)
#'     
#' \dontrun{
#' heinz.py <- "/usr/local/lib/heinz/heinz.py"
#' solver <- heinz.solver(heinz.py)
#' 
#' module.rn <- findModule(es.rn,
#'                         met.fdr=1e-6,
#'                         rxn.fdr=1e-6,
#'                         absent.met.score=-20,
#'                         solver=solver)
#' }
#' 
#' module.rn <- addMetabolitesForReactions(module.rn, es.rn)
#' module.rn <- addInterconnections(module.rn, es.rn)
#' module.rn <- addNormLogFC(module.rn)
#' module.rn <- removeHangingNodes(module.rn)
#' module.rn <- simplifyReactionNodes(module.rn, es.rn)
#' module.rn <- expandReactionNodeAttributesToEdges(module.rn)
#' @export
addInterconnections <- function(module, es) {
    stopifnot(all(V(module)$name %in% V(es$subnet)$name))
    met.nodes <- V(module)[nodeType == "met"]$name
    interconnects <- es$graph.raw$rxn[es$graph.raw$met.x %in% met.nodes & es$graph.raw$met.y %in% met.nodes]
    interconnects <- unique(es$rxn.de.origin.split[interconnects])
    
    new.nodes <- unique(c(V(module)$name, interconnects))
    res <- induced.subgraph(es$subnet, new.nodes)
    return(res)
}

#' Copy attributes from reaction node to adjacent edges
#' 
#' Copies reaction node attribues to adjacent edges. Useful
#' for visualization purposes after calling \code{simplifyReactionNodes}
#' function, so that all edges will have attributes.
#' 
#' @param module Module do work with
#' @return Module with attributes copied
#' @import igraph 
#' @examples
#' library("mouseMacrophages")
#' data("examplesGAM")
#' 
#' data("kegg.mouse.network")
#' 
#' es.rn <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M1.M2,
#'                            gene.de=gene.de.M1.M2,
#'                            reactions.as.edges=FALSE,
#'                            plot=FALSE)
#'     
#' \dontrun{
#' heinz.py <- "/usr/local/lib/heinz/heinz.py"
#' solver <- heinz.solver(heinz.py)
#' 
#' module.rn <- findModule(es.rn,
#'                         met.fdr=1e-6,
#'                         rxn.fdr=1e-6,
#'                         absent.met.score=-20,
#'                         solver=solver)
#' }
#' 
#' module.rn <- addMetabolitesForReactions(module.rn, es.rn)
#' module.rn <- addInterconnections(module.rn, es.rn)
#' module.rn <- addNormLogFC(module.rn)
#' module.rn <- removeHangingNodes(module.rn)
#' module.rn <- simplifyReactionNodes(module.rn, es.rn)
#' module.rn <- expandReactionNodeAttributesToEdges(module.rn)
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
#' 
#' Remove nodes with degree equal to one that don't have known p-value.
#' 
#' @param module Module to work with
#' @return Module with hanging nodes removed
#' @examples
#' library("mouseMacrophages")
#' data("examplesGAM")
#' 
#' data("kegg.mouse.network")
#' 
#' es.rn <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M1.M2,
#'                            gene.de=gene.de.M1.M2,
#'                            reactions.as.edges=FALSE,
#'                            plot=FALSE)
#'     
#' \dontrun{
#' heinz.py <- "/usr/local/lib/heinz/heinz.py"
#' solver <- heinz.solver(heinz.py)
#' 
#' module.rn <- findModule(es.rn,
#'                         met.fdr=1e-6,
#'                         rxn.fdr=1e-6,
#'                         absent.met.score=-20,
#'                         solver=solver)
#' }
#' 
#' module.rn <- addMetabolitesForReactions(module.rn, es.rn)
#' module.rn <- addInterconnections(module.rn, es.rn)
#' module.rn <- addNormLogFC(module.rn)
#' module.rn <- removeHangingNodes(module.rn)
#' module.rn <- simplifyReactionNodes(module.rn, es.rn)
#' module.rn <- expandReactionNodeAttributesToEdges(module.rn)
#' @export
removeHangingNodes <- function(module) {
    res <- module
    res <- delete.vertices(res, V(res)[degree(res) == 1 & is.na(pval)])
    return(res)
    
}

#' Add reaction trans- edges connecting metabolites into a module
#' 
#' Adds edges corresponding to trans-RPAIRs into the module.
#' 
#' @param module Module to work with
#' @param es Experiment set object
#' @return Modified module
#' @examples
#' data("kegg.mouse.network")
#' data("examplesGAM")
#' library("mouseMacrophages")
#' data("mmpData")
#' 
#' heinz.py <- "/usr/local/lib/heinz/heinz.py"
#' solver <- heinz.solver(heinz.py)
#' 
#' es.re <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M1.M2,
#'                            gene.de=gene.de.M1.M2,
#'                            reactions.as.edges=TRUE)
#' 
#' \dontrun{
#' module.re <- findModule(es.re,
#'                         met.fdr=3e-05,
#'                         rxn.fdr=3e-05,
#'                         absent.met.score=-20,
#'                         solver=solver)
#' }
#' 
#' module.re <- addTransEdges(module.re, es.re)
#' @export
addTransEdges <- function(module, es) {
    edges.keep <- es$net.edges.ext.all$rptype %in% c("main", "trans")
    net.with.trans <- graph.from.tables(
        node.table=es$met.de.ext,
        edge.table=es$net.edges.ext.all[edges.keep,],
        node.col="ID",
        edge.cols=c("met.x", "met.y"),
        directed=FALSE)
    return(induced.subgraph(net.with.trans, V(module)$name))
}
