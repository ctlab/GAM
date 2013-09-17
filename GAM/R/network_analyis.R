#' @import BioNet igraph0
NULL

#' Plot module with attribute as a node label
#' @param module Module to plot
#' @param attr Attribute to use as a label
#' @param layout Layout to use
#' @export
plotNetwork <- function(module, attr.label="label", layout=layout.kamada.kawai, ...) {
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
    plotModule(module, diff.expr=de, layout=layout, labels=all.labels[node.labels], ...)        
}


saveModuleToPdf <- function(module, outputFilePrefix) {
    pdf(paste(outputFilePrefix, "pdf", sep="."), width=15, height=15)
    plotNetwork(module, attr.label="label", vertex.size=2)
    dev.off()
}

#' Save module to different formats
#' @param module Module to save
#' @param outputFilePrefix Path to save to (without extension)
#' @param types Vector of file types, "pdf" or one of the supported by BioNet::saveNetwork function
#' @export
saveModule <- function(module, outputFilePrefix, types=c("pdf", "XGMML")) {
    # :ToDO: fix saving to XGMML (NAs, escaping, trivial modules)
    outdir <- dirname(outputFilePrefix)
    
    if (!file.exists(outdir)) {
        dir.create(outdir)
    }
    
    
    #module <- subNetwork(nodes(module), network=subnet)
    module    

    
    for (type in types) {
        if (type == "pdf") {
            saveModuleToPdf(module, outputFilePrefix)
        } else {
            saveNetwork(module,name=basename(outputFilePrefix),file=outputFilePrefix,type=type)
        }
        
    }
    
    
}

#' Add node attribute with normalized log-foldchange inside node groups
#' @param module Module to add attribut to
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
preprocessPvalAndMetDE <- function(es, met.ids, gene.ids) {
    if (!is.null(es$gene.de)) {
        print("Processing gene p-values...")
        es$gene.de$logFC <- fixInf(es$gene.de$logFC)
        if (!is.null(gene.ids)) {
            data("gene.id.map")
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
            data("met.id.map")
            es$met.de <- convertPval(es$met.de, 
                                      from=met.id.map[,met.ids], 
                                      to=met.id.map[,es$network$met.ids])
        }
        
        es$met.pval <- es$met.de$pval
        names(es$met.pval) <- es$met.de$ID
        es$fb.met <- fitBumModel(es$met.pval, plot = TRUE)
    } else {
        es$met.pval <- NULL
    }
    
    return(es)
}

#' @export
makeExperimentSet <- function(network, 
                              met.de=NULL, gene.de=NULL, rxn.de=NULL,
                              met.ids=NULL, gene.ids=NULL,
                              reactions.as.edges=F
) {
    es <- newEmptyObject()
    es$network <- network
    es$met.de <- met.de
    es$gene.de <- gene.de
    es$rxn.de <- rxn.de 
    
    
    es$graph.raw <- es$network$graph.raw
    es$reactions.as.edges <- reactions.as.edges
    
    
    es <- preprocessPvalAndMetDE(es, met.ids=met.ids, gene.ids=gene.ids)
    
    if (!is.null(es$rxn.de)) {
        print("Processing reaction p-values...")
        es$rxn.de$logFC <- fixInf(es$rxn.de$logFC)                        
        if ("origin" %in% colnames(es$rxn.de)) {
            
            es$rxn.de <- es$rxn.de[es$rxn.de$ID %in% es$graph.raw$rxn, ]            
            
            matched <- es$rxn.de$ID %in% es$network$rxn2name$rxn
            
            es$network$rxn2name$name[match(es$rxn.de$ID[matched], es$network$rxn2name$rxn)] <-
                gene.id.map$name[match(es$rxn.de$origin[matched], gene.id.map[,es$network$gene.ids])]            
            
            es$network$rxn2name <- rbind(es$network$rxn2name, 
                                         cbind(rxn=es$rxn.de$ID[!matched],
                                               name=gene.id.map$name[match(es$rxn.de$origin[!matched], gene.id.map[,es$network$gene.ids])]
                                               ))
            
        }
        
        es$rxn.pval <- es$rxn.de$pval
        names(es$rxn.pval) <- es$rxn.de$ID
        
        es$fb.rxn <- fitBumModel(es$rxn.pval, plot = TRUE)
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
        
        
        net.edges.pval <- (aggregate(pval ~ met.x * met.y, data=net.edges.ext, min))
        net.edges.ext <- merge(net.edges.pval, net.edges.ext)
        net.edges.ext <- net.edges.ext[!duplicated(net.edges.ext[,c("met.x", "met.y")]),]        
        net.edges.ext <- net.edges.ext[net.edges.ext$met.x != net.edges.ext$met.y,]    
        
        es$net.edges.ext <- net.edges.ext
        
        
        net1 <- graphNEL.from.tables(node.table=met.de.ext, edge.table=net.edges.ext,
                                     node.col="ID", edge.cols=c("met.x", "met.y"),
                                     directed=F, ignore.solitary.nodes=T)
        
        es$subnet <- net1
    } else {    
        
        edges <- rbind(rename(es$graph.raw[, c("met.x", "rxn")], c("met.x" = "met")),
                       rename(es$graph.raw[, c("met.y", "rxn")], c("met.y" = "met")))
        
        net1 <- graphNEL.from.tables(node.table=list(met=met.de.ext, rxn=rxn.de.ext), edge.table=edges,
                                     node.col="ID",
                                     directed=F, ignore.solitary.nodes=T)
        
        es$subnet <- net1
    }
    
    

    
    
    
    
    es$met.de$origin <- NULL
    es$rxn.de$origin <- NULL
    es$all.de <- rbind(es$met.de, es$rxn.de)
    es$all.pval <- c(es$rxn.pval, es$met.pval)
    
    es$fb.all <- fitBumModel(es$all.pval, plot = TRUE)
    return(es)
}

appendModule <- function(res, module.graph) {        
    res[[length(res)+1]] <- module.graph
    res
}


#' @export
findModules <- function(es,                         
                         fdr=NULL, met.fdr=NULL, gene.fdr=NULL,
                         absent.met.score=NULL,
                         score.separately=T,
                         heinz.py, heinz.nModules=1,
                         heinz.tolerance=10,
                         heinz.subopt_diff=100) {
    
    if (!is.null(fdr)) {
        met.fdr <- fdr
        gene.fdr <- fdr
    }
    
            
    met.scores <- NULL    
    if (!is.null(es$met.de) && !is.null(met.fdr)) {            
        if (score.separately) {
            met.scores <- scoreFunction(es$fb.met, met.fdr)             
        } else {
            met.scores <- scoreFunction(es$fb.all, met.fdr)[names(es$met.pval)]                
        }
    }
    
    
    rxn.scores <- NULL
    if (!is.null(es$rxn.de) && !is.null(gene.fdr)) {
        if (score.separately) {             
            rxn.scores <- scoreFunction(es$fb.rxn, gene.fdr)
        } else {                
            rxn.scores <- scoreFunction(es$fb.all, gene.fdr)[names(es$rxn.pval)]
        }
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
    
    
    res <- runHeinz(
        subnet=es$subnet, 
        heinz.py=heinz.py, 
        score.edges="score" %in% names(edgeDataDefaults(es$subnet)),
        score.nodes="score" %in% names(nodeDataDefaults(es$subnet)),
        heinz.nModules=heinz.nModules, 
        heinz.tolerance=heinz.tolerance,
        heinz.subopt_diff=heinz.subopt_diff)        
        
    return(res)
}



runHeinz <- function(subnet,
                      heinz.py, 
                      score.edges=F,
                      score.nodes=T,                      
                      heinz.nModules=1, 
                      heinz.tolerance=10,
                      heinz.subopt_diff=100) {
    tmpdir <- tempdir()        
    tmpdir <- "/tmp"
    edges_file <- paste(tmpdir, "edges.txt", sep="/")
    nodes_file <- paste(tmpdir, "nodes.txt", sep="/")
    
    writeHeinzEdges(subnet, file=edges_file, use.score=score.edges)
    
    if (!score.nodes) {
        # Hack to make writeHeinzeNodes working
        nodeDataDefaults(subnet, "score")  <- 0
        nodeData(subnet, attr="score") <- 0
    }
    writeHeinzNodes(subnet, file=nodes_file, use.score=T)
    
    
    system2(heinz.py,
            c("-n", nodes_file,
              "-e", edges_file,
              "-N", if (score.nodes) "True" else "False",
              "-E", if (score.edges) "True" else "False",              
              "-s", heinz.nModules,
              "--tolerance", heinz.tolerance,
              "--subopt_diff", heinz.subopt_diff,
              "-v"));
    
    
    res <- list()
    for (i in 0:(heinz.nModules-1)) {
        module.graph = readHeinzGraph(node.file = paste(nodes_file, i, "hnz", sep="."), 
                                      network = subnet)
        res <- appendModule(res, module.graph)
        
    }
    return(res)
}

makeExperimentSetSq <- function(network, 
                                met.de=NULL, gene.de=NULL, rxn.de=NULL,
                                met.ids=NULL, gene.ids=NULL,
                                collapse.reactions=F
) {
    es <- newEmptyObject()
    es$network <- network
    es$met.de <- met.de
    es$gene.de <- gene.de
    es$rxn.de <- rxn.de 
    
    
    es$graph <- es$network$graph.sq
    
    es <- preprocessPvalAndMetDE(es, met.ids=met.ids, gene.ids=gene.ids)
    
    if (!is.null(es$rxn.de)) {
        print("Processing reaction p-values...")
        es$rxn.de$logFC <- fixInf(es$rxn.de$logFC)                        
        if ("origin" %in% colnames(es$rxn.de)) {
            if (collaps.reactions) {
                print("Collapsing reactions by common most significant enzymes")
                #rxn.de.origin.split <- es$rxn.de$origin
                es$rxn.de.origin.split <- split.mapping.by.connectivity(es$graph, es$rxn.de$ID, es$rxn.de$origin)
                
                es$graph <- convert.node.names(es$graph, es$rxn.de$ID, es$rxn.de.origin.split)
                
                
                es$rxn.de.orig <- es$rxn.de
                es$rxn.de$ID <- es$rxn.de.origin.split            
            }
            
            es$rxn.de <- es$rxn.de[es$rxn.de$ID %in% nodes(es$graph), ]            
            nodeData(es$graph, es$rxn.de$ID, "shortName") <- 
                gene.id.map$name[match(es$rxn.de$origin, gene.id.map[,es$network$gene.ids])]
            nodeData(es$graph, es$rxn.de$ID, "nodeType") <- "rxns"                
            es$rxn.de$origin <- NULL
        }
        
        es$rxn.pval <- es$rxn.de$pval
        names(es$rxn.pval) <- es$rxn.de$ID
        
        es$fb.rxn <- fitBumModel(es$rxn.pval, plot = TRUE)
    } else {
        es$rxn.pval <- NULL
    }
    es$met.de$origin <- NULL
    es$rxn.de$origin <- NULL
    es$all.de <- rbind(es$met.de, es$rxn.de)
    es$all.pval <- c(es$rxn.pval, es$met.pval)
    
    print("Processing all p-values together")    
    es$subnet <- subNetwork(es$all.de$ID, es$graph)
    
    # :ToDo: change Inf value for genes and metabolites
    dm <- es$all.de[, "logFC"]   
    
    
    names(dm) <- es$all.de$ID
    nodeDataDefaults(es$subnet, "logFC") <- 0
    nodeData(es$subnet, nodes(es$subnet), "logFC") <- dm[nodes(es$subnet)]
    
    
    es$fb.all <- fitBumModel(es$all.pval, plot = TRUE)
    return(es)
}

findModulesSq <- function(es,                         
                            fdr=NULL,
                            met.fdr=NULL, gene.fdr=NULL,
                            score.separately=F,
                            heinz.py=NULL, heinz.nModules=1,
                            heinz.tolerance=10,
                            heinz.subopt_diff=100) {
    
    if (!is.null(fdr)) {
        met.fdr <- fdr
        gene.fdr <- fdr
    }
    
    met.scores <- NULL
    if (!is.null(es$met.de)) {            
        
        print(paste("met.fdr =", met.fdr))
        if (score.separately) {
            met.scores <- scoreFunction(es$fb.met, met.fdr)             
        } else {
            met.scores <- scoreFunction(es$fb.all, met.fdr)[names(es$met.pval)]                
        }
    }
    
    rxn.scores <- NULL
    if (!is.null(es$rxn.de)) {            
        print(paste("gene.fdr =", gene.fdr))
        if (score.separately) {             
            rxn.scores <- scoreFunction(es$fb.rxn, gene.fdr)
        } else {                
            rxn.scores <- scoreFunction(es$fb.all, gene.fdr)[names(es$rxn.pval)]
        }
    }
    
    all.scores <- c(met.scores, rxn.scores)        
    scores <- all.scores[nodes(es$subnet)]
    
    nodeDataDefaults(es$subnet, "score") <- -Inf
    nodeData(es$subnet, nodes(es$subnet), "score") <- scores
    
    if (!is.null(heinz.py)) {
        res <- runHeinz(
            subnet=es$subnet, 
            heinz.py=heinz.py, 
            score.edges=F, 
            score.nodes=T,
            heinz.nModules=heinz.nModules, 
            heinz.tolerance=heinz.tolerance,
            heinz.subopt_diff=heinz.subopt_diff)        
    } else {
        res <- runFastHeinz(es$subnet, scores)
        
    }
    res <- lapply(res, addNormLogFC)
    
    
    return(res)
}
