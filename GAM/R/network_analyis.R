#' @import BioNet 
NULL


# Preprocess experiment set's differential expression data for genes and metabolites
# 
# Converts metabolite and gene IDs to be the same as in network. Computes
# reaction differential expression data.
# 
# @param es Experiment set with DE data
# @param met.ids Type of IDs used in metabolite DE data (see met.id.map for possible values)
# @param gene.ids Type of IDs used in gene DE data (see colnames(es$gene.id.map) for possible values)
# @param plot If TRUE plot BUM fit
preprocessPvalAndMetDE <- function(es, met.ids, gene.ids, plot=T) {
    if (!is.null(es$gene.de)) {
        print("Processing gene p-values...")
        
        if (!is.null(es$gene.de$log2FC)) {
            es$gene.de$log2FC <- fixInf(es$gene.de$log2FC)
        }
        if (!is.null(gene.ids)) {
            gene.id.map <- es$network$gene.id.map
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
        if (!is.null(es$met.de$log2FC)) {
            es$met.de$log2FC <- fixInf(es$met.de$log2FC)        
        }
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

#' Make preprocessing necessary for \code{findModule} function
#' 
#' Converts gene and metabolite ids, makes an actual network to work with.
#' 
#' @param network Network object
#' @param met.de Differential expression data for metabolites
#' @param gene.de Differential expression data for genes
#' @param rxn.de Differential expression data for reactions (can be supplied instead of gene DE)
#' @param met.ids Type of IDs used in met.de, if NULL it will be determined automatically
#' @param gene.ids Type of IDs used in gene.de, if NULL it will be determined automatically
#' @param reactions.as.edges If TRUE, represent reaction as edges betwen metabolites,
#'                          otherwise represent them as nodes with connections to
#'                          compounds
#' @param collapse.reactions If TRUE, collapse reaction nodes if they share an enzyme
#'                           and at least one metabolite (only when reactions are nodes)
#' @param use.rpairs If TRUE, only rpairs will be used as reaction edges (only when reactions are edges)
#' @param plot If TRUE plot BUM-models
#' @importFrom plyr rename
#' @import data.table
#' @examples 
#' data(kegg.mouse.network)
#' data(examplesGAM)
#' es.re <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M1.M2,
#'                            gene.de=gene.de.M1.M2,
#'                            reactions.as.edges=TRUE)
#' @export
makeExperimentSet <- function(network, 
                              met.de=NULL, gene.de=NULL, rxn.de=NULL,
                              met.ids=NULL, gene.ids=NULL,
                              reactions.as.edges=TRUE,
                              collapse.reactions=TRUE,
                              use.rpairs=TRUE,
                              plot=TRUE) {
    if (is.null(met.ids) && !is.null(met.de)) {
        lazyData("met.id.map")
        met.id.map <- get("met.id.map")
        met.ids <- getIdType(met.de$ID, met.id.map)
    }
    
    if (is.null(gene.ids) && !is.null(gene.de)) {
        gene.ids <- getIdType(gene.de$ID, network$gene.id.map)
    }
    
    es <- newEmptyObject()
    es$network <- network
    es$met.de <- met.de
    es$gene.de <- gene.de
    es$rxn.de <- rxn.de 
    
    
    es$graph.raw <- es$network$graph.raw
    es$reactions.as.edges <- reactions.as.edges
    es$use.rpairs <- use.rpairs
    
    
    es <- preprocessPvalAndMetDE(es, met.ids=met.ids, gene.ids=gene.ids, plot=plot)
    
    if (!is.null(es$rxn.de)) {
        print("Processing reaction p-values...")
        if (!is.null(es$rxn.de$log2FC)) {
            es$rxn.de$log2FC <- fixInf(es$rxn.de$log2FC)                        
        }
        if ("origin" %in% colnames(es$rxn.de)) {
            
            es$rxn.de <- es$rxn.de[es$rxn.de$ID %in% es$graph.raw$rxn, ]            
            
            gene.id.map <- es$network$gene.id.map
            
            unknown.rxn <- setdiff(es$rxn.de$ID, es$network$rxn2name$rxn)
            unknown.rxn2name <- do.call(cbind,  c(
                list(unknown.rxn), 
                rep(list(rep("", length(unknown.rxn))), ncol(es$network$rxn2name) - 1)))
            unknown.rxn2name <- as.data.frame(unknown.rxn2name, stringsAsFactors=F)
            colnames(unknown.rxn2name) <- colnames(es$network$rxn2name)
            
            es$network$rxn2name <- rbind(es$network$rxn2name, unknown.rxn2name)
            
            es$network$rxn2name$name[match(es$rxn.de$ID, es$network$rxn2name$rxn)] <-
                gene.id.map$Symbol[match(es$rxn.de$origin, gene.id.map[,es$network$gene.ids])]            
            
            
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
            log2FC=NA,
            stringsAsFactors=F)
    }
    
    met.de.ext <- data.frame(ID=unique(c(es$graph.raw$met.x, es$graph.raw$met.y)), stringsAsFactors=F)
    met.de.ext <- merge(met.de.ext, es$met.de, all.x=T) # all.x=T — we keep mets if there is no MS data
    met.de.ext <- merge(met.de.ext, es$network$met2name, by.x="ID", by.y="met", all.x=T)
    met.de.ext$logPval <- log(met.de.ext$pval)
    es$met.de.ext <- met.de.ext
    
    
    if (is.null(es$rxn.de)) {
        es$rxn.de <- data.frame(ID=es$graph.raw$rxn, pval=NA, log2FC=NA, stringsAsFactors=F)
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
            es$rxn.de.origin.split <- splitMappingByConnectivity(rxn.net, rxn.de.ext$ID, rxn.de.ext$origin)
            
            t <- data.frame(from=es$rxn.de.ext$ID, to=es$rxn.de.origin.split, stringsAsFactors=F)                
            #from.table <- aggregate(from ~ to, t, function(x) paste(x, collapse="+"))
            
            es$rxn.de.ext$rxns <- es$rxn.de.ext$ID
            es$rxn.de.ext$ID <- es$rxn.de.origin.split            
            
            uniqueJoin <- function(x) {
                u  <-  unique(unlist(strsplit(as.character(x), split="\\+")))
                return(paste(u, collapse="+"))
            }
            
            cols.to.fix <- setdiff(c(colnames(es$network$rxn2name), "rxns"), c("rxn", "name"))
            f <- as.formula(paste0("cbind(", paste0(cols.to.fix, collapse=","), ") ~ ID"))
            cols.fixed <- aggregate(f , es$rxn.de.ext, uniqueJoin)
            
            es$rxn.de.ext <- unique(es$rxn.de.ext[, !names(es$rxn.de.ext) %in% cols.to.fix])
            es$rxn.de.ext <- merge(es$rxn.de.ext, cols.fixed)
            
            edges$rxn <- es$rxn.de.origin.split[edges$rxn]
            edges <- unique(edges)
        }
        
        es$edges <- edges
        net1 <- graph.from.tables(node.table=list(met=es$met.de.ext, rxn=es$rxn.de.ext), edge.table=edges,
                                     node.col="ID",
                                     directed=F)
        
        es$subnet <- net1
    }
    
    
    return(es)
}



# Scores p-value by fitted BUM-model
# This is a helper function based on BioNet::scoreFunction()
scoreValue <- function (fb, pval, fdr = 0.01) 
{
    return((fb$a - 1) * (log(pval) - log(fdrThreshold(fdr, fb))))
}

#' Assign scores to network's nodes and edges
#' @param es Experiment set object
#' @param fdr FDR for both metabolites and genes/reactions, if set met.fdr and gene.fdr aren't used
#' @param met.fdr FDR for metabolites only
#' @param gene.fdr FDR for genes/reactions only
#' @param absent.met.score Score for metabolites absent from data
#' @param absent.rxn.score Score for reactions when there is no genomic data
#' @return Experiment set object with scored network subnet.scored field
#' @import igraph 
#' @examples
#' data(kegg.mouse.network)
#' data(examplesGAM)
#' es.re <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M1.M2,
#'                            gene.de=gene.de.M1.M2,
#'                            reactions.as.edges=TRUE)
#' es.re.scored <- scoreNetwork(es.re, met.fdr=3e-5, gene.fdr=3e-5, absent.met.score=-20)
#' @export 
scoreNetwork <- function(es,                         
                       fdr=NULL, met.fdr=NULL, gene.fdr=NULL,
                       absent.met.score=NULL,
                       absent.rxn.score=0) {
    
    net <- es$subnet
    
    if (!is.null(fdr)) {
        met.fdr <- fdr
        gene.fdr <- fdr
    }
    
            
    met.scores <- NULL    
    if (!is.null(es$fb.met) && !is.null(met.fdr)) {            
        fb <- es$fb.met
        met.scores <- scoreValue(fb, na.omit(es$met.de.ext$pval), met.fdr)
        names(met.scores) <- es$met.de.ext$ID[!is.na(es$met.de.ext$pval)]
    }
    
    
    rxn.scores <- NULL
    if (!is.null(es$fb.rxn) && !is.null(gene.fdr)) {
        fb <- es$fb.rxn 
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
    
    
    
    es$subnet.scored <- net
    return(es)
}


#' Assign scores to network's nodes and edges without using BUM-model
#' Score for node or edge is log(pval)-log(pval.threshold)
#' @param es Experiment set object
#' @param met.pval.threshold Threshold so metabolic p-values
#' @param met.pval.default P-value for metabolites without provided p-value
#' @param rxn.pval.threshold Threshold so reaction p-values
#' @param rxn.pval.default P-value for reactions without provided p-value
#' @return Experiment set object with scored network subnet.scored field
#' @import igraph 
#' @examples
#' data(kegg.mouse.network)
#' data(examplesGAM)
#' es.re <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M1.M2,
#'                            gene.de=gene.de.M1.M2,
#'                            reactions.as.edges=TRUE)
#' es.re.scored <- scoreNetworkWithoutBUM(es.re,
#'                                        met.pval.threshold=1e-5,
#'                                        met.pval.default=1,
#'                                        rxn.pval.threshold=1e-5)
#' @export 
scoreNetworkWithoutBUM <- function(es,
                                   met.pval.threshold=1e-5,
                                   met.pval.default=NULL,
                                   rxn.pval.threshold=1e-5,
                                   rxn.pval.default=NULL
                                   ) {
    
    net <- es$subnet
    
            
    met.scores <- NULL    
    if (!is.null(es$fb.met) && !is.null(met.pval.threshold)) {            
        met.scores <- -log(es$met.de.ext$pval) + log(met.pval.threshold)
        names(met.scores) <- es$met.de.ext$ID
        met.scores <- na.omit(met.scores)
    }
    
    
    rxn.scores <- NULL
    if (!is.null(es$fb.rxn) && !is.null(rxn.pval.threshold)) {
        rxn.scores <- -log(es$rxn.de.ext$pval) + log(rxn.pval.threshold)
        names(rxn.scores) <- es$rxn.de.ext$ID
        rxn.scores <- na.omit(rxn.scores)
    }
    
    
    if (is.null(met.pval.default)) {
        absent.met.score <- mean(met.scores)
    } else {
        absent.met.score <- -log(met.pval.default) + log(met.pval.threshold)
    }
        
    absent.met.scores <- sapply(V(net)[nodeType == "met"]$name, function(x) absent.met.score)
    met.scores <- c(met.scores, absent.met.scores[!names(absent.met.scores) %in% names(met.scores)])
    
    met.scores <- met.scores[names(met.scores) %in% V(net)$name]
    V(net)[names(met.scores)]$score <- met.scores
    
    if (is.null(rxn.pval.default)) {
        absent.rxn.score <- mean(rxn.scores)
    } else {
        absent.rxn.score <- -log(rxn.pval.default) + log(rxn.pval.threshold)
    }
    
    absent.rxn.scores <- sapply(es$rxn.de$ID, function(x) absent.rxn.score)
    rxn.scores <- c(rxn.scores, absent.rxn.scores[!names(absent.rxn.scores) %in% names(rxn.scores)])
    
    if (es$reactions.as.edges) {
        
        E(net)$score <- rxn.scores[E(net)$rxn]
    } else {
        rxn.scores <- rxn.scores[names(rxn.scores) %in% V(net)$name]
        V(net)[names(rxn.scores)]$score <- rxn.scores
    }
    
    
    
    es$subnet.scored <- net
    return(es)
}


#' Find significant module in a network
#' @param es Experiment set object (scored or not)
#' @param solver Solver function of MWCS problem to use, first argument should be a network,
#'               result should be a module or list of modules
#' @param score.function Function to score network (is applied only if the network wasn't scored previously)
#' @param simplify If TRUE and only one module was found return just the module, not a list.
#' @param ... Additional arguments for scoring function
#' @return List of most significant modules
#' @import igraph 
#' @examples
#' data(kegg.mouse.network)
#' data(examplesGAM)
#' es.re <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M1.M2,
#'                            gene.de=gene.de.M1.M2,
#'                            reactions.as.edges=TRUE)
#' solver <- heinz.solver("/usr/local/lib/heinz/heinz.py")
#' \dontrun{
#' module.re <- findModule(es.re, solver, met.fdr=3e-5, gene.fdr=3e-5, absent.met.score=-20)
#' }
#' @export 
findModule <- function(es,                         
                       solver = fastHeinz.solver,
                       simplify=TRUE,
                       score.function=scoreNetwork,
                       ...) {
    
    if (is.null(es$subnet.scored)) {
        es <- score.function(es, ...)
    }
    
    res <- solver(es$subnet.scored)
    
    if (simplify && is(res, "list") && length(res) == 1) {
        return(res[[1]])
    }
    
    return(res)
}
