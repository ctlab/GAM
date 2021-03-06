#' @import BioNet 
NULL

# Preprocess experiment set's differential expression data for genes
# 
# Converts gene IDs to be the same as in network
# 
# @param es Experiment set with gene DE data
# @param gene.ids Type of IDs used in gene DE data (see colnames(es$gene.id.map) for possible values)
# @param plot If TRUE plot BUM fit
preprocessGeneDE <- function(es, gene.ids) {
    if (is.null(es$gene.de)) {
        return(es)
    }
    message("Preprocessing gene p-values...")
    
    if (!is.null(es$gene.de$log2FC)) {
        es$gene.de$log2FC <- fixInf(es$gene.de$log2FC)
    }
    if (!is.null(gene.ids)) {
        gene.id.map <- es$network$gene.id.map
        if (!is(gene.id.map, "data.table")) {
            gene.id.map <- as.data.table(gene.id.map)
        }
        message(sprintf("Converting gene ids from %s to %s", gene.ids, es$network$gene.ids))
        es$gene.de <- convertPval(es$gene.de, 
                                  from=gene.id.map[,get(gene.ids)], 
                                  to=gene.id.map[,get(es$network$gene.ids)])
    }    
    es$gene.de$origin <- NULL
    message("Converting gene p-values to reactions...")
    es$rxn.de <- convertPval(es$gene.de, from=es$network$rxn2gene$gene, to=es$network$rxn2gene$rxn)
        
    es
}

# Preprocess experiment set's differential expression data for genes and metabolites
# 
# Converts metabolite and gene IDs to be the same as in network. Computes
# reaction differential expression data.
# 
# @param es Experiment set with DE data
# @param met.ids Type of IDs used in metabolite DE data (see met.id.map for possible values)
# @param gene.ids Type of IDs used in gene DE data (see colnames(es$gene.id.map) for possible values)
# @param plot If TRUE plot BUM fit
preprocessMetDE <- function(es, met.ids, plot=TRUE) {    
    if (is.null(es$met.de)) {
        es$met.de <- data.frame(
            ID=unique(c(es$graph.raw$met.x, es$graph.raw$met.y)), 
            pval=NA,
            log2FC=NA,
            stringsAsFactors=FALSE)
        
        es$met.pval <- NULL        
        return(es)
    }
    
    message("Processing metabolite p-values...")
    if (!is.null(es$met.de$log2FC)) {
        es$met.de$log2FC <- fixInf(es$met.de$log2FC)        
    }
    if (!is.null(met.ids)) {
        lazyData("met.id.map")
        met.id.map <- get("met.id.map")
        if (!is(met.id.map, "data.table")) {
            met.id.map <- as.data.table(met.id.map)
        }
        message(sprintf("Converting metabolite ids from %s to %s", met.ids, es$network$met.ids))
        es$met.de <- convertPval(es$met.de, 
                                 from=met.id.map[, get(met.ids)], 
                                 to=met.id.map[, get(es$network$met.ids)])
    }
    
    es$met.pval <- es$met.de$pval
    names(es$met.pval) <- es$met.de$ID
    es$fb.met <- fitBumModel(es$met.pval[!duplicated(es$met.pval)], plot=FALSE)
    if (plot) {
        hist(es$fb.met, main="Histogram of metabolite p-values")
        plot(es$fb.met, main="QQ-Plot for metabolite BUM-model")
    }
    
    es
}

#' Process reaction DE
#' 
#' Sets reaction names, fits BUM-model.
#' 
#' @param es Experiment set
#' @param plot If TRUE plot BUM-models
#' @return Experiment set object whith updated network and fit BUM-model
#' @import data.table
processReactionDE <- function(es, plot) {
    if (is.null(es$rxn.de)) {
        es$rxn.de <- data.frame(ID=unique(es$graph.raw$rxn), pval=NA, log2FC=NA, stringsAsFactors=FALSE)
        es$rxn.pval <- NULL
        return(es)
    }
    
    message("Processing reaction p-values...")
    if (!is.null(es$rxn.de$log2FC)) {
        es$rxn.de$log2FC <- fixInf(es$rxn.de$log2FC)                        
    }
    
    if ("symbol" %in% colnames(es$rxn.de)) {
        symbol <- with(es$rxn.de, { x <- symbol; names(x) <- ID; x[intersect(ID, es$network$rxn2name$rxn)] } )
        es$network$rxn2name$name[match(names(symbol), es$network$rxn2name$rxn)] <- symbol
    } else if ("origin" %in% colnames(es$rxn.de)) {
        
        es$rxn.de <- es$rxn.de[es$rxn.de$ID %in% es$graph.raw$rxn, ]            
        
        gene.id.map <- es$network$gene.id.map
        
        unknown.rxn <- setdiff(es$rxn.de$ID, es$network$rxn2name$rxn)
        unknown.rxn2name <- do.call(cbind,  c(
            list(unknown.rxn), 
            rep(list(rep("", length(unknown.rxn))), ncol(es$network$rxn2name) - 1)))
        unknown.rxn2name <- as.data.frame(unknown.rxn2name, stringsAsFactors=FALSE)
        colnames(unknown.rxn2name) <- colnames(es$network$rxn2name)
        
        es$network$rxn2name <- rbind(es$network$rxn2name, unknown.rxn2name)
        
        es$network$rxn2name$name[match(es$rxn.de$ID, es$network$rxn2name$rxn)] <-
            gene.id.map$Symbol[match(es$rxn.de$origin, gene.id.map[[es$network$gene.ids]])]
        
        
    }
    
    es$gene.pval <- es$rxn.de$pval[!duplicated(es$rxn.de$origin)]
    names(es$gene.pval) <- unique(es$rxn.de$origin)
    
    es$fb.gene <- fitBumModel(es$gene.pval, plot=FALSE)
    if (plot) {
        hist(es$fb.gene, main="Histogram of gene p-values")
        plot(es$fb.gene, main="QQ-Plot for gene BUM-model")
    }            
    
    es
}

makeSubnetWithAtoms <- function(es) {
    lazyData("kegg.db")
    es$rxn.de.ext$rxn <- es$rxn.de.ext$ID
    rpair.de <- convertPval(es$rxn.de.ext, kegg.db$rpairs$rxn, kegg.db$rpairs$rpair)
    
    edges <- merge(kegg.db$rpaligns, rpair.de, by.x="rpair", by.y="ID")
    #edges
    #edges[atom.x > atom.y,] list(atom.x, atom.y)]  <- edges[atom.x > atom.y, list(atom.y, atom.x)]
    edges <- edges[!duplicated(edges[, c("atom.x", "atom.y")]), ]
    edges <- edges[edges$atom.x != edges$atom.y, ]
    edges$label <- edges$symbol
    
    nodes <- data.frame(ID=with(kegg.db$rpaligns, unique(c(atom.x, atom.y))), stringsAsFactors=F)
    nodes$KEGG <- gsub("_.*$", "", nodes$ID)
    
    nodes <- merge(nodes, es$met.de.ext, by.x="KEGG", by.y="ID", all.x=T, suffixes=c("", ".orig"))
    
    #nodes <- rename(nodes, c("name"="label"))
    
    
    atom.graph <- GAM:::graph.from.tables(node.table=list(met=nodes), node.col="ID", name.as.label=T,
                                          edge.table=edges, 
                                          edge.cols=c("atom.x", "atom.y"), directed=F)
    
    
    res <- delete.vertices(atom.graph, V(atom.graph)[KEGG %in% kegg.db$mets2mask]$name)
    gc <- clusters(res)
    res <- induced.subgraph(res, gc$membership == which.max(gc$csize))
    
    res 
    es$subnet <- res
    es    
}


makeSubnetWithReactionsAsEdges <- function(es) {
    net.edges.ext <- merge(es$graph.raw, es$rxn.de.ext, by.x="rxn", by.y="ID")
    edges2rev <- net.edges.ext$met.x > net.edges.ext$met.y
    net.edges.ext[edges2rev, c("met.x", "met.y")] <- 
        net.edges.ext[edges2rev, c("met.y", "met.x")]        
    
    
    net.edges.pval <- (aggregate(pval ~ met.x + met.y + rptype, data=net.edges.ext, min, na.action=na.pass))
    net.edges.ext <- merge(net.edges.pval, net.edges.ext)
    net.edges.ext <- net.edges.ext[
        order(
            factor(
                net.edges.ext$rptype, 
                levels=unique(c(c("main", "trans"), net.edges.ext$rptype)))),
                ]
    net.edges.ext <- net.edges.ext[!duplicated(net.edges.ext[,c("met.x", "met.y")]),]        
    net.edges.ext <- net.edges.ext[net.edges.ext$met.x != net.edges.ext$met.y,]    
    
        
    if (es$use.rpairs) {
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
                              directed=FALSE)
    
    es$subnet <- net1
    es    
    
}

makeSubnetWithReactionsAsNodes <- function(es) {
    edges <- rbind(rename(es$graph.raw[, c("met.x", "rxn")], c("met.x" = "met")),
                   rename(es$graph.raw[, c("met.y", "rxn")], c("met.y" = "met")))
    edges <- unique(edges)
    
    if (es$collapse.reactions && "origin" %in% names(es$rxn.de.ext)) {
        message("Collapsing reactions by common most significant enzymes")
        
        edges.dt <- data.table(edges, key="met")
        rxn.net <- merge(edges.dt, edges.dt, by="met", allow.cartesian=TRUE)
        rxn.net <- unique(rxn.net[, list(rxn.x, rxn.y)])
        
        #rxn.de.origin.split <- es$rxn.de$origin
        es$rxn.de.origin.split <- splitMappingByConnectivity(rxn.net, es$rxn.de.ext$ID, es$rxn.de.ext$origin)
        
        t <- data.frame(from=es$rxn.de.ext$ID, to=es$rxn.de.origin.split, stringsAsFactors=FALSE)                
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
                              directed=FALSE)
    
    es$subnet <- net1
    es
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
#' @param graph.typ type of the network :TODO:
#' @param collapse.reactions If TRUE, collapse reaction nodes if they share an enzyme
#'                           and at least one metabolite (only when reactions are nodes)
#' @param use.rpairs If TRUE, only rpairs will be used as reaction edges (only when reactions are edges)
#' @param drop.nonexpressed If TRUE, reactions that doesn't have expressed enzymes are removed from the network
#' @param plot If TRUE plot BUM-models
#' @importFrom plyr rename
#' @import data.table
#' @examples 
#' library("GAM.db")
#' library("GAM.networks")
#' data(kegg.mouse.network)
#' data(examplesGAM)
#' es.re <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M0.M1,
#'                            gene.de=gene.de.M0.M1,
#'                            graph.type="reactions_as_edges")
#' @export
makeExperimentSet <- function(network, 
                              met.de=NULL, gene.de=NULL,
                              met.ids=NULL, gene.ids=NULL,
                              graph.type=c("atom", "reactions_as_edges", "reactions_as_nodes"),
                              collapse.reactions=TRUE,
                              use.rpairs=TRUE,
                              drop.nonexpressed=TRUE,
                              plot=TRUE) {
    graph.type <- match.arg(graph.type)
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
  
    
    
    es$graph.raw <- es$network$graph.raw
    es$graph.type <- graph.type
    es$use.rpairs <- use.rpairs
    es$collapse.reactions <- collapse.reactions
    es$drop.nonexpressed <- drop.nonexpressed
    
    es <- preprocessGeneDE(es, gene.ids=gene.ids)    
    es <- preprocessMetDE(es, met.ids=met.ids, plot=plot)
    es <- processReactionDE(es, plot)
    
    
        
    message("Building network")    

    lazyData("kegg.db")    
    allMets <- unique(c(es$graph.raw$met.x, 
                        es$graph.raw$met.y,
                        kegg.db$rpairs$met.x,
                        kegg.db$rpairs$met.y))
    met.de.ext <- data.frame(ID=allMets, stringsAsFactors=FALSE)
    met.de.ext <- merge(met.de.ext, es$met.de, all.x=TRUE) # all.x=TRUE — we keep mets if there is no MS data
    met.de.ext <- merge(met.de.ext, es$network$met2name, by.x="ID", by.y="met", all.x=TRUE)
    met.de.ext$logPval <- log(met.de.ext$pval)
    if (!"ion" %in% colnames(met.de.ext)) {
        met.de.ext$ion <- sprintf("%s_%.2f", met.de.ext$formula, -met.de.ext$logPval)
    }
    
    es$met.de.ext <- met.de.ext
    
    
    rxn.de.ext <- data.frame(ID=unique(es$graph.raw$rxn), stringsAsFactors=FALSE)
    rxn.de.ext <- merge(rxn.de.ext, es$rxn.de, all.x=!es$drop.nonexpressed)
    rxn.de.ext <- merge(rxn.de.ext, es$network$rxn2name, by.x="ID", by.y="rxn", all.x=!es$drop.nonexpressed)
    rxn.de.ext$logPval <- log(rxn.de.ext$pval)
    es$rxn.de.ext <- rxn.de.ext
    
    es$graph.raw <- es$graph.raw[es$graph.raw$rxn %in% rxn.de.ext$ID,]
    es$graph.raw <- es$graph.raw[es$graph.raw$met.x %in% met.de.ext$ID,]
    es$graph.raw <- es$graph.raw[es$graph.raw$met.y %in% met.de.ext$ID,]
    
    
    if (graph.type == "reactions_as_edges") {    
        es <- makeSubnetWithReactionsAsEdges(es)
    } else  if (graph.type == "reactions_as_nodes") {    
        es <- makeSubnetWithReactionsAsNodes(es)
    } else if (graph.type == "atom") {
        es <- makeSubnetWithAtoms(es)
    } else {
        stop(sprintf("Unknown graph type: %s", graph.type))
    }
    
    
    return(es)
}



# Scores p-value by fitted BUM-model
# This is a helper function based on BioNet::scoreFunction()
#' @import BioNet
scoreValue <- function (fb, pval, fdr = 0.01) 
{
    return((fb$a - 1) * (log(pval) - log(fdrThreshold(fdr, fb))))
}

recommendedFDR <- function(fb, pval, num.positive=100) {
    lo <- -1    
    hi <- -200
    
    for (i in 1:10) {
        mid <- (lo + hi) / 2
        val <- sum(scoreValue(fb, pval, exp(mid)) > 0)
        if (val < num.positive) {
            # we want higher FDR
            hi <- mid
        } else {
            lo <- mid
        }
    }
    exp(mid)    
}

#' Assign scores to network's nodes and edges
#' @param es Experiment set object
#' @param fdr FDR for both metabolites and genes/reactions, if set met.fdr and rxn.fdr aren't used
#' @param met.fdr FDR for metabolites only, NA makes FDR to be set automatically, NULL makes metabolites not be scored
#' @param rxn.fdr FDR for genes/reactions only, NA makes FDR to be set automatically, NULL makes reactions not be scored
#' @param absent.met.score Score for metabolites absent from data
#' @param absent.rxn.score Score for reactions when there is no genomic data
#' @param met.score Score for metabolites if met.fdr is NULL or there is no metabolic data
#' @param rxn.score Score for reactions if met.fdr is NULL or there is no transcriptional data
#' @param num.positive Desired number of positevely scored metabolites/reactions (used to automatically set FDRs)
#' @param rxn.bias Bias towards the reactions, on the log2 scale ratio of positive scores for reactions and metabolites
#' @return Experiment set object with scored network subnet.scored field
#' @import igraph 
#' @examples
#' library("GAM.networks")
#' library("GAM.db")
#' data(kegg.mouse.network)
#' data(examplesGAM)
#' es.re.scored <- scoreNetwork(es.re, met.fdr=3e-5, rxn.fdr=3e-5, absent.met.score=-20)
#' @export 
scoreNetwork <- function(es,                         
                       fdr=NULL,
                         met.fdr=NA, 
                         rxn.fdr=NA,
                       absent.met.score=NULL,                       
                       absent.rxn.score=NULL,
                       met.score=0,
                       rxn.score=0,
                       num.positive=150,
                       discount.repeats=FALSE,
                       rxn.bias=0) {
    
    net <- es$subnet
    
    if (!is.null(fdr)) {
        met.fdr <- fdr
        rxn.fdr <- fdr
    }
    
    V(net)$score <- 0
    E(net)$score <- 0
    
    ions.to.score <- unique(V(net)[nodeType == "met"]$ion)
    met.origin <- V(net)[nodeType == "met"]$ion
    scoreMets <- function(scores) { V(net)[nodeType == "met"]$score <<- scores[met.origin] }
    if (es$graph.type == "reactions_as_nodes") {
        genes.to.score <- unique(V(net)[nodeType == "rxn"]$origin)
        rxn.origin <- V(net)[nodeType == "rxn"]$origin
        scoreRxns <- function(scores) { V(net)[nodeType == "rxn"]$score <<- scores[rxn.origin] }
    } else {
        genes.to.score <- unique(E(net)$origin)
        rxn.origin <- E(net)$origin
        scoreRxns <- function(scores) { E(net)$score <<- scores[rxn.origin] }
    }
                    
    if (!is.null(es$fb.met) && !is.null(met.fdr)) {            
        fb <- es$fb.met
        pvals <- with(es$met.de.ext, { 
            x <- setNames(pval, ion); 
            x <- na.omit(x)
            x <- x[!duplicated(names(x))]
        })
            
        if (is.na(met.fdr)) {
            met.fdr <- recommendedFDR(fb, pvals, num.positive=num.positive)
            message(sprintf("Using FDR of %.1e for metabolites", met.fdr))
        }
        
        ion.scores <- scoreValue(fb, pvals, met.fdr)
                
        if (is.null(absent.met.score)) {
            absent.met.score <- min(ion.scores)
            message(sprintf("Set score of absent metabolites to %.1f", absent.met.score))            
        }
        
        absent.met.scores <- sapply(ions.to.score, function(x) absent.met.score)
        
        ion.scores <- c(ion.scores, absent.met.scores)
        ion.scores <- ion.scores[names(ion.scores) %in% ions.to.score]
        
        if (discount.repeats) {
            positive <- names(which(ion.scores > 0))
            counts <- table(met.origin)
            counts <- setNames(as.vector(counts), names(counts))
            positiveCounts <- counts[names(counts) %in% positive]
            toMask <- names(positiveCounts[positiveCounts >= 2])
            ion.scores[toMask] <- ion.scores[toMask] / counts[toMask]
        }
        
        
        es$met.fdr <- met.fdr
    } else {
        ion.scores <- rep(met.score, length(ions.to.score))
        names(ion.scores) <- ions.to.score
        es$met.score <- met.score
    }
    
    if (!is.null(es$fb.gene) && !is.null(rxn.fdr)) {
        fb <- es$fb.gene 
        pvals <- es$gene.pval
                
        if (is.na(rxn.fdr)) {
            rxn.fdr <- GAM:::recommendedFDR(fb, pvals, num.positive=num.positive)
            message(sprintf("Using FDR of %.1e for reactions", rxn.fdr))
        }
        
        gene.scores <- GAM:::scoreValue(fb, pvals, rxn.fdr)
        
        if (is.null(absent.rxn.score)) {
            absent.rxn.score <- min(gene.scores)
            if (length(setdiff(genes.to.score, names(gene.scores))) > 0) {            
                message(sprintf("Set score of absent reactions to %.1f", absent.rxn.score))            
            }            
        }
        
        absent.gene.scores <- sapply(genes.to.score, function(x) absent.rxn.score)
        
        gene.scores <- c(gene.scores, absent.gene.scores)
        gene.scores <- gene.scores[names(gene.scores) %in% genes.to.score]
        gene.scores <- gene.scores[genes.to.score]

        if (discount.repeats) {
            positive <- names(which(gene.scores > 0))
            counts <- table(rxn.origin)
            counts <- setNames(as.vector(counts), names(counts))
            positiveCounts <- counts[names(counts) %in% positive]
            toMask <- names(positiveCounts[positiveCounts >= 2])
            gene.scores[toMask] <- gene.scores[toMask] / counts[toMask]
        }
        
        es$rxn.fdr <- rxn.fdr
    } else {
        gene.scores <- rep(rxn.score, length(genes.to.score))
        names(gene.scores) <- genes.to.score
        es$rxn.score <- rxn.score
    }
   
    # :ToDo: check if needed 
    # if (!is.null(rxn.bias)) {
    #     met.pos.sum <- sum(met.scores[met.scores > 0])
    #     rxn.pos.sum <- sum(rxn.scores[rxn.scores > 0])
    #     if (met.pos.sum > 1e-10 && rxn.pos.sum > 1e-10) {
    #         met.scores <- met.scores / met.pos.sum * 100 * (2 ^ (-rxn.bias/2))
    #         rxn.scores <- rxn.scores / rxn.pos.sum * 100 * (2 ^ (rxn.bias/2))
    #     }
    # }
    # 
    scoreMets(ion.scores)
    scoreRxns(gene.scores)
    
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
#' library("GAM.networks")
#' library("GAM.db")
#' data(kegg.mouse.network)
#' data(examplesGAM)
#' es.re.scored <- scoreNetworkWithoutBUM(es.re,
#'                                        met.pval.threshold=1e-5,
#'                                        met.pval.default=1,
#'                                        rxn.pval.threshold=1e-5)
#' @export 
scoreNetworkWithoutBUM <- function(es,
                                   met.pval.threshold=1e-5,
                                   met.pval.default=1,
                                   rxn.pval.threshold=1e-5,
                                   rxn.pval.default=1
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
    
    if (es$graph.type == "reactions_as_nodes") {
        rxn.scores <- rxn.scores[names(rxn.scores) %in% V(net)$name]
        V(net)[names(rxn.scores)]$score <- rxn.scores
    } else {
        E(net)$score <- rxn.scores[E(net)$rxn]
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
#' library("GAM.networks")
#' library("GAM.db")
#' data(kegg.mouse.network)
#' data(examplesGAM)
#' es.re <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M0.M1,
#'                            gene.de=gene.de.M0.M1,
#'                            graph.type="reactions_as_edges")
#' solver <- heinz.solver("/usr/local/lib/heinz/heinz.py")
#' \dontrun{
#' module.re <- findModule(es.re, solver, met.fdr=3e-5, rxn.fdr=3e-5, absent.met.score=-20)
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
