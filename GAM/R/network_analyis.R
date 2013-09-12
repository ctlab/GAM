plotModuleWithAttr <- function(module, attr="label", layout=layout.kamada.kawai, ...) {
    if (!(attr %in% names(nodeDataDefaults(module)))) {
        return()        
    }
    # Hack for coloring
    
    if ("norm.logFC" %in% nodeDataDefaults(module)) {
        de <- unlist(nodeData(module, nodes(module), "norm.logFC"))
        de <- de * 10
    } else {
        de <- unlist(nodeData(module, nodes(module), "logFC"))    
    }
    de[is.na(de)] <- 0
    
    de[de < 0] <- de[de < 0] - 1
    de[de > 0] <- de[de > 0] + 1
    
    attr.labels = na.omit(unlist(nodeData(module, nodes(module), attr=attr)))
    node.labels = nodes(module)
    names(node.labels) <-node.labels
    all.labels = c(attr.labels, node.labels[!names(node.labels) %in% names(attr.labels)])
    plotModule(module, diff.expr=de, layout=layout, labels=all.labels[node.labels], ...)        
}

save_module <- function(module, outputFilePrefix) {
    outdir <- dirname(outputFilePrefix)
    
    if (!file.exists(outdir)) {
        dir.create(outdir)
    }
    
    
    #module <- subNetwork(nodes(module), network=subnet)
    module    

        
    pdf(paste(outputFilePrefix, "pdf", sep="."), width=15, height=15)
    
    plotModuleWithAttr(module, "label", vertex.size=2)
        
    dev.off()
    
    saveNetwork(module,name=outputFilePrefix,file=outputFilePrefix,type="sif")
    #saveNetwork(module,name=basename(outputFilePrefix),file=outputFilePrefix,type="XGMML")
}

add_norm.diff.expr <- function(module.graph) {
    nodeDataDefaults(module.graph, "norm.logFC") <- NA
    module.nodes <- nodes(module.graph)
    node.types <- unique(na.omit(unlist(nodeData(module.graph, module.nodes, "nodeType"))))
        
    for (node.type in node.types) {
        module.type.nodes <- module.nodes[which(unlist(nodeData(module.graph, module.nodes, "nodeType")) == node.type)]
        module.type.de <- unlist(nodeData(module.graph, module.type.nodes, "logFC"))        
        module.type.norm.de <- module.type.de / max(abs(module.type.de))        
        nodeData(module.graph, module.type.nodes, "norm.logFC") <- unname(module.type.norm.de)
    }
    
    module.graph
}

scoreValue <- function (fb, pval, fdr = 0.01) 
{
    return((fb$a - 1) * (log(pval) - log(fdrThreshold(fdr, fb))))
}


# :ToDo: change name (and function)
.process_pval_and_met_de <- function(es, met.ids, gene.ids) {
    if (!is.null(es$gene.de)) {
        print("Processing gene p-values...")
        es$gene.de$logFC <- fix_inf(es$gene.de$logFC)
        if (!is.null(gene.ids)) {
            data("gene.id.map")
            es$gene.de <- convert.pval(es$gene.de, 
                                       from=gene.id.map[,gene.ids], 
                                       to=gene.id.map[,es$network$gene.ids])
        }    
        es$gene.de$origin <- NULL
        print("Converting gene p-values to reactions...")
        es$rxn.de <- convert.pval(es$gene.de, from=es$network$rxn2gene$gene, to=es$network$rxn2gene$rxn)
        
    }
    
    if (!is.null(es$met.de)) {
        print("Processing metabolite p-values...")
        es$met.de$logFC <- fix_inf(es$met.de$logFC)        
        if (!is.null(met.ids)) {
            data("met.id.map")
            es$met.de <- convert.pval(es$met.de, 
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

graphNEL.from.tables <- function(node.table, edge.table, node.col=1, edge.cols=c(1,2), directed=T, ignore.solitary.nodes=T) {    
    if (is.character(node.col)) {
        node.col <- match(node.col, colnames(node.table))
    }
    if (is.character(edge.cols)) {
        edge.cols <- match(edge.cols, colnames(edge.table))
    }
    
    net1 <- graph.edgelist(as.matrix(edge.table[,edge.cols]), directed=directed)
    net1 <- simplify(net1, remove.multiple=T)
    net1 <- igraph.to.graphNEL(net1)
    
    if (ignore.solitary.nodes) {        
        node.table <- node.table[node.table[,node.col] %in% nodes(net1),]
    }
    
    for (node_attr in colnames(node.table)[-node.col]) {
        if (node_attr == "name") {
            next
        }
        nodeDataDefaults(net1, node_attr) <- NA
        nodeData(net1, n=node.table[,node.col], attr=node_attr) <- node.table[,node_attr]
    }
    
    for (edge_attr in colnames(edge.table)[-edge.cols]) {
        if (edge_attr == "name") {
            next
        }
        edgeDataDefaults(net1, edge_attr) <- NA
        edgeData(net1, from=edge.table[,edge.cols[1]], to=edge.table[,edge.cols[2]], attr=edge_attr) <- 
            edge.table[,edge_attr]        
    }
    return(net1)
}

make_experiment_set <- function(network, 
                                met.de=NULL, gene.de=NULL, rxn.de=NULL,
                                met.ids=NULL, gene.ids=NULL                                
) {
    es <- newEmptyObject()
    es$network <- network
    es$met.de <- met.de
    es$gene.de <- gene.de
    es$rxn.de <- rxn.de 
    
    
    es$graph.raw <- es$network$graph.raw
    
    
    es <- .process_pval_and_met_de(es, met.ids=met.ids, gene.ids=gene.ids)
    
    if (!is.null(es$rxn.de)) {
        print("Processing reaction p-values...")
        es$rxn.de$logFC <- fix_inf(es$rxn.de$logFC)                        
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
        
    print("Processing all p-values together")    
    
    net.edges.ext <- merge(es$graph.raw, es$rxn.de, by.x="rxn", by.y="ID")
    edges2rev <- net.edges.ext$met.x > net.edges.ext$met.y
    net.edges.ext[edges2rev, c("met.x", "met.y")] <- 
        net.edges.ext[edges2rev, c("met.y", "met.x")]
    
    
    
    net.edges.pval <- (aggregate(pval ~ met.x * met.y, data=net.edges.ext, min))
    net.edges.ext <- merge(net.edges.pval, net.edges.ext)
    net.edges.ext <- net.edges.ext[!duplicated(net.edges.ext[,c("met.x", "met.y")]),]
    net.edges.ext <- merge(net.edges.ext, es$network$rxn2name, all.x=T)
    net.edges.ext <- net.edges.ext[net.edges.ext$met.x != net.edges.ext$met.y,]    
    net.edges.ext$logPval <- log(net.edges.ext$pval)
    es$net.edges.ext <- net.edges.ext
    
    met.de.ext <- merge(es$met.de, es$network$met2name, by.x="ID", by.y="met")
    met.de.ext$logPval <- log(met.de.ext$pval)
    
    net1 <- graphNEL.from.tables(node.table=met.de.ext, edge.table=net.edges.ext,
                                 node.col="ID", edge.cols=c("met.x", "met.y"),
                                 directed=F, ignore.solitary.nodes=T)



    nodeDataDefaults(net1, "label") <- NA
    nodeData(net1, attr="label") <- es$network$met2name$name[match(nodes(net1), es$network$met2name$met)]
    edgeDataDefaults(net1, "label") <- NA
    edgeData(net1, from=net.edges.ext$met.x, to=net.edges.ext$met.y, attr="label") <- net.edges.ext$name
    
    es$subnet <- net1
    
    
    es$met.de$origin <- NULL
    es$rxn.de$origin <- NULL
    es$all.de <- rbind(es$met.de, es$rxn.de)
    es$all.pval <- c(es$rxn.pval, es$met.pval)
    
    es$fb.all <- fitBumModel(es$all.pval, plot = TRUE)
    return(es)
}

append_module <- function(res, module.graph) {        
    res[[length(res)+1]] <- module.graph
    res
}


find_modules <- function(es,                         
                         fdr=NULL, met.fdr=NULL, gene.fdr=NULL,
                         absent.met.score=NULL,
                         score.separately=F,
                         heinz.py, heinz.nModules=1,
                         heinz.tolerance=10,
                         heinz.subopt_diff=100) {
    
    if (!is.null(fdr)) {
        met.fdr <- fdr
        gene.fdr <- fdr
    }
    
            
    met.scores <- NULL    
    if (!is.null(es$met.de)) {            
        if (score.separately) {
            met.scores <- scoreFunction(es$fb.met, met.fdr)             
        } else {
            met.scores <- scoreFunction(es$fb.all, met.fdr)[names(es$met.pval)]                
        }
    }
    
    
    rxn.scores <- NULL
    if (!is.null(es$rxn.de)) {                    
        if (score.separately) {             
            rxn.scores <- scoreFunction(es$fb.rxn, gene.fdr)
        } else {                
            rxn.scores <- scoreFunction(es$fb.all, gene.fdr)[names(es$rxn.pval)]
        }
    }
    
    
        
    if (is.null(absent.met.score)) {
        absent.met.score <- mean(met.scores[met.scores < 0])
        print(paste0("absent.met.score <- ", absenet.met.score))
    }
    
    nodeDataDefaults(es$subnet, "score") <- absent.met.score
    met.scores <- met.scores[names(met.scores) %in% nodes(es$subnet)]
    nodeData(es$subnet, names(met.scores), "score") <- met.scores
    
    edgeDataDefaults(es$subnet, "score") <- NA
    edgeData(es$subnet, from=es$net.edges.ext$met.x, to=es$net.edges.ext$met.y, attr="score") <- 
        rxn.scores[es$net.edges.ext$rxn]
    
    
    res <- run_heinz(
        subnet=es$subnet, 
        heinz.py=heinz.py, 
        score.edges=T, 
        score.nodes=T,
        heinz.nModules=heinz.nModules, 
        heinz.tolerance=heinz.tolerance,
        heinz.subopt_diff=heinz.subopt_diff)        
        
    return(res)
}


make_experiment_set.sq <- function(network, 
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
    
    es <- .process_pval_and_met_de(es, met.ids=met.ids, gene.ids=gene.ids)
    
    if (!is.null(es$rxn.de)) {
        print("Processing reaction p-values...")
        es$rxn.de$logFC <- fix_inf(es$rxn.de$logFC)                        
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

run_heinz <- function(subnet,
                      heinz.py, 
                      score.edges=F,
                      score.nodes=T,                      
                      heinz.nModules=1, 
                      heinz.tolerance=10,
                      heinz.subopt_diff=100) {
    tmpdir <- tempdir()        
    edges_file <- paste(tmpdir, "edges.txt", sep="/")
    nodes_file <- paste(tmpdir, "nodes.txt", sep="/")
    
    writeHeinzEdges(subnet, file=edges_file, use.score=score.edges)
    writeHeinzNodes(subnet, file=nodes_file, use.score=score.nodes)
    
    
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
        res <- append_module(res, module.graph)
        
    }
    return(res)
}

find_modules.sq <- function(es,                         
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
        res <- run_heinz(
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
    res <- lapply(res, add_norm.diff.expr)
    
    
    return(res)
}
