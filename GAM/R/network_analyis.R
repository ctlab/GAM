plotModuleWithAttr <- function(module, attr, ...) {
    if (!(attr %in% names(nodeDataDefaults(module)))) {
        return()        
    }
    attr.labels = na.omit(unlist(nodeData(module, nodes(module), attr=attr)))
    node.labels = nodes(module)
    names(node.labels) <-node.labels
    all.labels = c(attr.labels, node.labels[!names(node.labels) %in% names(attr.labels)])
    plotModule(module, layout=layout.kamada.kawai, labels=all.labels[node.labels], ...)        
}

save_module <- function(module, outputFilePrefix) {
    
    
    #module <- subNetwork(nodes(module), network=subnet)
    module
    

        
    pdf(paste(outputFilePrefix, "pdf", sep="."), width=10, height=10)
    
    plotModuleWithAttr(module, "shortName", vertex.size=5)
        
    dev.off()
    
    saveNetwork(module,name=outputFilePrefix,file=outputFilePrefix,type="sif")
    #saveNetwork(module,name=basename(outputFilePrefix),file=outputFilePrefix,type="XGMML")
}

find_modules <- function(met.de, rxn.de, network, met.fdrs=NULL, rxn.fdrs=NULL, nModules=1, heinz.py=NULL) {
    met.de$origin <- NULL
    rxn.de$origin <- NULL
    data.pval <- rbind(met.de, rxn.de)
    
    rxn.pval <- rxn.de$pval
    names(rxn.pval) <- rxn.de$ID
    met.pval <- met.de$pval
    names(met.pval) <- met.de$ID

        
    network
    network2 <- rmSelfLoops(network)
    network2
    
    subnet<-subNetwork(data.pval$ID,network2)
    subnet
    
    dm <- data.pval[, "logFC"]
    dm[dm == -Inf] = min(dm[dm != -Inf]) - 1
    dm[dm == Inf] = max(dm[dm != Inf]) + 1
    
    # Hack for coloring
    dm[dm < 0] = dm[dm < 0] - 1
    dm[dm > 0] = dm[dm > 0] + 1
    
    names(dm) <- data.pval$ID
    nodeDataDefaults(subnet, "diff.expr") <- 0
    nodeData(subnet, nodes(subnet), "diff.expr") <- dm[nodes(subnet)]
    fb.met <- fitBumModel(met.pval, plot = TRUE)
    fb.rxn <- fitBumModel(rxn.pval, plot = TRUE)
    res <- list()
        
    
    for (j in 1:length(met.fdrs)) {    
        met.fdr <- met.fdrs[j]
        rxn.fdr <- rxn.fdrs[j] 
        met.scores <- scoreFunction(fb.met, met.fdr)
        rxn.scores <- scoreFunction(fb.rxn, rxn.fdr)
        all.scores <- c(met.scores, rxn.scores)
        
        scores <- all.scores[nodes(subnet)]
        nodeDataDefaults(subnet, "score") <- -Inf
        nodeData(subnet, nodes(subnet), "score") <- scores
        
        append_module <- function(res, graph, n=NULL) {
            module=newEmptyObject()
            module$graph = graph
            module$met.fdr = met.fdr
            module$rxn.fdr = rxn.fdr
            module$n = n            
            res[[length(res)+1]] <- module
            res
        }    
        
        if (!is.null(heinz.py)) {
            tmpdir <- tempdir()
            tmpdir <- "/tmp"
            edges_file <- paste(tmpdir, "edges.txt", sep="/")
            nodes_file <- paste(tmpdir, "nodes.txt", sep="/")
            writeHeinzEdges(subnet, file=edges_file)
            writeHeinzNodes(subnet, file=nodes_file, use.score=T)
            
            
            system2(heinz.py,
                    c("-n", nodes_file,
                      "-e", edges_file,
                      "-N", "True",
                      "-E", "False",
                      "-s", nModules,
                      "--tolerance", "10",
                      "--subopt_diff", "100",
                      "-v")); date()
            
            
            for (i in 0:(nModules-1)) {
                graph = readHeinzGraph(node.file = paste(nodes_file, i, "hnz", sep="."), 
                                       network = subnet)
                res <- append_module(res, graph, i+1)
                
            }
            
        } else {
            graph <- runFastHeinz(subnet, scores)
            res <- append_module(res, graph)
        }
    }
    return(res)
}
