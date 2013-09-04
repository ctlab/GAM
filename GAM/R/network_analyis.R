library(BioNet)

plotModuleWithAttr <- function(module, attr) {
    if (!(attr %in% names(nodeDataDefaults(module)))) {
        return()        
    }
    attr.labels = na.omit(unlist(nodeData(module, nodes(module), attr=attr)))
    node.labels = nodes(module)
    names(node.labels) <-node.labels
    all.labels = c(attr.labels, node.labels[!names(node.labels) %in% names(attr.labels)])
    plotModule(module, labels=all.labels[node.labels])        
}

save_module <- function(module, outputFilePrefix) {
    
    
    #module <- subNetwork(nodes(module), network=subnet)
    module
    

    
    pdf(paste(outputFilePrefix, "pdf", sep="."))
    plotModule(module)
    
    
    
    plotModuleWithAttr(module, "Cobra")
    plotModuleWithAttr(module, "Common_name")
    
    dev.off()
    
    saveNetwork(module,name=outputFilePrefix,file=outputFilePrefix,type="sif")
    #saveNetwork(module,name=basename(outputFilePrefix),file=outputFilePrefix,type="XGMML")
}

find_modules <- function(data.pval, network, fdrs, nModules=1, heinz.py=NULL) {
    gene_names<-as.vector(data.pval$ID)
    
    pval<-as.numeric(data.pval$pval)
    names(pval)<-gene_names
    
    network
    network2 <- rmSelfLoops(network)
    network2
    
    subnet<-subNetwork(names(pval),network2)
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
    fb <- fitBumModel(pval, plot = TRUE)
    res <- list()
    j=1
    
    for (fdr in fdrs) {        
        
        scores <- scoreNodes(subnet, fb, fdr=fdr)
        nodeDataDefaults(subnet, "score") <- -Inf
        nodeData(subnet, nodes(subnet), "score") <- scores
        
        
        
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
                      "-v"))
            
            
            for (i in 0:(nModules-1)) {
                module=newEmptyObject()
                module$graph = readHeinzGraph(node.file = paste(nodes_file, i, "hnz", sep="."), 
                                        network = subnet)
                module$fdr = fdr
                module$n = i+1
                res[[j]] <- module
                j = j + 1
            }
            
        } else {
            module=newEmptyObject()
            module$graph <- runFastHeinz(subnet, scores)
            module$fdr = fdr            
            res[[j]] <- module
            j = j + 1
        }
    }
    return(res)
}
