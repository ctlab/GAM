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
    outdir <- dirname(outputFilePrefix)
    
    if (!file.exists(outdir)) {
        dir.create(outdir)
    }
    
    
    #module <- subNetwork(nodes(module), network=subnet)
    module    

        
    pdf(paste(outputFilePrefix, "pdf", sep="."), width=10, height=10)
    
    plotModuleWithAttr(module, "shortName", vertex.size=5)
        
    dev.off()
    
    saveNetwork(module,name=outputFilePrefix,file=outputFilePrefix,type="sif")
    #saveNetwork(module,name=basename(outputFilePrefix),file=outputFilePrefix,type="XGMML")
}

# :ToDo: extract function
find_modules <- function(network, 
                         met.de=NULL, gene.de=NULL, rxn.de=NULL,
                         met.ids=NULL, gene.ids=NULL,
                         fdrs=NULL, met.fdrs=NULL, gene.fdrs=NULL,                         
                         score.separately=F,
                         heinz.py=NULL, heinz.nModules=1,
                         heinz.tolerance=10,
                         heinz.subopt_diff=100) {
    graph <- network$graph    
    if (!is.null(fdrs)) {
        met.fdrs <- fdrs
        gene.fdrs <- fdrs
    }
    
    if (!is.null(gene.de)) {
        if (!is.null(gene.ids)) {
            data("gene.id.map")
            gene.de <- convert.pval(gene.de, 
                                    from=gene.id.map[,gene.ids], 
                                    to=gene.id.map[,network$gene.ids])
        }    
        gene.de$origin <- NULL
        print("Converting gene p-values to reactions...")
        rxn.de <- convert.pval(gene.de, from=network$rxn2gene$gene, to=network$rxn2gene$rxn)
    }
    
    if (!is.null(met.de)) {
        print("Processing metabolite p-values...")
        if (!is.null(met.ids)) {
            data("met.id.map")
            met.de <- convert.pval(met.de, 
                                   from=met.id.map[,met.ids], 
                                   to=met.id.map[,network$met.ids])
        }
        met.de$origin <- NULL
        
        met.pval <- met.de$pval
        names(met.pval) <- met.de$ID
        fb.met <- fitBumModel(met.pval, plot = TRUE)
    } else {
        met.pval <- NULL
    }
    
    
    if (!is.null(rxn.de)) {
        print("Processing reaction p-values...")
        if ("origin" %in% colnames(rxn.de)) {
            graph <- convert.node.names(graph, rxn.de$ID, rxn.de$origin)
            
            rxn.de.orig <- rxn.de
            rxn.de$ID <- rxn.de$origin
            
            rxn.de <- rxn.de[rxn.de$ID %in% nodes(graph), ]
            nodeData(graph, rxn.de$ID, "shortName") <- 
                gene.id.map$name[match(rxn.de$origin, gene.id.map[,network$gene.ids])]
            rxn.de$origin <- NULL
        }
        
        rxn.pval <- rxn.de$pval
        names(rxn.pval) <- rxn.de$ID
        
        fb.rxn <- fitBumModel(rxn.pval, plot = TRUE)
    } else {
        rxn.pval <- NULL
    }
    
    all.de <- rbind(met.de, rxn.de)
    all.pval <- c(rxn.pval, met.pval)

    print("Processing all p-values together")    
    subnet <- subNetwork(all.de$ID, graph)
        
    dm <- all.de[, "logFC"]
    dm[dm == -Inf] <- min(dm[dm != -Inf]) - 1
    dm[dm == Inf] <- max(dm[dm != Inf]) + 1
    
    # Hack for coloring
    dm[dm < 0] <- dm[dm < 0] - 1
    dm[dm > 0] <- dm[dm > 0] + 1
    
    names(dm) <- all.de$ID
    nodeDataDefaults(subnet, "diff.expr") <- 0
    nodeData(subnet, nodes(subnet), "diff.expr") <- dm[nodes(subnet)]
    
    
    fb.all <- fitBumModel(all.pval, plot = TRUE)
    
    res <- list()    
    
    L <- max(length(met.fdrs), length(gene.fdrs))
    for (j in 1:L) {    
        print(paste("Searching for", j, "out of", L, "module groups"))
        
        met.scores <- NULL
        met.fdr <- NULL
        if (!is.null(met.de)) {            
            met.fdr <- met.fdrs[j]    
            print(paste("met.fdr =", met.fdr))
            if (score.separately) {
                met.scores <- scoreFunction(fb.met, met.fdr)             
            } else {
                met.scores <- scoreFunction(fb.all, met.fdr)[names(met.pval)]                
            }
        }
        
        rxn.scores <- NULL
        gene.fdr <- NULL
        if (!is.null(rxn.de)) {            
            gene.fdr <- gene.fdrs[j] 
            print(paste("gene.fdr =", gene.fdr))
            if (score.separately) {             
                rxn.scores <- scoreFunction(fb.rxn, gene.fdr)
            } else {                
                rxn.scores <- scoreFunction(fb.all, gene.fdr)[names(rxn.pval)]
            }
        }
        
        all.scores <- c(met.scores, rxn.scores)        
        scores <- all.scores[nodes(subnet)]
        
        nodeDataDefaults(subnet, "score") <- -Inf
        nodeData(subnet, nodes(subnet), "score") <- scores
        
        append_module <- function(res, graph, n=NULL) {
            module=newEmptyObject()
            module$graph = graph
            module$met.fdr = met.fdr
            module$gene.fdr = gene.fdr
            module$n = n            
            res[[length(res)+1]] <- module
            res
        }    
        
        if (!is.null(heinz.py)) {
            tmpdir <- tempdir()        
            edges_file <- paste(tmpdir, "edges.txt", sep="/")
            nodes_file <- paste(tmpdir, "nodes.txt", sep="/")
            writeHeinzEdges(subnet, file=edges_file)
            writeHeinzNodes(subnet, file=nodes_file, use.score=T)
            
            
            system2(heinz.py,
                    c("-n", nodes_file,
                      "-e", edges_file,
                      "-N", "True",
                      "-E", "False",
                      "-s", heinz.nModules,
                      "--tolerance", heinz.tolerance,
                      "--subopt_diff", heinz.subopt_diff,
                      "-v")); date()
            
            
            for (i in 0:(heinz.nModules-1)) {
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
