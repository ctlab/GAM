plotModuleWithAttr <- function(module, attr, ...) {
    if (!(attr %in% names(nodeDataDefaults(module)))) {
        return()        
    }
    # Hack for coloring
    
    de <- unlist(nodeData(module, nodes(module), "norm.diff.expr"))
    de <- de * 10
    de[de < 0] <- de[de < 0] - 1
    de[de > 0] <- de[de > 0] + 1
    
    attr.labels = na.omit(unlist(nodeData(module, nodes(module), attr=attr)))
    node.labels = nodes(module)
    names(node.labels) <-node.labels
    all.labels = c(attr.labels, node.labels[!names(node.labels) %in% names(attr.labels)])
    plotModule(module, diff.expr=de, layout=layout.kamada.kawai, labels=all.labels[node.labels], ...)        
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

add_norm.diff.expr <- function(module.graph) {
    nodeDataDefaults(module.graph, "norm.diff.expr") <- NA
    module.nodes <- nodes(module.graph)
    node.types <- unique(na.omit(unlist(nodeData(module.graph, module.nodes, "nodeType"))))
        
    for (node.type in node.types) {
        module.type.nodes <- module.nodes[which(unlist(nodeData(module.graph, module.nodes, "nodeType")) == node.type)]
        module.type.de <- unlist(nodeData(module.graph, module.type.nodes, "diff.expr"))        
        module.type.norm.de <- module.type.de / max(abs(module.type.de))        
        nodeData(module.graph, module.type.nodes, "norm.diff.expr") <- unname(module.type.norm.de)
    }
    
    module.graph
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
    
    
    es$graph <- es$network$graph
    
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
        es$met.de$origin <- NULL
        
        es$met.pval <- es$met.de$pval
        names(es$met.pval) <- es$met.de$ID
        es$fb.met <- fitBumModel(es$met.pval, plot = TRUE)
    } else {
        es$met.pval <- NULL
    }
    
    
    if (!is.null(es$rxn.de)) {
        print("Processing reaction p-values...")
        es$rxn.de$logFC <- fix_inf(es$rxn.de$logFC)        
        if ("origin" %in% colnames(es$rxn.de)) {
            print("Collapsing reactions by common most significant enzymes")
            #rxn.de.origin.split <- es$rxn.de$origin
            rxn.de.origin.split <- split.mapping.by.connectivity(es$graph, es$rxn.de$ID, es$rxn.de$origin)
            
            es$graph <- convert.node.names(es$graph, es$rxn.de$ID, rxn.de.origin.split)
            
            
            es$rxn.de.orig <- es$rxn.de
            es$rxn.de$ID <- rxn.de.origin.split
            
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
    
    es$all.de <- rbind(es$met.de, es$rxn.de)
    es$all.pval <- c(es$rxn.pval, es$met.pval)
    
    print("Processing all p-values together")    
    es$subnet <- subNetwork(es$all.de$ID, es$graph)
    
    # :ToDo: change Inf value for genes and metabolites
    dm <- es$all.de[, "logFC"]   
    
    
    names(dm) <- es$all.de$ID
    nodeDataDefaults(es$subnet, "diff.expr") <- 0
    nodeData(es$subnet, nodes(es$subnet), "diff.expr") <- dm[nodes(es$subnet)]
    
    
    es$fb.all <- fitBumModel(es$all.pval, plot = TRUE)
    return(es)
}

# :ToDo: extract function
find_modules <- function(es,                         
                         fdrs=NULL, met.fdrs=NULL, gene.fdrs=NULL,
                         score.separately=F,
                         heinz.py=NULL, heinz.nModules=1,
                         heinz.tolerance=10,
                         heinz.subopt_diff=100) {
    
    res <- list()    
    if (!is.null(fdrs)) {
        met.fdrs <- fdrs
        gene.fdrs <- fdrs
    }
    
    L <- max(length(met.fdrs), length(gene.fdrs))
    for (j in 1:L) {    
        print(paste("Searching for", j, "out of", L, "module groups"))
        
        met.scores <- NULL
        met.fdr <- NULL
        if (!is.null(met.de)) {            
            met.fdr <- met.fdrs[j]    
            print(paste("met.fdr =", met.fdr))
            if (score.separately) {
                met.scores <- scoreFunction(es$fb.met, met.fdr)             
            } else {
                met.scores <- scoreFunction(es$fb.all, met.fdr)[names(es$met.pval)]                
            }
        }
        
        rxn.scores <- NULL
        gene.fdr <- NULL
        if (!is.null(rxn.de)) {            
            gene.fdr <- gene.fdrs[j] 
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
        
        append_module <- function(res, module.graph, n=NULL) {
            
            module=newEmptyObject()
            module$graph = add_norm.diff.expr(module.graph)
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
            writeHeinzEdges(es$subnet, file=edges_file)
            writeHeinzNodes(es$subnet, file=nodes_file, use.score=T)
            
            
            system2(heinz.py,
                    c("-n", nodes_file,
                      "-e", edges_file,
                      "-N", "True",
                      "-E", "False",
                      "-s", heinz.nModules,
                      "--tolerance", heinz.tolerance,
                      "--subopt_diff", heinz.subopt_diff,
                      "-v"));
            
            
            for (i in 0:(heinz.nModules-1)) {
                module.graph = readHeinzGraph(node.file = paste(nodes_file, i, "hnz", sep="."), 
                                       network = es$subnet)
                res <- append_module(res, module.graph, i+1)
                
            }
            
        } else {
            module.graph <- runFastHeinz(es$subnet, scores)
            res <- append_module(res, module.graph)
        }
    }
    return(res)
}
