#!/usr/bin/env Rscript	

library(BioNet)


save_module <- function(module, outputFilePrefix) {
    
    
    #module <- subNetwork(nodes(module), network=subnet)
    module
    
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
    
    pdf(paste(outputFilePrefix, "pdf", sep="."))
    plotModule(module)
    
    
    
    plotModuleWithAttr(module, "Cobra")
    plotModuleWithAttr(module, "Common_name")
    
    dev.off()
    
    saveNetwork(module,name=outputFilePrefix,file=outputFilePrefix,type="sif")
    saveNetwork(module,name=basename(outputFilePrefix),file=outputFilePrefix,type="XGMML")
}

find_module <- function(data.pval, network, fdr, nModules=1, heinz.py=NULL) {
    gene_names<-as.vector(data.pval$ID)
    
    pval<-as.numeric(data.pval$pval)
    names(pval)<-gene_names
    
    network
    network2 <- rmSelfLoops(network)
    network2
    
    subnet<-subNetwork(names(pval),network2)
    subnet
    
    fb <- fitBumModel(pval, plot = TRUE)
    scores <- scoreNodes(subnet, fb, fdr=fdr)
    nodeDataDefaults(subnet, "score") <- -Inf
    nodeData(subnet, nodes(subnet), "score") <- scores
    
    dm <- data.pval[, "logFC"]
    dm[dm == -Inf] = min(dm[dm != -Inf]) - 1
    dm[dm == Inf] = max(dm[dm != Inf]) + 1
    
    # Hack for coloring
    dm[dm < 0] = dm[dm < 0] - 1
    dm[dm > 0] = dm[dm > 0] + 1
    
    names(dm) <- data.pval$ID
    nodeDataDefaults(subnet, "diff.expr") <- 0
    nodeData(subnet, nodes(subnet), "diff.expr") <- dm[nodes(subnet)]
    
    res <- list()
    
    if (!is.null(opt$heinz.py)) {
        tmpdir <- tempdir()
        tmpdir <- "/tmp"
        edges_file <- paste(tmpdir, "edges.txt", sep="/")
        nodes_file <- paste(tmpdir, "nodes.txt", sep="/")
        writeHeinzEdges(subnet, file=edges_file)
        writeHeinzNodes(subnet, file=nodes_file, use.score=T)
        
        
         system2(opt$heinz.py,
                 c("-n", nodes_file,
                   "-e", edges_file,
                   "-N", "True",
                   "-E", "False",
                   "-s", nModules,
                   "--tolerance", "10",
                   "--subopt_diff", "100",
                       "-v"))
        
        
        for (i in 0:(nModules-1)) {
            module = readHeinzGraph(node.file = paste(nodes_file, i, "hnz", sep="."), 
                                    network = subnet)
            res <- c(res, module)            
        }
    
    } else {
        module <- runFastHeinz(subnet, scores)
        res <- c(res, module)
    }
    return(res)
}



network_analysis.main <-function() {
    library(optparse)
    option_list <- list(
        #    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
        #                help="Print extra output [default]"),
        #    make_option(c("-q", "--quietly"), action="store_false",
        #                dest="verbose", help="Print little output"),
        make_option(c("-i", "--pvals"),
                    dest="pvals_file",
                    help="File to read p-values from",
                    metavar="file"),
        make_option(c("-n", "--network"),
                    dest="network_file",
                    help="File to read p-values from",
                    metavar="file"),
        make_option(c("-o", "--output"),
                    dest="output_file",
                    help="Output file prefix to write network to",
                    metavar="file"),
        make_option(c("--heinz"),
                    dest="heinz.py",
                    help="Path to file heinz.py",
                    metavar="file"),
        make_option(c("--fdr"),
                    dest="fdr",
                    help="False discovery rate to coontrol component size",
                    metavar="number")
    )
    opt <- newEmptyObject()
    opt$network_file <- "./networks//mouse1415/mouse1415.nogene.masked.squared.nocomp.noex.hmdb.esc"
    opt$pvals_file <- "./data_new_MandLPSandIFNg-MandIL4/combined.pval"
    
    opt <- parse_args(OptionParser(option_list=option_list))
    
    
    
    data.pval<-read.table(file=opt$pvals_file,head=TRUE,sep="\t")
    network <- loadNetwork.sif(
        paste(opt$network_file, "sif", sep="."),
        list.files(dirname(opt$network_file), paste(basename(opt$network_file), "_\\w+.NA", sep=""), full.names=T)
    )
    
    
    
    modules <- find_module(data.pval=data.pval, network=network, nModules=2, fdr=as.numeric(opt$fdr), heinz.py=opt$heinz.py)
    
    for (i in 1:length(modules)) {
        module <- modules[[i]]
        save_module(module, paste(opt$output_file, opt$fdr, paste("#", i, sep=""), sep="."))
    }
}