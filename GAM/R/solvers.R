appendModule <- function(res, module.graph) {        
    res[[length(res)+1]] <- module.graph
    res
}

runHeinz <- function(subnet,
                      heinz.py, 
                      score.edges=F,
                      score.nodes=T,                      
                      nModules=1, 
                      tolerance=10,
                      subopt.diff=100,
                      cplexTimeLimit=1e+75) {
    
    graph.dir <- tempfile("graph")
    dir.create(graph.dir)
    edges.file <- file.path(graph.dir, "edges.txt")
    nodes.file <- file.path(graph.dir, "nodes.txt")
    
    writeHeinzEdges(subnet, file=edges.file, use.score=score.edges)
    
    if (!score.nodes) {
        # Hack to make writeHeinzeNodes working
        V(subnet)$score <- 0
    }
    writeHeinzNodes(subnet, file=nodes.file, use.score=T)
    
    
    wd.bak <- getwd()
    heinz.dir <- dirname(heinz.py)
    setwd(heinz.dir)
    heinz.tmpdir <- tempfile("heinztmp")
    system2(file.path(".", basename(heinz.py)),
            c("-n", nodes.file,
              "-e", edges.file,
              "-N", if (score.nodes) "True" else "False",
              "-E", if (score.edges) "True" else "False",              
              "-s", nModules,
              "--tolerance", tolerance,
              "--subopt_diff", subopt.diff,
              "--heinztmp", heinz.tmpdir,
              "--additional", paste0("'cplexTimeLimit ", cplexTimeLimit, "'"),
              "-v"),
            env="ILOG_LICENSE_FILE=access.ilm")
    setwd(wd.bak)
    
    
    res <- list()
    for (i in 0:(nModules-1)) {
        sol.file <- paste(nodes.file, i, "hnz", sep=".")
        if (!file.exists(sol.file)) {
            warning("Solution file not found")
            return(NULL)
        }
        module.graph <- readHeinzGraph(node.file = sol.file,
                                      network = subnet, format="igraph")
        res <- appendModule(res, module.graph)
        
    }
    return(res)
}

#' Solves MWCS instance using mwcs solver (under development) 
#' @param mwcs Path to mwcs executable
#' @param timeLimit Time limit for execution
#' @return solver function
#' @import igraph
#' @exampls 
#' solver <- mwcs.solver("/usr/local/bin/mwcs")
#' @export
mwcs.solver <- function(mwcs, timeLimit=-1) {
    function(network) {
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        
        graph.dir <- tempfile("graph")
        dir.create(graph.dir)
        edges.file <- file.path(graph.dir, "edges.txt")
        nodes.file <- file.path(graph.dir, "nodes.txt")
        
        writeHeinzEdges(network, file=edges.file, use.score=score.edges)
        
        if (!score.nodes) {
            # Hack to make writeHeinzeNodes working
            V(subnet)$score <- 0
        }
        writeHeinzNodes(network, file=nodes.file, use.score=T)
        
        solution.file <- file.path(graph.dir, "sol.txt")
        
        system2(paste0(mwcs),
                c("-n", nodes.file,
                  "-e", edges.file,
                  "-o", solution.file,
                  "-v", 1,
                  "-t", timeLimit))
        
        
        if (!file.exists(solution.file)) {
            warning("Solution file not found")
            return(NULL)
        }
        res <- readHeinzGraph(node.file = solution.file,
                              network = network, format="igraph")
        return(res)
    }
}

#' Solves MWCS instance using heinz
#' @param heinz.py Path to heinz.py executable
#' @param nModules Number of modules to search for
#' @param tolerance tolerance parameter for heinz
#' @param subopt.diff subopt_diff parameter for heinz
#' @param timeLimit Time limit for execution
#' @return solver function
#' @exampls 
#' solver <- heinz.solver("/usr/local/lib/heinz/heinz.py")
#' @export
heinz.solver <- function(heinz.py,
                         nModules=1,
                         tolerance=10,
                         subopt.diff=100,
                         timeLimit=1e+75
                         ) {
    
    
    function(network) {
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        
        res <- runHeinz(
            subnet=network, 
            heinz.py=heinz.py, 
            score.edges=score.edges,
            score.nodes=score.nodes,
            nModules=nModules, 
            tolerance=tolerance,
            subopt.diff=subopt.diff,
            cplexTimeLimit=timeLimit
            )        
    }
}

#' Solves MWCS instance with BioNet::runFastHeinz algorithm
#' @param network Netowrk to find module in
#' @return Module
#' @exampls 
#' data(kegg.mouse.network)
#' data(examples)
#' es.rn <- makeExperimentSet(network=kegg.mouse.network,
#'                            met.de=met.de.M1.M2,
#'                            gene.de=gene.de.M1.M2,
#'                            reactions.as.edges=F)
#' module.rn <- findModule(es.rn, solver=fastHeinz.solver, met.fdr=1e-3, gene.fdr=1e-3, absent.met.score=-20)
#' @export
fastHeinz.solver <- function(network) {
    score.edges <- "score" %in% list.edge.attributes(network) && !all(E(network)$score == 0)
    if (score.edges) {
        stop("Can't run fast heinz on network with scored edges")
    }
    scores <- V(network)$score
    names(scores) <- V(network)$name
    res <- list(runFastHeinz(network, scores))
}
