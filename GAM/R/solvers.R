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
                      subopt_diff=100,
                      cplexTimeLimit=1e+75) {
    
    graph.dir <- tempfile("graph")
    dir.create(graph.dir)
    edges_file <- paste(graph.dir, "edges.txt", sep="/")
    nodes_file <- paste(graph.dir, "nodes.txt", sep="/")
    
    writeHeinzEdges(subnet, file=edges_file, use.score=score.edges)
    
    if (!score.nodes) {
        # Hack to make writeHeinzeNodes working
        V(subnet)$score <- 0
    }
    writeHeinzNodes(subnet, file=nodes_file, use.score=T)
    
    
    wd.bak <- getwd()
    heinz.dir <- dirname(heinz.py)
    setwd(heinz.dir)
    heinz.tmpdir <- tempfile("heinztmp")
    system2(paste0("./", basename(heinz.py)),
            c("-n", nodes_file,
              "-e", edges_file,
              "-N", if (score.nodes) "True" else "False",
              "-E", if (score.edges) "True" else "False",              
              "-s", nModules,
              "--tolerance", tolerance,
              "--subopt_diff", subopt_diff,
              "--heinztmp", heinz.tmpdir,
              "--additional", paste0("'cplexTimeLimit ", cplexTimeLimit, "'"),
              "-v"),
            env="ILOG_LICENSE_FILE=./access.ilm")
    setwd(wd.bak)
    
    
    res <- list()
    for (i in 0:(nModules-1)) {
        sol.file <- paste(nodes_file, i, "hnz", sep=".")
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
#' @export
mwcs.solver <- function(mwcs, timeLimit=-1) {
    function(network) {
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        
        graph.dir <- tempfile("graph")
        dir.create(graph.dir)
        edges.file <- paste(graph.dir, "edges.txt", sep="/")
        nodes.file <- paste(graph.dir, "nodes.txt", sep="/")
        
        writeHeinzEdges(network, file=edges.file, use.score=score.edges)
        
        if (!score.nodes) {
            # Hack to make writeHeinzeNodes working
            V(subnet)$score <- 0
        }
        writeHeinzNodes(network, file=nodes.file, use.score=T)
        
        solution.file <- paste(graph.dir, "sol.txt", sep="/")
        
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
#' @param subopt_diff subopt_diff parameter for heinz
#' @param timeLimit Time limit for execution
#' @return solver function
#' @export
heinz.solver <- function(heinz.py,
                         nModules=1,
                         tolerance=10,
                         subopt_diff=100,
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
            subopt_diff=subopt_diff,
            cplexTimeLimit=timeLimit
            )        
    }
}

#' Solves MWCS instance with BioNet::runFastHeinz algorithm
#' @param network Netowrk to find module in
#' @return Module
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