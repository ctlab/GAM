appendModule <- function(res, module.graph) {        
    res[[length(res)+1]] <- module.graph
    res
}

runHeinz <- function(subnet,
                      heinz.py, 
                      score.edges=FALSE,
                      score.nodes=TRUE,                      
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
    writeHeinzNodes(subnet, file=nodes.file, use.score=TRUE)
    
    
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
#' @examples 
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
        writeHeinzNodes(network, file=nodes.file, use.score=TRUE)
        
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
#' @examples 
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
#' @examples 
#' data(kegg.mouse.network)
#' data(examplesGAM)
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


#' Solves MWCS instance with randomized heuristic algorithm
#' @param nruns Number of algorithm runs
#' @return solver function
#' @examples 
#' data(kegg.mouse.network)
#' data(examplesGAM)
#' solver <- randHeur.solver()
#' module.re <- findModule(es.re, solver=solver, met.fdr=1e-3, gene.fdr=1e-3, absent.met.score=-20)
#' @export
randHeur.solver <- function(nruns=10) {
    res <- function(network) {        
        jobs <- list()
        for (i in seq_len(nruns)) {
            jobs <- c(jobs, list(mcparallel(solveMwcsRandHeur(network))))
        }
        solutions <- mccollect(jobs)
        best <- which.max(sapply(solutions, function(s) s$score))        
        solution <- induced.subgraph(network, solutions[[best]]$nodes)        
        message("best score: ", solutions[[best]]$score)
        solution
    }
    res
}

solveMwcsRandHeur <- function(network) {
    n <- length(V(network))
    m <- length(E(network))
    vscores <- V(network)$score
    escores <- E(network)$score
    
    score <- function(nodes, edges) {
        sum(vscores[nodes]) + sum(escores[edges])
    }
    
    best.solutions <- list()
    best.solutions.edges <- list()
    best.scores <- rep(0, n)
    edges <- get.edgelist(network, names=FALSE)
    
    solution <- list()
    solution$score <- -Inf    
    
    makeSolution <- function(id, nodes, edges, score=NULL) {
        res <- newEmptyObject()
        res$id <- id
        res$nodes <- nodes
        res$edges <- edges
        res$score <- score       
        
        
        if (is.null(res$score)) {
            res$score <- score(nodes, edges)
        }
        res
        
    }
    
    updateBest <- function(newsolution) {
        if (newsolution$score > solution$score) {
            solution <<- newsolution
            #             message("new best score: ", solution$score)
        }
    }
    
    for (v in V(network)) {
        newsolution <- makeSolution(v, v, NULL)
        updateBest(newsolution);
        best.solutions[[v]] <- newsolution
    }
    
    
    checked <- rep(FALSE, m)
    
    max.iterations <- 10000
    
    for (i in seq_len(max.iterations)) {                
        edges.unchecked <- which(!checked)        
        if (length(edges.unchecked) == 0) {
            #             message("exhausted")            
            break
        }
        edge <- max(ceiling(runif(1, max=length(edges.unchecked))), 1)
        checked[edge] <- TRUE
        
        start <- edges[edge, 1]
        end <- edges[edge, 2]
        
        startsolution <- best.solutions[[start]]
        endsolution <- best.solutions[[end]]
        
        score.start <- startsolution$score
        score.end <- endsolution$score
        
        if (startsolution$id == endsolution$id && edge %in% startsolution$edges) {
            next
        }
        
        newnodes <- union(startsolution$nodes, endsolution$nodes)
        newedges <- union(edge, union(startsolution$edges, endsolution$edges))
        
        
        newsolution <- makeSolution(n + i, newnodes, newedges)
        
        
        to_update <- c()
        
        
        
        if (newsolution$score > score.start) {
            start.nodes <- startsolution$nodes
            to_update <- union(to_update, start.nodes[best.scores[start.nodes] < newsolution$score])
        }
        
        if (newsolution$score > score.end) {
            end.nodes <- endsolution$nodes
            to_update <- union(to_update, end.nodes[best.scores[end.nodes] < newsolution$score])
        }
        
        if (length(to_update) > 0) {
            updateBest(newsolution)
            best.solutions[to_update] <- list(newsolution)             
            to_check <- sapply(E(network)[adj(to_update)], identity)
            checked[to_check] <- FALSE
        }        
    }    
    solution
}
