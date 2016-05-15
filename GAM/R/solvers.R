appendModule <- function(res, module.graph) {        
    res[[length(res)+1]] <- module.graph
    res
}

readGraph <- function(node.file, edge.file, network) {
    nodes <- as.matrix(read.table(file = node.file, 
                                  na.strings = "n/a"))
    nodes2 <- which(!is.na(as.numeric(nodes[, 2])))
    
    
    edges <- as.matrix(read.table(file = edge.file,
                                  na.strings = "n/a"))
    edges2 <- which(!is.na(as.numeric(edges[, 3])))
    
    res <- subgraph.edges(network, eids = edges2, delete.vertices = T)
    stopifnot(setequal(V(network)[nodes2]$name, V(res)$name))
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
        module.graph <- readGraph(node.file = sol.file,
                                  edge.file = paste(edges.file, i, "hnz", sep="."),
                                  network = subnet)
        res <- appendModule(res, module.graph)
        
    }
    return(res)
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


MWCSize <- function(g) {         
    was.connected <- is.connected(g)
    et <- as.data.table(get.edge.attributes(g, include.ends=TRUE))
    et <- et[, list(from, to, score)]
    et[, from.weight := pmax(-V(g)[et$from]$score, 0)+1e-3]
    et[, to.weight := pmax(-V(g)[et$to]$score, 0)+1e-3]
    et <- et[score > 0, ]
    
    if (nrow(et) > 0) {
        et[, from.d := score * from.weight / (from.weight + to.weight)]
        et[, to.d := score * to.weight / (from.weight + to.weight)]
        
        ds <- rbind(
            et[, list(v=from, d=from.d)],
            et[, list(v=to, d=to.d)])
        
        ds <- aggregate(d ~ v, data=ds, sum)
        
        V(g)[ds$v]$score <- V(g)[ds$v]$score + ds$d
        E(g)$origEdge <- E(g)
        E(g)[score > 0]$score <- 0        
    }
    
    
    neg.edges <- E(g)[score < 0]
    if (length(neg.edges) > 0) {
        neg.edges.table <- get.edgelist(g)[neg.edges, ]
        neg.edges.names <- sprintf("%s_%s_%s", neg.edges.table[,1], neg.edges.table[,2], neg.edges)
        V(g)$wasEdge <- FALSE
        g <- add.vertices(g, length(neg.edges), name=neg.edges.names, score=E(g)[neg.edges]$score, wasEdge=TRUE, origEdge=neg.edges)
        new.edges <- rbind(
            cbind(neg.edges.table[,1], neg.edges.names),
            cbind(neg.edges.table[,2], neg.edges.names)
        )
        g <- add.edges(g, t(as.matrix(new.edges)), score=0)
        g <- delete.edges(g, E(g)[score < 0])
    }
    
    E(g)$score <- 0
    stopifnot(is.connected(g) == was.connected)
    g
}

deMWCSize <- function(m, g) {
    E(g)$origEdge <- E(g)
    orig.vertices <- V(g)[name %in% V(m)$name]
    edges.to.add <- unique(c(
        E(induced.subgraph(g, orig.vertices))[score >= 0]$origEdge,
        na.omit(c(V(m)$origEdge, E(m)$origEdge))))        
    
    m.x <- subgraph.edges(g, edges.to.add)    
    m.x <- remove.edge.attribute(m.x, "origEdge")
    m.x
}

#' Solves MWCS instance using heinz2 solver
#' @param heinz2 Path to heinz2 executable
#' @param nthreads Number of threads to use
#' @param timeLimit Time limit for execution
#' @return solver function
#' @import igraph
#' @examples 
#' solver <- heinz2.solver("/usr/local/lib/heinz2/heinz")
#' @export
heinz2.solver <- function(heinz2, nthreads=1, timeLimit=-1) {
    function(network) {
        network.orig <- network
        
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        
        
        graph.dir <- tempfile("graph")
        dir.create(graph.dir)
        edges.file <- file.path(graph.dir, "edges.txt")
        nodes.file <- file.path(graph.dir, "nodes.txt")
        
                
        if (!score.nodes) {
            # Hack to make writeHeinzeNodes working
            V(network)$score <- 0
        }
        
        if (score.edges) {
            network <- MWCSize(network.orig)            
        }        
                
        writeHeinzNodes(network, file=nodes.file, use.score=TRUE)
        writeHeinzEdges(network, file=edges.file, use.score=score.edges)
        
        solution.file <- file.path(graph.dir, "sol.txt")
        
        system2(paste0(heinz2),
                c("-n", nodes.file,
                  "-e", edges.file,
                  "-o", solution.file,
                  "-m", nthreads,
                  "-p",
                  "-v", 1,
                  "-t", timeLimit))
        
        
        if (!file.exists(solution.file)) {
            warning("Solution file not found")
            return(NULL)
        }
        res <- readHeinzGraph(node.file = solution.file,
                              network = network, format="igraph")
        
        if (score.edges) {
            res <- deMWCSize(res, network.orig)
        }
        
        return(res)
    }
}


#' Solves MWCS instance with randomized heuristic algorithm
#' @param nruns Number of algorithm runs
#' @param ... Additional arguments for solveMwcsRandHeur
#' @return solver function
#' @examples 
#' library("GAM.networks")
#' library("GAM.db")
#' data(kegg.mouse.network)
#' data(examplesGAM)
#' solver <- randHeur.solver()
#' module.re <- findModule(es.re, solver=solver, met.fdr=1e-3, rxn.fdr=1e-3, absent.met.score=-20)
#' @import parallel
#' @export
randHeur.solver <- function(nruns=10, ...) {
    res <- function(network) {        
        jobs <- list()
        RNGkind("L'Ecuyer-CMRG")
        for (i in seq_len(nruns)) {
            mc.reset.stream()
            jobs <- c(jobs, list(mcparallel({                    
                solveMwcsRandHeur(network, ...)})))
        }
        solutions <- mccollect(jobs)
        best <- which.max(sapply(solutions, function(s) s$score))        
        if (is.null(solutions[[best]]$edges)) {
            solution <- induced.subgraph(network, solutions[[best]]$nodes)
        } else {
            solution <- subgraph.edges(network, solutions[[best]]$edges)        
        }
        message("best score: ", solutions[[best]]$score)
        solution
    }
    res
}

solveMwcsRandHeur <- function(network, max.iterations=10000) {
    n <- length(V(network))
    m <- length(E(network))
    vscores <- V(network)$score
    
    escores <- E(network)$score
    
    score <- function(nodes, edges) {
        sum(vscores[nodes]) + sum(escores[edges])
    }
    
    best.solutions <- list()
    best.solutions.edges <- list()
    
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
        }
    }
    
    for (v in V(network)) {
        newsolution <- makeSolution(v, v, NULL)
        updateBest(newsolution);
        best.solutions[[v]] <- newsolution
    }
    best.scores <- vscores
    
    checked <- rep(FALSE, m)
    
    
    
    for (i in seq_len(max.iterations)) {                
        edges.weights <- pmax(best.scores[edges[, 1]],
                              best.scores[edges[, 2]]) + escores
        edges.unchecked <- which(!checked)         
        
        if (length(edges.unchecked) == 0) {
            break
        }
        
        candidates <- sapply(1:2, function(x) {
            max(ceiling(runif(1, max=length(edges.unchecked))), 1)
        })
        
        edge <- candidates[which.max(edges.weights[candidates])]
        
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
            best.scores[to_update] <- newsolution$score
            to_check <- sapply(E(network)[adj(to_update)], identity)
            checked[to_check] <- FALSE
        }        
    }    
    solution    
   
}

#' Solves GMWCS instance using gmwcs solver
#' @param gmwcs Path to gmwcs executable
#' @param nthreads Number of threads to use
#' @param timeLimit Time limit for execution
#' @return solver function
#' @import igraph
#' @examples 
#' solver <- gmwcs("/usr/local/bin/gmwcs")
#' @export
gmwcs.solver <- function (gmwcs, nthreads = 1, timeLimit = -1) {
  function(network) {
    network.orig <- network
    score.edges <- "score" %in% list.edge.attributes(network)
    score.nodes <- "score" %in% list.vertex.attributes(network)
    graph.dir <- tempfile("graph")
    dir.create(graph.dir)
    edges.file <- file.path(graph.dir, "edges.txt")
    nodes.file <- file.path(graph.dir, "nodes.txt")
    if (!score.nodes) {
      V(network)$score <- 0
    }
    
    BioNet::writeHeinzNodes(network, file = nodes.file, use.score = TRUE)
    BioNet::writeHeinzEdges(network, file = edges.file, use.score = score.edges)
    system2(gmwcs, c("-n", nodes.file, "-e", edges.file, 
                     "-m", nthreads, "-t", timeLimit
                      ,             "-b"
    ))
    solution.file <- paste0(nodes.file, ".out")
    if (!file.exists(solution.file)) {
      warning("Solution file not found")
      return(NULL)
    }
    res <- GAM:::readGraph(node.file = solution.file,
                     edge.file = paste0(edges.file, ".out"),
                     network = network)
    return(res)
  }
}

