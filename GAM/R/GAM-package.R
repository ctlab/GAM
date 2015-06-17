#' Package for integrated analyses of genomic and metabolomic data.
#' 
#' This package provides functions for integrated metabolic and
#' transcriptional data analisys based on a reaction network.
#' @name GAM-package
#' @aliases GAM
#' @docType package
#' @author Alexey Sergushichev \email{asergushichev@@wustl.edu}
NULL

#' Gene differential expression data for mouse macrophages between M1 and M2 states
#' @docType data
#' @name gene.de.M0.M1
NULL

#' Metabolite differential expression data for mouse macrophages between M1 and M2 states
#' @docType data
#' @name met.de.M0.M1
NULL

#' Example of an experiment set object for M1 vs. M2 comparison whith reactions as edges
#'
#' @docType data
#' @name es.re
NULL

#' Example of an experiment set object for M1 vs. M2 comparison whith reactions as nodes
#'
#' @docType data
#' @name es.rn
NULL

#' Example of a module for M1 vs. M2 comparison whith reactions as edges
#' 
#' met.fdr = 3e-5
#' rxb,fdr = 3e-5
#' absent.met.score = -20
#' 
#' @docType data
#' @name module.re
NULL



#' Example of a module for M1 vs. M2 comparison whith reactions as nodes
#' 
#' met.fdr = 2e-7
#' rxb,fdr = 2e-7
#' absent.met.score = -20
#' 
#' @docType data
#' @name module.rn
NULL
