#' Package for integrated analyses of genomic and metabolomic data.
#' 
#' This package provides functions for integrated metabolic and
#' transcriptional data analisys based on a reaction network.
#' @name GAM-package
#' @aliases GAM
#' @docType package
#' @author Alexey Sergushichev \email{asergushichev@@wustl.edu}
NULL

#' Table with mapping between metabolite IDs
#' For now just HMDB, KEGG and KEGG name
#' @docType data
#' @name met.id.map
NULL

#' Object with tables from KEGG database
#' @docType data
#' @name kegg.db
NULL

#' Mouse reactions network based on KEGG
#' @docType data
#' @name kegg.mouse.network
NULL

#' Human reactions network based on KEGG
#' @docType data
#' @name kegg.human.network
NULL

#' Gene differential expression data for mouse macrophages between M1 and M2 states
#' @docType data
#' @name gene.de.M1.M2
NULL

#' Metabolite differential expression data for mouse macrophages between M1 and M2 states
#' @docType data
#' @name met.de.M1.M2
NULL

#' Example of a module for M1 vs. M2 comparison whith reactions as edges
#' 
#' met.fdr = 3e-5
#' gene.fdr = 3e-5
#' absent.met.score = -20
#' 
#' @docType data
#' @name module.re
NULL



#' Example of a module for M1 vs. M2 comparison whith reactions as nodes
#' 
#' met.fdr = 2e-7
#' gene.fdr = 2e-7
#' absent.met.score = -20
#' 
#' @docType data
#' @name module.rn
NULL
