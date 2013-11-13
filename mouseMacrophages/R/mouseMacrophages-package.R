#' Package with transcription and metabolic experimental data for mouse macrophages
#'
#' This package provides RNA-seq based transcript expression data and MS-data
#' for metabolites for mouse macrophages under different treatments.
#' 
#' @name mouseMacrophages-package
#' @aliases mouseMacrophages
#' @docType package
#' @author Abhishek Jha, Ching-Cheng Huang, Alexey Sergushichev, Yulia Ivanova, Karina Chmielewski, Kelly Stewart, Bart Everts, Juliet Ashall, Edward J. Pearce, Edward M. Driggers, Maxim N. Artyomov
#' @examples 
#' library("mouseMacrophages")
#' data("mmpData")
#' mmpGeneSet
#' mmpMetSet
NULL

#' ExpressionSet with gene expressions calculated with DESeq from 3' RNA-seq data
#' @docType data
#' @name mmpGeneSet
NULL

#' Mass spectrometry data for metabolites, log-scaled and normalized.
#' @docType data
#' @name mmpMetSet
NULL
