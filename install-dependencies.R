#!/usr/bin/env Rscript

source("http://bioconductor.org/biocLite.R")
biocLite()

installCRAN <- function(list.of.packages) {
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

    if (length(new.packages) > 0) {
        install.packages(new.packages, repos="http://cran.us.r-project.org")
    }
}

installBioc <- function(list.of.packages) {
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

    if (length(new.packages) > 0) {
        biocLite(new.packages)
    }
}

installBioc("RCurl")

installCRAN(c("devtools", "roxygen2", "igraph", "data.table", 
              "plyr", "knitr", "testthat", "rjson"))
              

installBioc(c("BioNet", "org.Mm.eg.db", "org.Hs.eg.db", "biomaRt"))

