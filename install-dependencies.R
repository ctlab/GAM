#!/usr/bin/env Rscript

list.of.packages <- c("devtools", "roxygen2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if (length(new.packages) > 0) {
    install.packages(new.packages)
}


list.of.packages <- c("DESeq")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if (length(new.packages) > 0) {
    biocLite(new.packages)
}


