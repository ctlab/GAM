#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
other.name <- paste(sep="/", script.dirname, "transformations.R")
source(other.name)

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
                metavar="sif file"),
    make_option(c("-o", "--output"),
                dest="output_file",
                help="Output file prefix to write network to",
                metavar="file")
)

opt <- parse_args(OptionParser(option_list=option_list))


filename <- opt$network_file
matrsq <- sif2matrix(filename)
pval_file <- opt$pvals_file

nodes_pvals<-read.table(pval_file,header=TRUE)
to_keep<-nodes_pvals[,1]
matrsq_specific<-matrix_subset(matrsq,to_keep)
length(which(matrsq_specific!=0))

#outputfile <- paste(pval_file,'specific','tab', 'x', sep=".")

outputfile <- opt$output_file
coarsed_network<-matrix2tab(matrsq_specific,outputfile)



