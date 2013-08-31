#!/usr/bin/env Rscript

library(optparse)
option_list <- list(
#    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
#                help="Print extra output [default]"),
#    make_option(c("-q", "--quietly"), action="store_false",
#                dest="verbose", help="Print little output"),
    make_option(c("-n", "--network-base"),
                dest="network_base",
                help="File basename (without extension) with network",
                metavar="base"),
    make_option(c("-m", "--mask"),
                dest="mets2mask",
                help="Metabolites to mask from network",
                metavar="sif file"),
    make_option(c("-o", "--output"),
                dest="output_file",
                help="Output file prefix to write network to",
                metavar="file")
)

opt <- parse_args(OptionParser(option_list=option_list))


in_network_file <- paste(opt$network_base, "sif", sep=".")
net <- read.table(in_network_file, header=FALSE, colClasses=c("character", "factor", "character"), comment.char="#", col.names=c("u", "e", "v"))


#sort such that rows and columns are diagonal and create input matrix
 #expM<-expM[order(row.names(expM)),order(row.names(expM))]

rm_comp <- function(s) { gsub("(\\[.\\]|\\(.\\))$", "" , s) }

mask_info<-read.table(opt$mets2mask,header=FALSE,colClasses="character")
node_info<-read.table(paste(opt$network_base, "nodeType.noa", sep="_"),header=FALSE,colClasses="character",skip=1,comment.char="#")

to_keep<-node_info[which(node_info[,3]!='gene'),1]
to_keep<-to_keep[!(rm_comp(to_keep) %in% t(mask_info))]

net <- net[(net$u %in% to_keep) & (net$v %in% to_keep),]
net_rev <- cbind(u=net$v, e=as.character(net$e), v=net$u)
net <- unique(rbind(net, net_rev))

net2 <- net
colnames(net2) <- c("v", "e", "w")
net_sq <- merge(net, net2, by="v")
net_sq <- cbind(u=net_sq$u, e="act", v=net_sq$w)

net_sq_plus_1 <- unique(rbind(net, net_sq))

write.table(net_sq_plus_1, opt$output_file, col.names=F, row.names=F, quote=F, sep="\t")
