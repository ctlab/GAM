#!/usr/bin/env Rscript
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
source(paste(sep="/", script.dirname, "io.R"))
source(paste(sep="/", script.dirname, "modify_net.R"))



library(optparse)
option_list <- list(
    make_option(c("-i", "--network"),
                dest="network_file",
                help="Name prefix of a file with network (without .sif)",
                metavar="sif file"),
    make_option(c("--cobra2hmdb"),
                dest="cobra2hmdb",
                help="File with mapping between Cobra and HMDB IDs",
                metavar="tsv file"),
    make_option(c("-m", "--mask"),
                dest="mets2mask",
                help="File with metabolite names to be masked out of network",
                metavar="file"),
    make_option(c("-d", "--metabolite-ids"),
                dest="met.ids",
                help="Dictionary with metabolite IDs",
                metavar="tsv file")
    
)

opt<- newEmptyObject()
opt$network_file <- "networks//mouse1415/mouse1415"
opt$cobra2hmdb <- "./misc//HMDB2cobraID_dictionary.noncomp.txt"
opt$mets2mask <- "./misc/mets2mask.noncomp.hmdb.txt"
opt$met.ids <- "./misc/metabolite_ids.tsv"

opt <- parse_args(OptionParser(option_list=option_list))


cobra2hmdb.map <- read.table(opt$cobra2hmdb, header=T, colClasses="character")
to_mask<- read.table(opt$mets2mask,header=FALSE,colClasses="character")[,1]
metabolite.ids <- read.csv(opt$met.ids, header=T, colClasses="character", sep="\t")

net <- network_from_sif(opt$network_file); suffix = ""
net <- rm.nodes(net, net$meta$v[net$meta$nodeType == "gene"]); suffix=paste(suffix, "nogene", sep=".")
net <- rm.comp.network(net); suffix=paste(suffix, "nocomp", sep=".")
net <- rm.nodes(net, net$meta$v[!is.na(match(sapply(net$meta$v, rm.comp), to_mask))]); suffix=paste(suffix, "masked", sep=".")
net <- add.square.network(net); suffix=paste(suffix, "squared", sep=".")
net <- rm.ex.network(net); suffix=paste(suffix, "noex", sep=".")
net <- convert.ids.network(net, cobra2hmdb.map$Cobra, cobra2hmdb.map$HMDB); suffix=paste(suffix, "hmdb", sep=".")
net <- escape.names(net); suffix=paste(suffix, "esc", sep=".")
net$meta <- merge(net$meta, metabolite.ids, by.x = "v", by.y="HMDB", all.x=T)
net$meta$Cobra <- cobra2hmdb.map$Cobra[match(net$meta$v, cobra2hmdb.map$HMDB)]

network_to_sif(net, paste(opt$network_file, suffix, sep=""))
