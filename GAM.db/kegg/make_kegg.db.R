#!/usr/bin/env Rscript
library(methods)
library(plyr)

kegg.db <- newEmptyObject()

kegg.db$enz2gene <- read.table("enz2gene.all.tsv", header=T, colClasses="character")
kegg.db$rxn2enz <- read.table("rxn2enz.tsv", header=T, colClasses="character")

net <- read.table("net.sif", header=T, colClasses="character")
rpairs <- read.table("rpairs.tsv", header=T, colClasses="character")
rpairs <- rbind(rpairs, rename(rpairs, c(met.x="met.y", met.y="met.x")))
kegg.db$net <- merge(net, rpairs, all=T)

kegg.db$rxn2name <- read.table("rxn2name.tsv", header=T, colClasses="character")
kegg.db$met2name <- read.table("met2name.tsv", header=T, colClasses="character")

kegg.db$mets2mask <- read.table("mets2mask.lst", colClasses="character")[,1]
kegg.db$rxns2mask <- read.table("rxns2mask.lst", colClasses="character")[,1]
kegg.db$mets2collapse <- read.table("mets2collapse.tsv", header=T, colClasses="character")

save(kegg.db, file="kegg.db.rda", compress="xz")

