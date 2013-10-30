#!/usr/bin/env Rscript
reflink <- read.csv("./misc/reflink.txt", sep="\t", header=T, colClasses="character")

gene.id.map <- reflink[, c("mrnaAcc", "locusLinkId", "name")]
colnames(gene.id.map) <- c("RefSeq", "Entrez", "name")
# Not compressing now for faster loading
save(gene.id.map, file="./GAM//data/gene.id.map.rda", compress=F)

hmdb2kegg <- read.csv("./misc/hmdb2kegg.tsv", header=T, colClasses="character", sep="\t")
kegg2name <- read.table("./networks//kegg/met2name.tsv", header=T, colClasses="character")
kegg2name <- kegg2name[, c("met", "name")]
kegg2name$name <- gsub(";.*$", "", kegg2name$name)
colnames(kegg2name) <- c("KEGG", "name")
met.id.map <- merge(hmdb2kegg, kegg2name, all.x=T)
save(met.id.map, file="./GAM//data/met.id.map.rda", compress=F)
