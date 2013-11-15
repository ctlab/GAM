#!/usr/bin/env Rscript
#reflink <- read.csv("reflink.txt", sep="\t", header=T, colClasses="character")

#gene.id.map <- reflink[, c("mrnaAcc", "locusLinkId", "name")]
#colnames(gene.id.map) <- c("RefSeq", "Entrez", "name")
# Not compressing now for faster loading
#save(gene.id.map, file="./GAM//data/gene.id.map.rda", compress=F)

hmdb2kegg <- read.csv("metabocards-2.5.extr.tsv", header=T, colClasses="character", sep="\t")
hmdb2kegg <- na.omit(hmdb2kegg[, c("HMDB", "KEGG")])
nrow(hmdb2kegg)
head(hmdb2kegg)

kegg2name <- read.table("../networks/kegg/met2name.tsv", header=T, colClasses="character")
kegg2name <- kegg2name[, c("met", "name")]
kegg2name$name <- gsub(";.*$", "", kegg2name$name)
colnames(kegg2name) <- c("KEGG", "name")
nrow(kegg2name)
head(kegg2name)

met.id.map <- merge(hmdb2kegg, kegg2name, all.x=T)
head(met.id.map)
save(met.id.map, file="met.id.map.rda", compress=F)
