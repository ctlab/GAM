#!/usr/bin/env Rscript

#library(biomaRt)
all_exps <- read.csv("all.exp", row.names=1)
#all_conditions <- read.csv("conditions.csv")

#conditions <- all_conditions$condition[match(colnames(all_exps), all_conditions$sample)]
#all_exps <- rbind(as.character(conditions), all_exps)

ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

refseq2entrez_map <- getBM(
    attributes=c("refseq_mrna", "entrezgene"),
    filters="refseq_mrna",
    values=rownames(all_exps),
    mart=ensembl)

a <- merge(refseq2entrez_map, all_exps, by.x="refseq_mrna")
head(a)


#write.csv(all_exps, file="all.exp", row.names=T)
