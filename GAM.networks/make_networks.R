#!/usr/bin/env Rscript
library("devtools")
library(methods)
library(BioNet)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(org.Sc.sgd.db)
library(org.At.tair.db)

library(GAM.db)

load_all("../GAM")


data("kegg.db")

#kegg.db$net <- rbind(kegg.db$net, c("C00417", "RXIRG1", "C00490", "RPIRG1", "main"))
#kegg.db$rxn2enz <- rbind(kegg.db$rxn2enz, c("RXIRG1", "Irg1"))
#kegg.db$enz2gene <- rbind(kegg.db$enz2gene, c("Irg1", "16365", "MMU"))
#kegg.db$enz2gene <- rbind(kegg.db$enz2gene, c("Irg1", "730249", "HSA"))

kegg.mouse.network <- makeKeggNetwork(kegg.db, "MMU")
save(kegg.mouse.network, file="data/kegg.mouse.network.rda", compress="xz")

kegg.human.network <- makeKeggNetwork(kegg.db, "HSA")
save(kegg.human.network, file="data/kegg.human.network.rda", compress="xz")

kegg.arabidopsis.network <- makeKeggNetwork(kegg.db, "ATH")
save(kegg.arabidopsis.network, file="data/kegg.arabidopsis.network.rda", compress="xz")

kegg.yeast.network <- makeKeggNetwork(kegg.db, "SCE")
save(kegg.yeast.network, file="data/kegg.yeast.network.rda", compress="xz")

