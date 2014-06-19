#!/usr/bin/env Rscript
library(methods)
library(BioNet)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

library("GAM")

load("kegg.db.rda")


#kegg.db$net <- rbind(kegg.db$net, c("C00417", "RXIRG1", "C00490", "RPIRG1", "main"))
#kegg.db$rxn2enz <- rbind(kegg.db$rxn2enz, c("RXIRG1", "Irg1"))
#kegg.db$enz2gene <- rbind(kegg.db$enz2gene, c("Irg1", "16365", "MMU"))
#kegg.db$enz2gene <- rbind(kegg.db$enz2gene, c("Irg1", "730249", "HSA"))

kegg.mouse.network <- makeKeggNetwork(kegg.db, "MMU")
save(kegg.mouse.network, file="kegg.mouse.network.rda", compress="xz")

kegg.human.network <- makeKeggNetwork(kegg.db, "HSA")
save(kegg.human.network, file="kegg.human.network.rda", compress="xz")


