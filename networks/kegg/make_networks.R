#!/usr/bin/env Rscript
library(methods)
library(BioNet)
library("GAM")

load("kegg.db.rda")

kegg.db$net <- rbind(kegg.db$net, c("C00417", "RXIRG1", "C00490"))
kegg.db$rxn2enz <- rbind(kegg.db$rxn2enz, c("RXIRG1", "Irg1"))
kegg.db$enz2gene <- rbind(kegg.db$enz2gene, c("Irg1", "16365", "MMU"))
kegg.db$enz2gene <- rbind(kegg.db$enz2gene, c("Irg1", "730249", "HSA"))

kegg.mouse.network <- makeKeggNetwork(kegg.db, "MMU")
saveNetwork(kegg.mouse.network$graph.sq, file="mouse.net.sq", type="sif")
save(kegg.mouse.network, file="kegg.mouse.network.rda")

kegg.human.network <- makeKeggNetwork(kegg.db, "HSA")
saveNetwork(kegg.human.network$graph.sq, file="human.net.sq", type="sif")
save(kegg.human.network, file="kegg.human.network.rda")

