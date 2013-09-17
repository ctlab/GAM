kegg.db <- newEmptyObject()

kegg.db$enz2gene <- read.table("./networks//kegg/enz2gene.tsv", header=T, colClasses="character")
kegg.db$rxn2enz <- read.table("./networks//kegg/rxn2enz.tsv", header=T, colClasses="character")
kegg.db$net <- read.table("./networks//kegg/net.sif", header=T, colClasses="character")
kegg.db$rxn2name <- read.table("./networks//kegg/rxn2name.tsv", header=T, colClasses="character")
kegg.db$met2name <- read.table("./networks//kegg/met2name.tsv", header=T, colClasses="character")

kegg.db$mets2mask <- read.table("./networks//kegg/mets2mask.lst", colClasses="character")[,1]
kegg.db$rxns2mask <- read.table("./networks//kegg/rxns2mask.lst", colClasses="character")[,1]
kegg.db$mets2collapse <- read.table("./networks//kegg/mets2collapse.tsv", header=T, colClasses="character")

save(kegg.db, file="./GAM//data/kegg.db.rda")


kegg.db$net <- rbind(kegg.db$net, c("C00417", "RXIRG1", "C00490"))
kegg.db$rxn2enz <- rbind(kegg.db$rxn2enz, c("RXIRG1", "Irg1"))
kegg.db$enz2gene <- rbind(kegg.db$enz2gene, c("Irg1", "16365", "MMU"))
kegg.db$enz2gene <- rbind(kegg.db$enz2gene, c("Irg1", "730249", "HSA"))

kegg.mouse.network <- makeKeggNetwork(kegg.db, "MMU")
saveNetwork(kegg.mouse.network$graph.sq, file="./networks/kegg/mouse.net.sq", type="sif")
save(kegg.mouse.network, file="./GAM//data/kegg.mouse.network.rda")

kegg.human.network <- makeKeggNetwork(kegg.db, "HSA")
saveNetwork(kegg.human.network$graph.sq, file="./networks/kegg/human.net.sq", type="sif")
save(kegg.human.network, file="./GAM//data/kegg.human.network.rda")

