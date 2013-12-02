library(igraph)
library(BioNet)

context("Post-processing")

test_that("addNormLogFC works", {
    met.nt <- data.frame(ID=c("C01", "C02"), log2FC=c(1, Inf), stringsAsFactors=F)
    rxn.nt <- data.frame(ID=c("R01", "R02"), log2FC=c(3, 12), stringsAsFactors=F)
    
    et <- data.frame(
        u=c("C01", "C01", "C02", "C02"), 
        v=c("R01", "R02", "R01", "R02"), 
        stringsAsFactors=F)
    
    g <- GAM:::graph.from.tables(node.table=list(met=met.nt, rxn=rxn.nt), edge.table=et, directed=F)
    
    g1 <- addNormLogFC(g)
    
    expect_equal(V(g1)["C01"]$log2FC.norm, 0.5)
    expect_equal(V(g1)["C02"]$log2FC.norm, 1)
    expect_equal(V(g1)["R01"]$log2FC.norm, 0.25)
    expect_equal(V(g1)["R02"]$log2FC.norm, 1)
})

test_that("addInterconnections works", {
    met.nt <- data.frame(ID=c("C01", "C02"), log2FC=c(1, 2), stringsAsFactors=F)
    rxn.nt <- data.frame(ID=c("R01", "R02"), log2FC=c(3, 12), stringsAsFactors=F)
    
    et <- data.frame(
        u=c("C01", "C01", "C02", "C02"), 
        v=c("R01", "R02", "R01", "R02"), 
        stringsAsFactors=F)
    
    g <- GAM:::graph.from.tables(node.table=list(met=met.nt, rxn=rxn.nt), edge.table=et, directed=F)
    es <- newEmptyObject()
    es$subnet <- g
    es$graph.raw <- data.frame(met.x=c("C01", "C01"), rxn=c("R01", "R02"), met.y=c("C02", "C02"), stringsAsFactors=F)
    es$rxn.de.origin.split <- rxn.nt$ID
    names(es$rxn.de.origin.split) <- rxn.nt$ID
    module <- subNetwork(c("C01", "C02", "R01"), es$subnet)
    
    module1 <- addInterconnections(module, es)
    
    expect_true("R02" %in% V(module1)$name)
})

test_that("expandReactionNodeAttributesToEdges works", {
    met.nt <- data.frame(ID=c("C01", "C02"), log2FC=c(1, 2), stringsAsFactors=F)
    rxn.nt <- data.frame(ID=c("R01", "R02"), log2FC=c(3, 12), stringsAsFactors=F)
    
    et <- data.frame(
        u=c("C01", "C01", "C02", "C02"), 
        v=c("R01", "R02", "R01", "R02"), 
        stringsAsFactors=F)
    
    g <- GAM:::graph.from.tables(node.table=list(met=met.nt, rxn=rxn.nt), edge.table=et, directed=F)
    
    g1 <- expandReactionNodeAttributesToEdges(g)
    
    expect_true("log2FC" %in% list.edge.attributes(g1))
    expect_equal(E(g1, P=c("C01", "R01"))$log2FC, 3)
})

test_that("removeHangingNodes works", {
    met.nt <- data.frame(ID=c("C01", "C02"), pval=c(1e-5, NA), stringsAsFactors=F)
    rxn.nt <- data.frame(ID=c("R01", "R02"), pval=c(1e-12, 1e-42), stringsAsFactors=F) 
    
    et <- data.frame(
        u=c("C01", "C01", "C02"), 
        v=c("R01", "R02", "R01"), 
        stringsAsFactors=F)
    
    g <- GAM:::graph.from.tables(node.table=list(met=met.nt, rxn=rxn.nt), edge.table=et, directed=F)
    
    g1 <- removeHangingNodes(g)
    
    expect_true(!"C02" %in% V(g1)$name)
})

test_that("addTransEdges works", {
    met.nt <- data.frame(ID=c("C01", "C02", "C03"), pval=c(1e-5, NA, 1e-6), stringsAsFactors=F)
    et <- data.frame(
        met.x=c("C01", "C02", "C01"),
        met.y=c("C02", "C03", "C03"),
        rxn=c("R01", "R02", "R03"),
        rptype=c("main", "main", "trans"),
        pval=c(1e-12, 1e-42, 1e-4),
        stringsAsFactors=F)
    
    g <- GAM:::graph.from.tables(node.table=list(met=met.nt), edge.table=et[et$rptype == "main", ], directed=F)
    
    es <- newEmptyObject()
    es$met.de.ext <- met.nt
    es$net.edges.ext.all <- et
    
    g1 <- addTransEdges(g, es)
    expect_true("R03" %in% E(g1)$rxn)
    
})