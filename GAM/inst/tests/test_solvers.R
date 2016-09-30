context("Solvers")

library(data.table)
library(BioNet)
library(igraph)

gmwcs <- "/usr/local/bin/gmwcs"

test_that("gmwcs.solver support one-node output", {
    if (!file.exists(gmwcs)) {
        skip()
    }
    
    met.nt <- data.frame(ID=c("C01", "C02", "C03"), score=c(1, -1, -1), stringsAsFactors=F)
    
    et <- data.frame(
        met.x=c("C01", "C02"),
        met.y=c("C02", "C03"),
        rxn=c("R01", "R02"),
        score=c(0, 0),
        stringsAsFactors=F)
    
    g <- GAM:::graph.from.tables(node.table=list(met=met.nt), edge.table=et, directed=F)
    solve <- gmwcs.solver(gmwcs, 1)
    
    solve(g)    
})


df.from.rows <- function(names, ...) {    
    rows <- list(...)
    res <- do.call(data.frame, 
               c(stringsAsFactors=F,
                   lapply(seq_along(names), 
                          function(i) sapply(rows, function(row) row[[i]]))))
    colnames(res) <- names
    res
}

# R1: C1 + C2 <=> C3 + C4  (G1)
# R2: C3 <=> C5   (G2)
# R3: C4 <=> C6   (G3)

names(kegg.mouse.network)

test.mouse.network  <- newEmptyObject()
test.mouse.network$rxn2gene <- data.frame(stringsAsFactors=FALSE, 
    rxn=c("R1", "R2", "R3"),
    gene=c("G1", "G2", "G3"))

test.mouse.network$rxn2name <- df.from.rows(
    c("rxn", "name", "pathway"),
    c("R1", "reaction1", "rn01"),
    c("R2", "reaction2", "rn01"),
    c("R3", "reaction3", "rn01"))
    

test.mouse.network$met2name <- 
    data.frame(
        stringsAsFactors=F,
        met=paste0("C", 1:6),
        name=paste0("compound", 1:6),
        pathway="map00")

test.mouse.network$graph.raw <- data.frame(stringsAsFactors=FALSE, rbind(
    c(met.x="C1", rxn="R1", met.y="C3", rpair=NA, rptype=NA),
    c("C2", "R1", "C3", NA, NA),
    c("C1", "R1", "C4", NA, NA),
    c("C2", "R1", "C4", NA, NA),
    c("C3", "R2", "C5", NA, NA),
    c("C4", "R3", "C6", NA, NA)))

test.mouse.network$gene.id.map <- data.table(stringsAsFactors=FALSE, rbind(
    c("Entrez"="G1", "Symbol"="gene1"),
    c("G2", "gene2"),
    c("G3", "gene3")))

test.mouse.network$gene.ids <- "Entrez"

test.gene.de <- df.from.rows(
    c("ID", "pval", "log2FC", "baseMean"),
    list("G1", 1e-3, -0.8, 100),
    list("G2", 0.5, 0.1, 200),
    list("G3", 1e-5, -1.4, 459))
    

test_that("randHeur works with genes-only and reactions as nodes", {
    es <- makeExperimentSet(test.mouse.network, gene.de=test.gene.de, reactions.as.edges=F)
    es.scored <- scoreNetwork(es, rxn.fdr=1e-4)
    m <- findModule(es.scored, solver=randHeur.solver(2))    
    expect_true(length(V(m)) == 1)
    
    es.scored <- scoreNetwork(es, rxn.fdr=1e-2)
    set.seed(42)
    m <- findModule(es.scored, solver=randHeur.solver(2))    
    expect_true(length(V(m)) == 3)
})
