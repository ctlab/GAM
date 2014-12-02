context("Network analysis")

heinz.py <- "/usr/local/lib/heinz/heinz.py"
mwcs <- "/usr/local/bin/mwcs"

if (file.exists(heinz.py)) {
    test_that("runHeinz works", {
        met.nt <- data.frame(ID=c("C01", "C02", "C03"), score=c(1, 1, 1), stringsAsFactors=F)
        
        et <- data.frame(
            met.x=c("C01", "C02"),
            met.y=c("C02", "C03"),
            rxn=c("R01", "R02"),
            score=c(1, -10),
            stringsAsFactors=F)
        
        g <- GAM:::graph.from.tables(node.table=list(met=met.nt), edge.table=et, directed=F)
        
        module <- GAM:::runHeinz(g, heinz.py, score.nodes=T, score.edges=T)[[1]]
        
        expect_equivalent(V(module)$name, c("C01", "C02"))
    })
}

data("kegg.mouse.network")
library(mouseMacrophages)
data(examplesGAM)
data(mmpData)
library(igraph)

# :ToDo: check the result for makeExperimentSet tests

test_that("makeExperimentSet works with full data", {
    es.M1.M2.full.rn.cr <- makeExperimentSet(network=kegg.mouse.network, 
                                  met.de=met.de.M1.M2,
                                  gene.de=gene.de.M1.M2,
                                  reactions.as.edges=F, collapse.reactions=T, plot=F)
    
    es.M1.M2.full.rn <- makeExperimentSet(network=kegg.mouse.network, 
                                  met.de=met.de.M1.M2,
                                  gene.de=gene.de.M1.M2,
                                  reactions.as.edges=F, collapse.reactions=F, plot=F)
    
    es.M1.M2.full.re.rp <- makeExperimentSet(network=kegg.mouse.network,
                                  met.de=met.de.M1.M2,
                                  gene.de=gene.de.M1.M2,
                                  reactions.as.edges=T, use.rpairs=T, plot=F)
    
    es.M1.M2.full.re <- makeExperimentSet(network=kegg.mouse.network,
                                  met.de=met.de.M1.M2,
                                  gene.de=gene.de.M1.M2,
                                  reactions.as.edges=T, use.rpairs=F, plot=F)
})

test_that("makeExperimentSet works without metabolic data", {
    es.M1.M2 <- makeExperimentSet(network=kegg.mouse.network, 
                                  gene.de=gene.de.M1.M2,
                                  reactions.as.edges=F, collapse.reactions=T, plot=F)

    es.M1.M2 <- makeExperimentSet(network=kegg.mouse.network, 
                                  gene.de=gene.de.M1.M2,
                                  reactions.as.edges=F, collapse.reactions=F, plot=F)
    
    es.M1.M2 <- makeExperimentSet(network=kegg.mouse.network,
                                  gene.de=gene.de.M1.M2,
                                  reactions.as.edges=T, use.rpairs=T, plot=F)
    
    es.M1.M2 <- makeExperimentSet(network=kegg.mouse.network,
                                  gene.de=gene.de.M1.M2,
                                  reactions.as.edges=T, use.rpairs=F, plot=F)
})

test_that("makeExperimentSet works without genomic data", {
    es.M1.M2 <- makeExperimentSet(network=kegg.mouse.network, 
                                  met.de=met.de.M1.M2,
                                  reactions.as.edges=F, plot=F)
    
    expect_true(length(E(es.M1.M2$subnet)[adj("C05528")]) > 0)
    
    es.M1.M2 <- makeExperimentSet(network=kegg.mouse.network,
                                  met.de=met.de.M1.M2,
                                  reactions.as.edges=T, use.rpairs=T, plot=F)
    
    es.M1.M2 <- makeExperimentSet(network=kegg.mouse.network,
                                  met.de=met.de.M1.M2,
                                  reactions.as.edges=T, use.rpairs=F, plot=F)
})

test_that("findModule works", { # :ToDo:
})

test_that("makeKeggNetwork works", { # :ToDo:
})

test_that("scoreNetwork works", {
    es.M1.M2.scored <- scoreNetwork(es.M1.M2.full.rn.cr)
    es.M1.M2.scored <- scoreNetwork(es.M1.M2.full.rn)
    es.M1.M2.full.re.rp.scored <- scoreNetwork(es.M1.M2.full.re.rp)
    expect_true("score" %in% list.edge.attributes(es.M1.M2.full.re.rp.scored$subnet.scored))
    expect_true("score" %in% list.vertex.attributes(es.M1.M2.full.re.rp.scored$subnet.scored))
    es.M1.M2.scored <- scoreNetwork(es.M1.M2.full.re)
    expect_true("score" %in% list.edge.attributes(es.M1.M2.scored$subnet.scored))
})

test_that("scoreNetworkWithoutBUM works", { # :ToDo:
    es.M1.M2.scored <- scoreNetworkWithoutBUM(es.M1.M2.full.rn.cr)
    es.M1.M2.scored <- scoreNetworkWithoutBUM(es.M1.M2.full.rn)
    es.M1.M2.scored <- scoreNetworkWithoutBUM(es.M1.M2.full.re.rp)
    es.M1.M2.scored <- scoreNetworkWithoutBUM(es.M1.M2.full.re)
    
})

test_that("fastHeinz.solver works", { # :ToDo:
})

test_that("heinz.solver works", { # :ToDo:
})

test_that("mwcs.solver works", { # :ToDo:
})

test_that("randHeur.solver works", { # :ToDo:
})
