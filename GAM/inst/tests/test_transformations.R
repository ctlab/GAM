context("Transformations")

test_that("graphNEL.from.tables works", {
    met.nt <- data.frame(ID=c("C01", "C02"), mass=c(10, 20), name=c("Blah", "Blah"), stringsAsFactors=F)
    rxn.nt <- data.frame(ID=c("R01", "R02"), gene=c("Gene1", "Gene2"), stringsAsFactors=F)
    
    et <- data.frame(
        u=c("C01", "C01", "C02", "C02"), 
        v=c("R01", "R02", "R01", "R02"), 
        k=1:4,
        name=paste0("Connection", 1:4),
        stringsAsFactors=F)
    
    g <- graphNEL.from.tables(node.table=list(met=met.nt, rxn=rxn.nt), edge.table=et, directed=F)
    
    expect_equal(isDirected(g), F)
    expect_equal(length(nodes(g)), nrow(met.nt) + nrow(rxn.nt))
    expect_equal(nrow(getEdgeList(g)), nrow(et))
    
    expect_equal(nodeData(g, "C01", "nodeType")[[1]], "met")
    expect_equal(nodeData(g, "C02", "mass")[[1]], 20)
    expect_equal(nodeData(g, "C02", "label")[[1]], "Blah")
    
    expect_equal(nodeData(g, "R01", "gene")[[1]], "Gene1")
    expect_equal(nodeData(g, "R02", "nodeType")[[1]], "rxn")
    
    expect_equal(edgeData(g, from="C02", to="R01", "k")[[1]], 3)
    expect_equal(edgeData(g, from="C02", to="R01", "label")[[1]], "Connection3")
})