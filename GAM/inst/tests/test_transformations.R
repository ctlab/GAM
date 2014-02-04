library(igraph)

context("Graph transformations")

test_that("graph.from.tables works", {
    met.nt <- data.frame(ID=c("C01", "C02"), mass=c(10, 20), name=c("Blah", "Blah"), stringsAsFactors=F)
    rxn.nt <- data.frame(ID=c("R01", "R02"), gene=c("Gene1", "Gene2"), stringsAsFactors=F)
    
    et <- data.frame(
        u=c("C01", "C01", "C02", "C02"), 
        v=c("R01", "R02", "R01", "R02"), 
        k=1:4,
        name=paste0("Connection", 1:4),
        stringsAsFactors=F)
    
    g <- GAM:::graph.from.tables(node.table=list(met=met.nt, rxn=rxn.nt), edge.table=et, directed=F)
    
    expect_equal(is.directed(g), F)
    expect_equal(length(V(g)), nrow(met.nt) + nrow(rxn.nt))
    expect_equal(length(E(g)), nrow(et))
    
    expect_equal(V(g)["C01"]$nodeType, "met")
    expect_equal(V(g)["C02"]$mass, 20)
    expect_equal(V(g)["C02"]$label, "Blah")
    
    expect_equal(V(g)["R01"]$gene, "Gene1")
    expect_equal(V(g)["R02"]$nodeType, "rxn")
    
    expect_equal(E(g, P=c("C02", "R01"))$k, 3)
    expect_equal(E(g, P=c("C02", "R01"))$label, "Connection3")
})

test_that("module2list works", {
    met.nt <- data.frame(
        ID=c( "C01", "C02", "C03"), 
        pval=c(1e-5,    NA,  1e-6), stringsAsFactors=F)
    et <- data.frame(
        met.x=c("C01", "C02", "C01"),
        met.y=c("C02", "C03", "C03"),
        rxn=c("R01", "R02", "R03"),
        rptype=c("main", "main", "trans"),
        pval=c(1e-12, 1e-42, 1e-4),
        stringsAsFactors=F)
    
    g <- GAM:::graph.from.tables(node.table=list(met=met.nt), edge.table=et[et$rptype == "main", ], directed=F)
    l <- module2list(g)
    expect_equal(l$nodes[[1]]$pval, 1e-5)
    expect_equal(l$edges[[1]]$rxn, "R01")
})

test_that("plotNetwork works", { # :ToDo:
})


test_that("get.edge.attributes works", { # :ToDo:
})

test_that("get.vertex.attributes works", { # :ToDo:
})
