context("Solvers")


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
