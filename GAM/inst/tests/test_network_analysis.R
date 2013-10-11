context("Network analysis")

heinz.py <- "/usr/local/lib/heinz/heinz.py"

if (file.exists(heinz.py)) {
    test_that("runHeinz works", {
        met.nt <- data.frame(ID=c("C01", "C02", "C03"), score=c(1, 1, 1), stringsAsFactors=F)
        
        et <- data.frame(
            met.x=c("C01", "C02"),
            met.y=c("C02", "C03"),
            rxn=c("R01", "R02"),
            score=c(1, -10),
            stringsAsFactors=F)
        
        g <- graphNEL.from.tables(node.table=list(met=met.nt), edge.table=et, directed=F)
        
        module <- runHeinz(g, heinz.py, score.nodes=T, score.edges=T)[[1]]
        
        expect_equivalent(nodes(module), c("C01", "C02"))
    })
}