library("mouseMacrophages")
data("mmpData")

data("kegg.db")
library("org.Mm.eg.db")
kegg.mouse.network <- makeKeggNetwork(kegg.db, "MMU")

state1 <- "MandLPSandIFNg"
state2 <- "MandIL4"

gene.de.M1.M2 <- diffExpr(
    exprs=exprs(mmpGeneSet), conditions.vector=pData(mmpGeneSet)$condition,
    state1=state1, state2=state2, 
    use.deseq=TRUE, top=Inf, min.expr=5)
met.de.M1.M2 <- diffExpr(
    exprs=exprs(mmpMetSet), conditions.vector=pData(mmpMetSet)$condition,
    state1=state1, state2=state2, 
    log2=FALSE, quantile=FALSE, use.deseq=FALSE, top=Inf)


heinz.py <- "/usr/local/lib/heinz/heinz.py"
solver <- heinz.solver(heinz.py)

es.rn <- makeExperimentSet(network=kegg.mouse.network,
                           met.de=met.de.M1.M2,
                           gene.de=gene.de.M1.M2,
                           reactions.as.edges=FALSE,
                           plot=FALSE)

module.rn <- findModule(es.rn,
                        met.fdr=1e-6,
                        gene.fdr=1e-6,
                        absent.met.score=-20,
                        solver=solver)

module.rn <- addMetabolitesForReactions(module.rn, es.rn)
module.rn <- addInterconnections(module.rn, es.rn)
module.rn <- addNormLogFC(module.rn)
module.rn <- removeHangingNodes(module.rn)
module.rn <- removeSimpleReactions(module.rn, es.rn)
module.rn <- expandReactionNodeAttributesToEdges(module.rn)

plotNetwork(module.rn)