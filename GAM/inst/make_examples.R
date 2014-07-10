library("GAM")
library("igraph")
library("mouseMacrophages")
data("mmpData")

state1 <- "MandLPSandIFNg"
state2 <- "MandIL4"
gene.de.M1.M2 <- diffExpr(
    exprs=exprs(mmpGeneSet), conditions.vector=pData(mmpGeneSet)$condition,
    state1=state1, state2=state2, 
    use.deseq=T, top=Inf, min.expr=5)
met.de.M1.M2 <- diffExpr(
    exprs=exprs(mmpMetSet), conditions.vector=pData(mmpMetSet)$condition,
    state1=state1, state2=state2, 
    log2=F, quantile=F, use.deseq=F, top=Inf)
head(gene.de.M1.M2)
head(met.de.M1.M2)


data("kegg.mouse.network")

data(met.id.map)

mets <- unique(c(kegg.mouse.network$graph.raw$met.x, kegg.mouse.network$graph.raw$met.y))
mets.hmdb <- met.id.map$HMDB[met.id.map$KEGG %in% mets] 
rxns <- unique(kegg.mouse.network$graph.raw$rxn)
genes <- kegg.mouse.network$rxn2gene$gene[kegg.mouse.network$rxn2gene$rxn %in% rxns]
genes.refseq <- kegg.mouse.network$gene.id.map$RefSeq[kegg.mouse.network$gene.id.map$Entrez %in% genes] 

gene.de.M1.M2 <- gene.de.M1.M2[gene.de.M1.M2$ID %in% genes.refseq, ]
met.de.M1.M2 <- met.de.M1.M2[met.de.M1.M2$ID %in% mets.hmdb, ]

heinz.py <- "/usr/local/lib/heinz/heinz.py"
solver <- heinz.solver(heinz.py)

es.re <- makeExperimentSet(network=kegg.mouse.network,
                           met.de=met.de.M1.M2,
                           gene.de=gene.de.M1.M2,
                           reactions.as.edges=T)



met.fdr=c(3e-5)
gene.fdr=c(3e-5)
absent.met.score=c(-20)

module.re <- findModule(es.re,
                        met.fdr=met.fdr,
                        gene.fdr=gene.fdr,
                        absent.met.score=absent.met.score,
                        solver=solver)
module.re

es.rn <- makeExperimentSet(network=kegg.mouse.network,
                           met.de=met.de.M1.M2,
                           gene.de=gene.de.M1.M2,
                           reactions.as.edges=F,
                           plot=F)

met.fdr=c(2e-7)
gene.fdr=c(2e-7)
absent.met.score=c(-20)

module.rn <- findModule(es.rn,
                        met.fdr=met.fdr,
                        gene.fdr=gene.fdr,
                        absent.met.score=absent.met.score,
                        solver=solver)
module.rn

save(gene.de.M1.M2, met.de.M1.M2, module.re, module.rn, file="examplesGAM.rda")
