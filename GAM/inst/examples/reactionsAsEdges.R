data("kegg.mouse.network")
data("examples")
library("mouseMacrophages")
data("mmpData")

heinz.py <- "/usr/local/lib/heinz/heinz.py"
solver <- heinz.solver(heinz.py)

es.re <- makeExperimentSet(network=kegg.mouse.network,
                           met.de=met.de.M1.M2,
                           gene.de=gene.de.M1.M2,
                           reactions.as.edges=TRUE)

\dontrun{
module.re <- findModule(es.re,
                        met.fdr=3e-05,
                        rxn.fdr=3e-05,
                        absent.met.score=-20,
                        solver=solver)
}

module.re <- addTransEdges(module.re, es.re)

plotNetwork(module.re)

saveModule(module.re,
           paste0("module.M1.M2.re", 
                  ".mf=", met.fdr,
                  ".rf=", rxn.fdr,
                  ".ms=", absent.met.score),
           types=c("pdf", "XGMML")