library("GAM.db")
library("GAM.networks")
data("examples")

data("kegg.mouse.network")

es.rn <- makeExperimentSet(network=kegg.mouse.network,
                           met.de=met.de.M0.M1,
                           gene.de=gene.de.M0.M1,
                           reactions.as.edges=FALSE,
                           plot=FALSE)
    
\dontrun{
heinz.py <- "/usr/local/lib/heinz/heinz.py"
solver <- heinz.solver(heinz.py)

module.rn <- findModule(es.rn,
                        met.fdr=1e-6,
                        rxn.fdr=1e-6,
                        absent.met.score=-20,
                        solver=solver)
}

module.rn <- addMetabolitesForReactions(module.rn, es.rn)
module.rn <- addInterconnections(module.rn, es.rn)
module.rn <- addNormLogFC(module.rn)
module.rn <- removeHangingNodes(module.rn)
module.rn <- simplifyReactionNodes(module.rn, es.rn)
module.rn <- expandReactionNodeAttributesToEdges(module.rn)
