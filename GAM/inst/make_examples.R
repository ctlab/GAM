#!/usr/bin/env Rscript
library("devtools")
library("GAM.db")
load_all("GAM")
load_all("GAM.networks")

library("igraph")

gene.de.M0.M1 <- read.csv(text=getURL("https://artyomovlab.wustl.edu/publications/supp_materials/GAM_2015/Ctrl.vs.MandLPSandIFNg.gene.de.tsv"), sep="\t")
met.de.M0.M1 <- read.csv(text=getURL("https://artyomovlab.wustl.edu/publications/supp_materials/GAM_2015/Ctrl.vs.MandLPSandIFNg.met.de.tsv"), sep="\t")


data("kegg.mouse.network")

data(met.id.map)

mets <- unique(c(kegg.mouse.network$graph.raw$met.x, kegg.mouse.network$graph.raw$met.y))
mets.hmdb <- met.id.map$HMDB[met.id.map$KEGG %in% mets] 
rxns <- unique(kegg.mouse.network$graph.raw$rxn)
genes <- kegg.mouse.network$rxn2gene$gene[kegg.mouse.network$rxn2gene$rxn %in% rxns]
genes.refseq <- kegg.mouse.network$gene.id.map$RefSeq[kegg.mouse.network$gene.id.map$Entrez %in% genes] 

gene.de.M0.M1 <- gene.de.M0.M1[gene.de.M0.M1$ID %in% genes.refseq, ]
met.de.M0.M1 <- met.de.M0.M1[met.de.M0.M1$ID %in% mets.hmdb, ]

heinz.py <- "/usr/local/lib/heinz/heinz.py"
solver <- heinz.solver(heinz.py)

es.re <- makeExperimentSet(network=kegg.mouse.network,
                           met.de=met.de.M0.M1,
                           gene.de=gene.de.M0.M1,
                           reactions.as.edges=T)



met.fdr=c(3e-5)
rxn.fdr=c(3e-5)
absent.met.score=c(-20)

module.re <- findModule(es.re,
                        met.fdr=met.fdr,
                        rxn.fdr=rxn.fdr,
                        absent.met.score=absent.met.score,
                        solver=solver)
module.re

es.rn <- makeExperimentSet(network=kegg.mouse.network,
                           met.de=met.de.M0.M1,
                           gene.de=gene.de.M0.M1,
                           reactions.as.edges=F,
                           plot=F)

met.fdr=c(2e-7)
rxn.fdr=c(2e-7)
absent.met.score=c(-20)

solver2 <- heinz2.solver(heinz2="/usr/local/lib/heinz2/heinz", timeLimit=30)
module.rn <- findModule(es.rn,
                        met.fdr=met.fdr,
                        rxn.fdr=rxn.fdr,
                        absent.met.score=absent.met.score,
                        solver=solver2)
module.rn

save(gene.de.M0.M1, met.de.M0.M1, es.re, es.rn, module.re, module.rn, file="examplesGAM.rda")
