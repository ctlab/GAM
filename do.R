#!/usr/bin/env Rscript

state1="MandLPSandIFNg"
state2="MandIL4"
outdir=paste("./data_new1_", state1, "-", state2, sep="")

fdrs=c(1e-9, 1e-7, 1e-5, 1e-3)

network.base="./networks/mouse1415/mouse1415.nogene.nocomp.masked.squared.noex.hmdb.esc"
rxn2genes.map.file="./networks//mouse1415/mouse1415.rxn2genes.tsv"
gene.network.ids="entrezgene"
network.ensembl.dataset="mmusculus_gene_ensembl"

gene.exprs.file <- "./data_new_good/Gene.exp.tsv"
gene.conditions.file <- "./data_new_good/Gene.conditions.csv"
gene.data.ids="refseq_mrna"

met.exprs.file <- "./data_new_good/MetDataR2.csv"
met.conditions.file <- "./data_new_good/Met.conditions.csv"


#detach('package:GAM', unload=T)
library(GAM)

modules <- do.analysis(
    network.base=network.base,
    rxn2genes.map.file=rxn2genes.map.file,
    gene.exprs.file=gene.exprs.file,
    gene.conditions.file=gene.conditions.file,
    met.exrps.file=met.exprs.file,
    met.conditions.file=met.conditions.file, 
    state1=state1, state2=state2,
    fdrs=fdrs,
    gene.data.ids=gene.data.ids, 
    gene.network.ids=gene.network.ids,
    network.ensembl.dataset=network.ensembl.dataset,
    #nModules=2, heinz.py="./heinz.py"
    nModules=1, heinz.py=NULL
    )

gene.exprs.sep=","; if (grepl("tsv$", gene.exprs.file)) { gene.exprs.sep <- "\t" }
gene.exprs <- read.csv(file=gene.exprs.file, head=TRUE, row.names=1, sep=gene.exprs.sep)
gene.conditions <- read.csv(file=gene.conditions.file, head=TRUE)
gene.conditions.vector <- gene.conditions$condition[match(colnames(gene.exprs), gene.conditions$sample)]

met.exprs.sep=","; if (grepl("tsv$", met.exprs.file)) { met.exprs.sep <- "\t" }
met.exprs <- read.csv(file=met.exprs.file, head=TRUE, row.names=1, sep=met.exprs.sep)
met.conditions <- read.csv(file=met.conditions.file, head=TRUE)
met.conditions.vector <- met.conditions$condition[match(colnames(met.exprs), met.conditions$sample)]

gene.de <- diff.expr(
    exprs=gene.exprs, conditions.vector=gene.conditions.vector,
    state1=state1, state2=state2, 
    log2=F, quantile=F, use.deseq=T, top=10000)

met.de <- diff.expr(
    exprs=met.exprs, conditions.vector=met.conditions.vector,
    state1=state1, state2=state2, 
    log2=F, quantile=F, use.deseq=F, top=Inf)


if (!is.null(network.ensembl.dataset) && !is.null(gene.data.ids) && !is.null(gene.network.ids)) {
    library(biomaRt)
    mart = useMart("ensembl", dataset=network.ensembl.dataset)
    gene.de <- convert.pval.biomart(gene.de, from=gene.data.ids, to=gene.network.ids, mart=mart)    
}    

rxn2genes <- read.csv(rxn2genes.map.file, header= TRUE, check.names=F, sep="\t", colClasses="character")

rxn.de <- convert.pval(gene.de, from=rxn2genes$gene, to=rxn2genes$rxn)

combined.de <- rbind(met.de, rxn.de)

library(BioNet)
network <- loadNetwork.sif(
    paste(network.base, "sif", sep="."),
    list.files(dirname(network.base), paste(basename(network.base), "_\\w+.NA", sep=""), full.names=T)
)

dir.create(outdir)
for (fdr in fdrs) {
    modules <- find_modules(data.pval=combined.de, network=network, nModules=nModules, fdrs=fdrs, heinz.py=heinz.py)
    for (module in modules) {
        save_module(module$graph, 
                    paste0(outdir, "/module.", 
                           module$fdr, 
                           if (is.null(module$n)) "" else paste0("#", module$n)
                    )
        )
    }
}
