#!/usr/bin/env Rscript

set.seed(42)
state1="MandLPSandIFNg"
state2="MandIL4"
outdir=paste("./kegg_", state1, "-", state2, sep="")

met.fdrs=c(1e-9, 1e-9, 1e-7, 1e-7, 1e-5, 1e-3)
rxn.fdrs=c(1e-7, 1e-5, 1e-5, 1e-3, 1e-3, 1e-3)

network.base="./networks//kegg/net.sq"

gene.network.ids="entrezgene"
network.ensembl.dataset="mmusculus_gene_ensembl"

gene.exprs.file <- "./data_new_good/Gene.exp.tsv"
gene.conditions.file <- "./data_new_good/Gene.conditions.csv"
gene.data.ids="refseq_mrna"

met.exprs.file <- "./data_new_good/MetDataR2.csv"
met.conditions.file <- "./data_new_good/Met.conditions.csv"

#nModules=1; heinz.py=NULL
nModules=2; heinz.py="./heinz.py"

met.ids <- read.table("./misc/metabolite_ids.tsv", header=T, colClasses="character")

#detach('package:GAM', unload=T)
library(GAM)

enz2gene <- read.table("./networks//kegg/enz2gene.tsv", header=T, colClasses="character")
rxn2enz <- read.table("./networks//kegg/rxn2enz.tsv", header=T, colClasses="character")
rxn2gene = merge(rxn2enz, enz2gene)[, c("rxn", "gene")]




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
    gene.de$origin <- NULL
}    

enz.de <- convert.pval(gene.de, from=enz2gene$gene, to=enz2gene$enz)
met.de <- convert.pval(met.de, from=met.ids$HMDB, to=met.ids$KEGG)
rxn.de <- convert.pval(enz.de, from=rxn2enz$enz, to=rxn2enz$rxn)


library(BioNet)
network <- loadNetwork.sif(
    paste(network.base, "sif", sep="."),
    list.files(dirname(network.base), paste(basename(network.base), "_\\w+.NA", sep=""), full.names=T)
)



edgelist <- function(network) {
    edges <- edges(network)
    vs <- unlist(edges)
    us <- rep(names(edges), sapply(edges, length))
    res <- data.frame(u=us, v=vs, stringsAsFactors=F)
    return(res)
}

convert.node.names <- function(network, from, to) {
    edges <- edgelist(network)
    m1 <- match(edges$u, from)
    edges$u[!is.na(m1)] <- to[m1[!is.na(m1)]]
    m2 <- match(edges$v, from)
    edges$v[!is.na(m2)] <- to[m2[!is.na(m2)]]
    edges <- unique(edges)
    network2 <- graph.edgelist(as.matrix(edges), directed=F)
    network2 <- simplify(network2, remove.multiple=T)
    network2 <- igraph.to.graphNEL(network2)
    
    nodeDataDefaults(network2) <- nodeDataDefaults(network)
    
    old_nodes <- intersect(nodes(network), nodes(network2))
    
    for (attr in names(nodeDataDefaults(network))) {
        nodeData(network2, old_nodes, attr) <- nodeData(network, old_nodes, attr)
    }
    
    network2   
}

network <- convert.node.names(network, rxn.de$ID, rxn.de$origin)
rxn.de.orig <- rxn.de
rxn.de$ID <- rxn.de$origin

reflink <- read.csv("../reflink.txt", sep="\t", header=T, colClasses="character")
rxn.de <- rxn.de[rxn.de$ID %in% nodes(network), ]
nodeData(network, rxn.de$ID, "shortName") <- reflink$name[match(rxn.de$origin, reflink$locusLinkId)]


dir.create(outdir)
modules <- find_modules(met.de=met.de, rxn.de=rxn.de,
                        network=network, 
                        met.fdrs=met.fdrs, rxn.fdrs=rxn.fdrs, 
                        nModules=nModules, heinz.py=heinz.py)
for (module in modules) {
    save_module(module$graph, 
                paste0(outdir, "/module.", 
                       "mf=", module$met.fdr,
                       ".rf=", module$rxn.fdr,
                       if (is.null(module$n)) "" else paste0("#", module$n)
                )
    )
}
