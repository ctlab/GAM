
state1="MandLPSandIFNg"
state2="MandIL4"
fdr=c(1e-9, 1e-7, 1e-5, 1e-3)
fdr=1e-7

network.base="./networks/mouse1415/mouse1415.nogene.nocomp.masked.squared.noex.hmdb.esc"
rxn2genes.map.file="./networks//mouse1415/mouse1415.rxn2genes.tsv"

gene.exprs.file <- "./data_new_good/Gene.exp.tsv"
gene.conditions.file <- "./data_new_good/Gene.conditions.csv"

met.exprs.file <- "./data_new_good/MetDataR2.csv"
met.conditions.file <- "./data_new_good/Met.conditions.csv"

gene.exprs.sep=","; if (grepl("tsv$", gene.exprs.file)) { gene.exprs.sep <- "\t" }
gene.exprs <- read.csv(file=gene.exprs.file, head=TRUE, row.names=1, sep=gene.exprs.sep)
gene.conditions <- read.csv(file=gene.conditions.file, head=TRUE)
gene.conditions.vector <- gene.conditions$condition[match(colnames(gene.exprs), gene.conditions$sample)]

met.exprs.sep=","; if (grepl("tsv$", met.exprs.file)) { met.exprs.sep <- "\t" }
met.exprs <- read.csv(file=met.exprs.file, head=TRUE, row.names=1, sep=met.exprs.sep)
met.conditions <- read.csv(file=met.conditions.file, head=TRUE)
met.conditions.vector <- met.conditions$condition[match(colnames(met.exprs), met.conditions$sample)]

detach('package:GAM', unload=T)
library(GAM)

gene.de <- diff.expr(
    exprs=gene.exprs, conditions.vector=gene.conditions.vector,
    state1=state1, state2=state2, 
    log2=F, quantile=F, use.deseq=T, top=10000)

met.de <- diff.expr(
    exprs=met.exprs, conditions.vector=met.conditions.vector,
    state1=state1, state2=state2, 
    log2=F, quantile=F, use.deseq=F, top=Inf)


library(biomaRt)
mart = useMart("ensembl", dataset="mmusculus_gene_ensembl")
gene.de <- convert.pval.biomart(gene.de, from="refseq_mrna", to="entrezgene", mart=mart)

rxn2genes <- read.csv(rxn2genes.map.file, header= TRUE, check.names=F, sep="\t", colClasses="character")

rxn.de <- convert.pval(gene.de, from=rxn2genes$gene, to=rxn2genes$rxn)