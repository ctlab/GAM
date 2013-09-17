library(GAM)
gene.exprs.file <- "./data_new_good/Gene.exp.tsv"
gene.conditions.file <- "./data_new_good/Gene.conditions.csv"
gene.data.ids="refseq_mrna"

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

state0 <- "Ctrl"
state1 <- "MandLPSandIFNg"
state2 <- "MandIL4"

gene.de.M1.M2 <- diffExpr(
    exprs=gene.exprs, conditions.vector=gene.conditions.vector,
    state1=state1, state2=state2, 
    log2=F, quantile=F, use.deseq=T, top=20000)

met.de.M1.M2 <- diffExpr(
    exprs=met.exprs, conditions.vector=met.conditions.vector,
    state1=state1, state2=state2, 
    log2=F, quantile=F, use.deseq=F, top=Inf)

gene.de.M0.M1 <- diffExpr(
    exprs=gene.exprs, conditions.vector=gene.conditions.vector,
    state1=state0, state2=state1, 
    log2=F, quantile=F, use.deseq=T, top=20000)

met.de.M0.M1 <- diffExpr(
    exprs=met.exprs, conditions.vector=met.conditions.vector,
    state1=state0, state2=state1, 
    log2=F, quantile=F, use.deseq=F, top=Inf)

gene.de.M0.M2 <- diffExpr(
    exprs=gene.exprs, conditions.vector=gene.conditions.vector,
    state1=state0, state2=state2, 
    log2=F, quantile=F, use.deseq=T, top=20000)

met.de.M0.M2 <- diffExpr(
    exprs=met.exprs, conditions.vector=met.conditions.vector,
    state1=state0, state2=state2, 
    log2=F, quantile=F, use.deseq=F, top=Inf)

mouse.macrophages <- newEmptyObject()
mouse.macrophages$gene.exprs <- gene.exprs
mouse.macrophages$gene.conditions.vector <- gene.conditions.vector
mouse.macrophages$met.conditions.vector <- met.conditions.vector
mouse.macrophages$met.exprs <- met.exprs
mouse.macrophages$state0 <- state0
mouse.macrophages$state1 <- state1
mouse.macrophages$state2 <- state2
mouse.macrophages$gene.de.M1.M2 <- gene.de.M1.M2
mouse.macrophages$met.de.M1.M2 <- met.de.M1.M2
mouse.macrophages$gene.de.M0.M1 <- gene.de.M0.M1
mouse.macrophages$met.de.M0.M1 <- met.de.M0.M1
mouse.macrophages$gene.de.M0.M2 <- gene.de.M0.M2
mouse.macrophages$met.de.M0.M2 <- met.de.M0.M2
mouse.macrophages$met.ids <- "HMDB"
mouse.macrophages$gene.ids <- "RefSeq"
save(mouse.macrophages, file="./GAM/data/mouse.macrophages.rda")
