#!/usr/bin/env Rscript

set.seed(42)
library("GAM")

# Loading network
data("kegg.mouse.network")

# Loading dataset for mouse macrophages cells
data("mouse.macrophages")

# The dataset consists of gene expression data,
str(mouse.macrophages$gene.exprs)

# conditions vector for gene,
str(mouse.macrophages$gene.conditions.vector)

# metabolite mass-spectrometry data,
str(mouse.macrophages$met.exprs)

# conditions vector for metabolites,
str(mouse.macrophages$met.conditions.vector)

# gene and metabolite differential expression results
str(mouse.macrophages$gene.de)
str(mouse.macrophages$met.de)

# between states
mouse.macrophages$state1
mouse.macrophages$state2

# Results will be put into directory.
outdir=paste0("./kegg_mouse_", mouse.macrophages$state1, "-", mouse.macrophages$state2)
outdir

# If heinz.py is available and working we can set its path.
heinz.py <- "./heinz.py"; 
# And number of modules we'd like to find.
heinz.nModules <- 3

# If there is no heinz.py, heuristic search will be run instead.
#heinz.py <- NULL

es.M1.M2.raw <- make_experiment_set.raw(network=kegg.mouse.network,
                                met.de=mouse.macrophages$met.de,
                                gene.de=mouse.macrophages$gene.de,
                                met.ids=mouse.macrophages$met.ids, 
                                gene.ids=mouse.macrophages$gene.ids)

# Pair of vectors of corresponding FDR values for metabolites and genes.
# Vectors have to be of the same length!
met.fdrs  <- c(1e-9, 1e-9, 1e-7, 1e-7, 1e-5, 1e-5, 1e-5, 1e-5, 1e-3, 1e-3)
gene.fdrs <- c(1e-9, 1e-7, 1e-7, 1e-5, 1e-7, 1e-6, 1e-5, 1e-3, 1e-3, 1e-2)



met.fdrs=c(1e-2)
gene.fdrs=c(1e-3)
met.fdrs=c(1e-5)
gene.fdrs=c(1e-5)

# Getting list of most significant modules base both on gene and metabolite data
modules <- find_modules.raw(es.M1.M2.raw,
                        met.fdrs=met.fdrs,
                        gene.fdrs=gene.fdrs,                         
                        heinz.py=heinz.py, 
                        heinz.nModules=3,
                        score.separately=T
                        )


set.seed(42)
# Saving pdf- and sif- files for found modules
for (module in modules) {
    save_module(module$graph, 
                paste0(outdir, "/module.raw.", 
                       "mf=", module$met.fdr,
                       ".rf=", module$gene.fdr,
                       if (is.null(module$n)) "" else paste0("#", module$n)
                )
    )
}
exit(0)

state0 <- "Ctrl"

gene.de.M0.M1 <- diff.expr(
    exprs=mouse.macrophages$gene.exprs, conditions.vector=mouse.macrophages$gene.conditions.vector,
    state1=state0, state2=mouse.macrophages$state1, 
    log2=F, quantile=F, use.deseq=T, top=10000)
gene.de.M0.M1 <- convert.pval(gene.de.M0.M1, from=gene.id.map$RefSeq, to=gene.id.map$Entrez)
gene.de.M0.M1$logFC <- fix_inf(gene.de.M0.M1$logFC)

met.de.M0.M1 <- diff.expr(
    exprs=mouse.macrophages$met.exprs, conditions.vector=mouse.macrophages$met.conditions.vector,
    state1=state0, state2=mouse.macrophages$state1, 
    log2=F, quantile=F, use.deseq=F, top=Inf)
met.de.M0.M1 <- convert.pval(met.de.M0.M1, from=met.id.map$HMDB, to=met.id.map$KEGG)
all.de.M0.M1 <- rbind(gene.de.M0.M1, met.de.M0.M1)


gene.de.M0.M2 <- diff.expr(
    exprs=mouse.macrophages$gene.exprs, conditions.vector=mouse.macrophages$gene.conditions.vector,
    state1=state0, state2=mouse.macrophages$state2, 
    log2=F, quantile=F, use.deseq=T, top=10000)
gene.de.M0.M2 <- convert.pval(gene.de.M0.M2, from=gene.id.map$RefSeq, to=gene.id.map$Entrez)
gene.de.M0.M2$logFC <- fix_inf(gene.de.M0.M2$logFC)

met.de.M0.M2 <- diff.expr(
    exprs=mouse.macrophages$met.exprs, conditions.vector=mouse.macrophages$met.conditions.vector,
    state1=state0, state2=mouse.macrophages$state2, 
    log2=F, quantile=F, use.deseq=F, top=Inf)
met.de.M0.M2 <- convert.pval(met.de.M0.M2, from=met.id.map$HMDB, to=met.id.map$KEGG)
all.de.M0.M2 <- rbind(gene.de.M0.M2, met.de.M0.M2)

set.seed(42)
for (module in modules) {    
    nodeData(module$graph, attr="diff.expr") <-
        all.de.M0.M1$logFC[match(nodes(module$graph), all.de.M0.M1$ID)]
    unmatched <- nodes(module$graph)[(!nodes(module$graph) %in% all.de.M0.M1$ID)]
    nodeData(module$graph, unmatched, attr="diff.expr") <- rep(0, length(unmatched))
    module$graph <- add_norm.diff.expr(module$graph)
    save_module(module$graph, 
                paste0(outdir, "/module.M0.M1", 
                       "mf=", module$met.fdr,
                       ".rf=", module$gene.fdr,
                       if (is.null(module$n)) "" else paste0("#", module$n)
                )
    )
}

set.seed(42)
for (module in modules) {    
    nodeData(module$graph, attr="diff.expr") <-
        all.de.M0.M2$logFC[match(nodes(module$graph), all.de.M0.M2$ID)]
    unmatched <- nodes(module$graph)[(!nodes(module$graph) %in% all.de.M0.M2$ID)]
    nodeData(module$graph, unmatched, attr="diff.expr") <- rep(0, length(unmatched))
    module$graph <- add_norm.diff.expr(module$graph)
    save_module(module$graph, 
                paste0(outdir, "/module.M0.M2", 
                       "mf=", module$met.fdr,
                       ".rf=", module$gene.fdr,
                       if (is.null(module$n)) "" else paste0("#", module$n)
                )
    )
}
outdir.M0.M1="./kegg_mouse_M0-M1"
modules.M0.M1 <- find_modules(network=kegg.mouse.network,
                        met.de.M0.M1,
                        gene.de.M0.M1,                        
                        met.fdrs=met.fdrs,
                        gene.fdrs=gene.fdrs,                         
                        heinz.py=heinz.py, 
                        heinz.nModules=heinz.nModules)
for (module in modules.M0.M1) {
    save_module(module$graph, 
                paste0(outdir.M0.M1, "/module.", 
                       "mf=", module$met.fdr,
                       ".rf=", module$gene.fdr,
                       if (is.null(module$n)) "" else paste0("#", module$n)
                )
    )
}

outdir.M0.M2="./kegg_mouse_M0-M2"
modules.M0.M2 <- find_modules(network=kegg.mouse.network,
                              met.de.M0.M2,
                              gene.de.M0.M2,                        
                              met.fdrs=met.fdrs,
                              gene.fdrs=gene.fdrs,                         
                              heinz.py=heinz.py, 
                              heinz.nModules=heinz.nModules)
for (module in modules.M0.M2) {
    save_module(module$graph, 
                paste0(outdir.M0.M2, "/module.", 
                       "mf=", module$met.fdr,
                       ".rf=", module$gene.fdr,
                       if (is.null(module$n)) "" else paste0("#", module$n)
                )
    )
}                        


met.fdrs=1e-9
gene.fdrs=1e-9

for (subopt_diff in c(0:5) * 20) {
    print(subopt_diff)
    modules <- find_modules(network=kegg.mouse.network,
                            met.de=mouse.macrophages$met.de,
                            gene.de=mouse.macrophages$gene.de,
                            met.ids=mouse.macrophages$met.ids, 
                            gene.ids=mouse.macrophages$gene.ids,
                            met.fdrs=met.fdrs,
                            gene.fdrs=gene.fdrs,                         
                            heinz.py=heinz.py, 
                            heinz.nModules=3,
                            heinz.subopt_diff=subopt_diff                            
    )
    
    
    
    # Saving pdf- and sif- files for found modules
    for (module in modules) {
        save_module(module$graph, 
                    paste0(outdir, "/module.", 
                           "sd=", subopt_diff,
                           ".mf=", module$met.fdr,
                           ".rf=", module$gene.fdr,
                           if (is.null(module$n)) "" else paste0("#", module$n)
                    )
        )
    }
    
}

exit(0)
# FDRs for gene-only and metabolite-only analysis
fdrs=c(1e-9, 1e-7, 1e-5, 1e-3, 1e-2)

# Finding metabolite modules. Just leaving gene data out of arguments.
modules.mets <- find_modules(network=kegg.mouse.network,
                        met.de=mouse.macrophages$met.de,                        
                        met.ids=mouse.macrophages$met.ids,                         
                        fdrs=fdrs,
                        heinz.py=heinz.py, 
                        heinz.nModules=heinz.nModules
)

# Saving them.
for (module in modules.mets) {
    save_module(module$graph, 
                paste0(outdir, "/module.mets.", 
                       "fdr=", module$met.fdr,
                       if (is.null(module$n)) "" else paste0("#", module$n)
                )
    )
}

# # Finding gene modules. Vice versa, leaving metabolite data out of arguments.
modules.genes <- find_modules(network=kegg.mouse.network,                        
                        gene.de=mouse.macrophages$gene.de,                        
                        gene.ids=mouse.macrophages$gene.ids,
                        fdrs=fdrs,                                 
                        heinz.py=heinz.py, 
                        heinz.nModules=heinz.nModules
)

# Saving them.
for (module in modules.genes) {
    save_module(module$graph, 
                paste0(outdir, "/module.genes.", 
                       "fdr=", module$gene.fdr,                       
                       if (is.null(module$n)) "" else paste0("#", module$n)
                )
    )
}

