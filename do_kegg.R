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
heinz.nModules <- 2

# If there is no heinz.py, heuristic search will be run instead.
#heinz.py <- NULL

# Pair of vectors of corresponding FDR values for metabolites and genes.
# Vectors have to be of the same length!
met.fdrs=c(1e-9, 1e-9, 1e-7, 1e-7, 1e-5)
gene.fdrs=c(1e-9, 1e-7, 1e-7, 1e-5, 1e-5)
 
# Getting list of most significant modules base both on gene and metabolite data
modules <- find_modules(network=kegg.mouse.network,
                        met.de=mouse.macrophages$met.de,
                        gene.de=mouse.macrophages$gene.de,
                        met.ids=mouse.macrophages$met.ids, 
                        gene.ids=mouse.macrophages$gene.ids,
                        met.fdrs=met.fdrs,
                        gene.fdrs=gene.fdrs,                         
                        heinz.py=heinz.py, 
                        heinz.nModules=heinz.nModules
                        )

# Saving pdf- and sif- files for found modules
for (module in modules) {
    save_module(module$graph, 
                paste0(outdir, "/module.", 
                       "mf=", module$met.fdr,
                       ".rf=", module$gene.fdr,
                       if (is.null(module$n)) "" else paste0("#", module$n)
                )
    )
}

# FDRs for gene-only and metabolite-only analysis
fdrs=c(1e-9, 1e-7, 1e-5, 1e-3)

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

