#' Get KEGG network for specified organism
#' @param kegg.db object with KEGG mappings
#' @param organism Organism
#' @return Network for the organism usable for an analysis
#' @export
makeKeggNetwork <- function(kegg.db, organism) {
    enz2gene <- kegg.db$enz2gene[kegg.db$enz2gene$organism == organism, c("enz", "gene")]
    rxn2gene <- merge(kegg.db$rxn2enz, enz2gene)[, c("rxn", "gene")]
    
    rxns2keep <- unique(rxn2gene$rxn)
    rxns2keep <- setdiff(rxns2keep, kegg.db$rxns2mask)
    
    net <- kegg.db$net
    m <- match(net$met.x, kegg.db$mets2collapse$from)
    net$met.x[!is.na(m)] <- kegg.db$mets2collapse$to[na.omit(m)]
    m <- match(net$met.y, kegg.db$mets2collapse$from)
    net$met.y[!is.na(m)] <- kegg.db$mets2collapse$to[na.omit(m)]
    net <- unique(net)
    
    net <- net[net$rxn %in% rxns2keep,]
    
    net <- net[!(net$met.x %in% kegg.db$mets2mask) & !(net$met.y %in% kegg.db$mets2mask),]
    
    rownames(net) <- NULL
    
    met2name <- kegg.db$met2name
    rxn2name <- kegg.db$rxn2name
    met2name$name <- gsub(";.*$", "", met2name$name)
    rxn2name$name <- gsub(";.*$", "", rxn2name$name)
    
    kegg.network <- newEmptyObject()
    kegg.network$enz2gene <- enz2gene
    kegg.network$rxn2enz <- kegg.db$rxn2enz
    kegg.network$rxn2gene <- rxn2gene
    kegg.network$graph.raw <- net
    kegg.network$met.ids <- "KEGG"
    kegg.network$gene.ids <- "Entrez"
    kegg.network$met2name <- met2name
    kegg.network$rxn2name <- rxn2name
    kegg.network$organism <- organism
    return(kegg.network)
}
