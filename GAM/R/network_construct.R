organism.id.map <- rbind(
    c("MMU",  "org.Mm.eg"),
    c("HSA",  "org.Hs.eg")
    )
colnames(organism.id.map) <- 
    c("KEGG", "Annotation")

organism.id.map <- data.frame(organism.id.map, stringsAsFactors=FALSE)

makeGeneIdMap <- function(org.annotation, ids = c("REFSEQ", "SYMBOL")) {
    annotation.package <- paste0(org.annotation, ".db")
    if (!do.call("require", list(annotation.package))) {
        stop(paste("Can't load package", annotation.package, "to make gene ID map"))
    }
#     annotation.namespace <- getNamespace(annotation.package)
    bimaps <- lapply(ids, 
                     function(id) AnnotationDbi::as.data.frame(get(paste0(org.annotation, id))))
    
    res <- bimaps[[1]]
    for (i in seq_len(length(ids) - 1)+1) {
        res <- merge(res, bimaps[[i]], all=TRUE)
    }
    res
}

#' Make KEGG reaction network for a specified organism
#' @param kegg.db object with KEGG mappings
#' @param organism.kegg Organism KEGG ID (only MMU or HSA by default)
#' @param organism.annotation Genome wide annotation package to use (without ".db" suffix)
#' @return Network for the organism usable for an analysis
#' data(kegg.db)
#' if (require(org.Mm.eg)) {
#'     kegg.mouse.network <- makeKeggNetwork(kegg.db, "MMU")
#' }
#' @export
makeKeggNetwork <- function(kegg.db, organism.kegg, organism.annotation=organism.id.map$Annotation[organism.id.map$KEGG == organism.kegg]) {
    enz2gene <- kegg.db$enz2gene[kegg.db$enz2gene$organism == organism.kegg, c("enz", "gene")]
    rxn2enz <- kegg.db$rxn2enz
    rxn2gene <- merge(rxn2enz, enz2gene)[, c("rxn", "gene")]
    
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
    
    gene.id.map <- makeGeneIdMap(organism.annotation, 
                                 c("REFSEQ", "SYMBOL"))
    colnames(gene.id.map) <- c("Entrez", "RefSeq", "Symbol")
    
    gene.id.map <- gene.id.map[gene.id.map$Entrez %in% enz2gene$gene, ]
    gene.id.map <- as.data.table(gene.id.map)
    
    
    mets <- unique(c(net$met.x, net$met.y))
    rxns <- unique(net$rxn)
    
    met2name <- met2name[met2name$met %in% mets, ]
    rxn2name <- rxn2name[rxn2name$rxn %in% rxns, ]
    
    
    kegg.network <- newEmptyObject()
#     kegg.network$enz2gene <- enz2gene
#     kegg.network$rxn2enz <- rxn2enz
    kegg.network$rxn2gene <- rxn2gene
    kegg.network$graph.raw <- net
    kegg.network$met.ids <- "KEGG"
    kegg.network$gene.ids <- "Entrez"
    kegg.network$met2name <- met2name
    kegg.network$rxn2name <- rxn2name
    kegg.network$organism <- organism.kegg
    kegg.network$gene.id.map <- gene.id.map
    return(kegg.network)
}
