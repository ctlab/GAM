#' Convert IDs in differential expression data using BioMart database
#' 
#' Converts IDs in differential expression data using BioMart database.
#' When there is multiple entries with the same after conversion 
#' only one with minimal p-value is kept. Expects \code{ID} and \code{pval}
#' columns.
#' 
#' @param pval Table to covnert
#' @param from Mart attribute for original IDs
#' @param to Mart attribute for result IDs
#' @param mart Mart to use
#' @return Table with IDs converted
#' @examples
#' \dontrun{
#' data(examplesGAM)
#' if (require(biomarRt)) {
#'     mart <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
#'     gene.de.M0.M1.entrez <- convertPvalBiomart(gene.de.M0.M1, "refseq_mrna", "entrezgene", mart)
#' }
#' }
#' @export
convertPvalBiomart <- function(pval, from, to, mart) {
    if (!require(biomaRt)) {
        stop("convertPvalBiomart needs biomaRt package to work")
    }
    map <- getBM(attributes=c(from, to), filters=from, values=pval$ID, mart=mart)
    colnames(map) <- c("from", "to")
    map$to <- as.character(map$to)
    map <- na.omit(map)
    convertPval(pval, map$from, map$to)
}

#' Get possible ID type  of IDs from a vector
#' 
#' Checks to which coulmn of id.map the biggest number 
#' of IDs is matched.
#' 
#' @param ids Vector of IDs
#' @param id.map Map between different type of IDs, better to be data.table
#' @return Type of IDs in a vector or NULL if no suitable type was found
#' @examples
#' library(GAM.db)
#' data(met.id.map)
#' id.type <- getIdType(c("HMDB00001", "HMDB02092"), met.id.map)
#' @export
getIdType <- function(ids, id.map) {
    if (!is(id.map, "data.table")) {
        warning("Converting id.map to data.table")
        id.map <- as.data.table(id.map)
    }
    counts <- sapply(names(id.map), function(id.type) {
        sum(ids %in% id.map[, get(id.type)])    
    })
    
    if (max(counts) == 0) {
        return(NULL)
    }
    
    names(counts[which.max(counts)])
}

#' Convert IDs in differential expression data
#' 
#' Converts IDs in differential expression data with provided mapping.
#' When there is multiple entries with the same after conversion 
#' only one with minimal p-value is kept. Expects \code{ID} and \code{pval}
#' columns.
#' @param pval Table to covnert
#' @param from Vector of IDs to convert from
#' @param to Vector of IDs to convert to 
#' @return Table with IDs converted
#' @importFrom plyr rename
#' @examples
#' library(GAM.db)
#' data(met.id.map)
#' data(examplesGAM)
#' met.de.M0.M1.kegg <- convertPval(met.de.M0.M1, met.id.map$HMDB, met.id.map$KEGG)
#' @export
convertPval <- function(pval, from, to) {
    map <- data.frame(ID=from, to=to, stringsAsFactors=FALSE)
    pval.ext <- merge(map, pval, by="ID")
    if ("origin" %in% colnames(pval)) {
        pval.ext$ID<- NULL
    } else {
        pval.ext <- rename(pval.ext, c("ID"="origin"))
    }
    
    res <- rename(pval.ext, c("to"="ID"))
    
    res <- merge(aggregate(pval ~ ID, res, min), res)
    res <- res[!duplicated(res$ID),]
    res
}

# Replace infinities with real numbers

# This function replaces +Inf with maximal value + 1 and
# -Inf with minimal value - 1
# @param dm Vector of nummbers
# @return Vector of numbers withoud infinities
fixInf <- function(dm) {    
    
    if (any(na.omit(dm == -Inf))) {
        dm[dm == -Inf] <- min(dm[dm != -Inf]) - 1
    }
    if (any(na.omit(dm == Inf))) {
        dm[dm == Inf] <- max(dm[dm != Inf]) + 1
    }
    dm
}
