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
#'     gene.de.M1.M2.entrez <- convertPvalBiomart(gene.de.M1.M2, "refseq_mrna", "entrezgene", mart)
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
#' data(met.id.map)
#' data(examplesGAM)
#' met.de.M1.M2.kegg <- convertPval(met.de.M1.M2, met.id.map$HMDB, met.id.map$KEGG)
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

# Normalize expression table
# @param exprs Table with expressions
# @param zero.rm If TRUE removes genes with zero expression in all samples
# @param log2 It TRUE applies log2 transform. Zeroes are replaced with minimal non-zero element for sample
# @param quantile If TRUE applies quantile normalization
normalizeExpressions <- function(exprs, zero.rm=TRUE, log2=TRUE, quantile=TRUE) {
    if (zero.rm) {
        # removing unexpressed genes
        keep <- apply(exprs, 1, function(x) max(x, na.rm = T)) > 0
        exprs <- exprs[keep,]
    }
    if (log2) {
        # adding pseudocount for zero expressions  
        min2s <- apply(exprs, 2, function(x) { min(x[x != 0], na.rm=T) })
        
        for (sample in colnames(exprs)) {
            t <- exprs[,sample]
            t[t == 0] <- min2s[sample];
            exprs[,sample] <- t
        }
        
        exprs <- log2(exprs)
    }
    
    if (quantile) {
        if (!require(limma)) {
            stop("quantile normalization needs limma package to work")
        }
        exprs2 <- as.data.frame(normalizeBetweenArrays(as.matrix(exprs), method="quantile"))
        colnames(exprs2) <- colnames(exprs)
        rownames(exprs2) <- rownames(exprs)
        exprs <- exprs2
    }
    return(exprs)
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

#' Make differential expression analysis using limma or DESeq
#' @param exprs Table or matrix with expressions
#' @param conditions.vector Vector of conditions for exprs columns
#' @param state1 First condition to compare
#' @param state2 Second condition to compare
#' @param top How many most expressed genes should be kept
#' @param log2 If TRUE apply log2 transformation (zeros replaced with minimal value in column)
#' @param quantile Apply quantile normalisation
#' @param use.deseq Use DESeq for analysis
#' @param min.expr Minimal mean expression for a feature to be kept
#' @return Table with p-values for differential expression and log-fold changes
#' @examples
#' \dontrun{
#' if (require(mouseMacrophages) && require(DESeq)) {
#'     data(mmpData)
#'     gene.de.M1.M2 <- diffExpr(
#'         exprs=exprs(mmpGeneSet), conditions.vector=pData(mmpGeneSet)$condition,
#'         state1="MandLPSandIFNg", state2="MandIL4", 
#'         use.deseq=TRUE, top=Inf, min.expr=5)
#' }
#' }
#' @export
diffExpr <- function(exprs, 
                     conditions.vector,
                     state1, state2,
                     top=Inf,
                     log2=FALSE, quantile=FALSE,
                     use.deseq=FALSE,
                     min.expr=0) {
    exprs <-as.matrix(exprs)
    
    
    conditions.vector <- as.character(conditions.vector)
    
    if (use.deseq) {
        if (!require(DESeq)) {
            stop("use.deseq=TRUE needs DESeq package to work")
        }
        cds <- newCountDataSet(round(exprs), conditions.vector)
        cds <- estimateSizeFactors(cds)
        cds <- estimateDispersions(cds)
        res <- nbinomTest(cds, state1, state2)
        res <- res[, c("id", "pval", "log2FoldChange", "baseMean")]
        colnames(res) <- c("ID", "pval", "log2FC", "baseMean")
    } else {
        if (!require(limma)) {
            stop("use.deseq=FALSE needs limma package to work")
        }
        
        exprs.normalized <- normalizeExpressions(exprs, zero.rm=TRUE, log2=log2, quantile=quantile)
        
        names(conditions.vector) <- conditions.vector
        conditions.vector <- factor(conditions.vector)
        
        design <- model.matrix(~0 + conditions.vector)
        colnames(design) <- levels(conditions.vector)
        
        colnames(exprs.normalized) <- conditions.vector
        
        conditions.levels <- unique(conditions.vector)
        
        fit <- lmFit(exprs.normalized, design)
        j <- which(conditions.levels == state1)
        if (length(j) != 1) {
            stop(c("invalied state '", state1, "'", sep=""))
        }
        jj <- which(conditions.levels == state2)
        if (length(jj) != 1) {
            stop(c("invalied state '", state2, "'", sep=""))
        }
        
        contr.str <- paste(conditions.levels[j], conditions.levels[jj], sep="-")
        contr.mat <- makeContrasts(contrasts=contr.str, levels=conditions.vector)
        
        fit2 <- contrasts.fit(fit,contr.mat)
        fit2 <- eBayes(fit2)
        
        f.top <- topTable(fit2, number=Inf)
        ids <- if ("ID" %in% names(f.top)) f.top$ID else rownames(f.top)
        res <- data.frame(ID=ids, pval=f.top$adj.P.Val, log2FC=-f.top$logFC, baseMean=f.top$AveExpr, stringsAsFactors=FALSE)        
    }
    res <- res[res$baseMean > min.expr,]
    res <- res[order(res$baseMean, decreasing=TRUE),]        
    res <- head(res, n=top)
    res <- res[order(res$pval),]
    res <- na.omit(res)
    rownames(res) <- seq_len(nrow(res))
    return(res)
}
