#' Convert IDs in differential expression data using BioMart database
#' When there is multiple entries with the same after conversion 
#' only one with minimal p-value is kept.
#' @param pval Table to covnert
#' @param from Mart attribute for original IDs
#' @param to Mart attribute for result IDs
#' @param mart Mart to use
#' @return Table with IDs converted
#' @export
convertPvalBiomart <- function(pval, from, to, mart) {
    map <- getBM(attributes=c(from, to), filters=from, values=pval$ID, mart=mart)
    colnames(map) <- c("from", "to")
    map$to <- as.character(map$to)
    map <- na.omit(map)
    convert.pval(pval, map$from, map$to)
}

#' @export
getIdsType <- function(ids, id.map) {
    res <- c()
    for (id.type in names(id.map)) {
        if (any(ids %in% id.map[, id.type])) {
            res <- c(res, id.type)
        }
    }
    res
}

#' Convert IDs in differential expression data
#' When there is multiple entries with the same after conversion 
#' only one with minimal p-value is kept.
#' @param pval Table to covnert
#' @param from Vector of IDs to convert from
#' @param to Vector of IDs to convert to 
#' @return Table with IDs converted
#' @importFrom plyr rename
#' @export
convertPval <- function(pval, from, to) {
    map <- data.frame(from=from, to=to, stringsAsFactors=F)
    pval.ext <- merge(map, pval, by.x = "from", by.y = "ID")
    
    if ("origin" %in% colnames(pval)) {
        pval.ext$from <- NULL
    } else {
        pval.ext <- rename(pval.ext, c("from"="origin"))
    }
    
    res <- rename(pval.ext, c("to"="ID"))
    
    get_min_pval <- function(id) {
        subset <- res[res$ID == id,]
        
        # :ToDo: more intelligent aggregate function?
        subset[which.min(subset$pval),]        
    }
    
    res <- merge(aggregate(pval ~ ID, res, min), res)
    res <- res[!duplicated(res$ID),]
    res <- res[order(res$pval),]    
    res
}

#' Normalize expression table
#' @param exprs Table with expressions
#' @param zero.rm If TRUE removes genes with zero expression in all samples
#' @param log2 It TRUE applies log2 transform. Zeroes are replaced with minimal non-zero element for sample
#' @param quantile If TRUE applies quantile normalization
#' @export
normalizeExpressions <- function(exprs, zero.rm=T, log2=T, quantile=T) {
    if (zero.rm) {
        # removing unexpressed genes
        keep <- apply(exprs, 1, max) > 0
        exprs <- exprs[keep,]
    }
    if (log2) {
        # adding pseudocount for zero expressions  
        min2s <- apply(exprs, 2, function(x) { min(x[x != 0]) })
        
        for (sample in colnames(exprs)) {
            t <- exprs[,sample]
            t[t == 0] <- min2s[sample];
            exprs[,sample] <- t
        }
        
        exprs <- log2(exprs)
    }
    
    if (quantile) {
        exprs2 <- as.data.frame(normalize.quantiles(as.matrix(exprs), copy=T))
        colnames(exprs2) <- colnames(exprs)
        rownames(exprs2) <- rownames(exprs)
        exprs <- exprs2
    }
    return(exprs)
}

#' Replace infinities with real numbers
#' This function replaces +Inf with maximal value + 1 and
#' -Inf with minimal value - 1
#' @param dm Vector of nummbers
#' @return Vector of numbers withoud infinities
#' @export
fixInf <- function(dm) {    
    dm[dm == -Inf] <- min(dm[dm != -Inf]) - 1
    dm[dm == Inf] <- max(dm[dm != Inf]) + 1
    dm
}

#' Make differential expression analysis using limma or DESeq
#' @param exprs Table with expressions
#' @param conditions.vector Vector of conditions for exprs columns
#' @param state1 First condition to compare
#' @param state2 Second condition to compare
#' @param top How many most expressed genes should be kept
#' @param log2 If TRUE apply log2 transformation (zeros replaced with minimal value in column)
#' @param quantile Apply quantile normalisation
#' @param use.deseq Use DESeq for analysis
#' @return Table with p-values for differential expression and log-fold changes
#' @import DESeq 
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @export
diffExpr <- function(exprs, conditions.vector, state1, state2, top=10000, log2=F, quantile=F, use.deseq=F) {
    expression_matrix<-as.matrix(exprs)
    
    
    classes_vector <- as.character(conditions.vector)
    
    
    unique_classes_vector<-unique(classes_vector)
    unique_classes_vector
    number_samples<-ncol(expression_matrix)
    number_classes<-length(unique_classes_vector)
    print(paste('Number of samples:',number_samples,'      Number of Cell Types:',number_classes))
    
    
    ####################
    
    if (use.deseq) {
        counts <- expression_matrix
                
        cds <- newCountDataSet(round(counts), classes_vector)
        cds <- estimateSizeFactors(cds)
        cds <- estimateDispersions(cds)
        res <- nbinomTest(cds, state1, state2)
        res <- res[res$baseMean > 5,]
        res <- res[order(res$baseMean, decreasing=T),]        
        res <- res[1:top,]
        res <- res[order(res$pval),]
        res <- na.omit(res)
        res <- res[, c("id", "pval", "log2FoldChange")]
        colnames(res) <- c("ID", "pval", "logFC")
    } else {
        log_expression_matrix<-normalizeExpressions(expression_matrix, zero.rm=T, log2=log2, quantile=quantile)
        group_names<-classes_vector
        
        names(group_names)<-group_names
        group_names<-factor(group_names)
        
        design<-model.matrix(~0+group_names)
        colnames(design)<-levels(group_names)
        
        colnames(log_expression_matrix)<-group_names
        
        fit<-lmFit(log_expression_matrix, design)
        j = which(unique_classes_vector == state1)
        if (length(j) != 1) {
            stop(c("invalied state '", state1, "'", sep=""))
        }
        jj = which(unique_classes_vector == state2)
        if (length(jj) != 1) {
            stop(c("invalied state '", state2, "'", sep=""))
        }
        
        contr.str<-paste(unique_classes_vector[j],unique_classes_vector[jj],sep="-")
        contr.mat<-makeContrasts(contrasts=contr.str, levels=group_names)
        
        fit2<-contrasts.fit(fit,contr.mat)
        fit2<-eBayes(fit2)
        
        f.top<-topTable(fit2, sort.by="B", number=top, resort.by="P")
        ids <- if ("ID" %in% names(f.top)) f.top$ID else rownames(f.top)
        res = data.frame(ID=ids, pval=f.top$adj.P.Val, logFC=-f.top$logFC, stringsAsFactors=F)        
    }
    return(res)
}
