# Normalize expression table
# @param exprs Table with expressions
# @param zero.rm If TRUE removes genes with zero expression in all samples
# @param log2 It TRUE applies log2 transform. Zeroes are replaced with minimal non-zero element for sample
# @param quantile If TRUE applies quantile normalization
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

diffExprM <- function(exprs, conditions.vector, states, top=Inf, log2=F, quantile=F, method=c("edgeR", "limma"), log2.min.expr=-Inf) {
    method <- match.arg(method)
    exprs <-as.matrix(exprs)
    
    
    conditions.vector <- as.character(conditions.vector)
    



    if (method == "edgeR") {
        if (!require(edgeR)) {
            stop("metdho=edgeR needs edgeR package to work")
        }
        exprs.normalized <- GAM:::normalizeExpressions(exprs, log2=T, quantile=T)
        y <- DGEList(counts=round(exprs), group=gene.conditions.vector)
        y <- calcNormFactors(y)
        y <- estimateCommonDisp(y)
        y <- estimateTagwiseDisp(y)
        design <- model.matrix(~group, data=y$samples)
        fit <- glmFit(y, design)
        lrt <- glmLRT(fit, coef=2:length(unique(y$samples$group)))
        res <- as.data.frame(topTags(lrt, n=Inf))
        res <- data.frame(ID=rownames(res), pval=res$FDR, stringsAsFactors=F)
    } else if (method == "limma") {
        if (!require(limma)) {
            stop("method=limma needs limma package to work")
        }
        
        exprs.normalized <- normalizeExpressions(exprs, zero.rm=T, log2=log2, quantile=quantile)
        
        names(conditions.vector) <- conditions.vector
        conditions.vector <- factor(conditions.vector)
        
        design <- model.matrix(~0 + conditions.vector)
        colnames(design) <- levels(conditions.vector)
        
        colnames(exprs.normalized) <- conditions.vector
        
        conditions.levels <- unique(conditions.vector)
        
        fit <- lmFit(exprs.normalized, design)
        
        contrasts <- c()
        for (i in seq_along(states)) {
            for (j in seq_len(length(states) - i) + i) {
                contrasts <- c(contrasts, paste0(states[i], "-", states[j]))
            }
        }
        contr.mat <- makeContrasts(contrasts=contrasts, levels=conditions.vector)
        
        fit2 <- contrasts.fit(fit, contr.mat)
        fit2 <- eBayes(fit2)
        
        pvals <- p.adjust(fit2$F.p.value, method="BH")
        
        res <- data.frame(ID=rownames(fit2), pval=pvals, stringsAsFactors=F)        
    }
    
    conditions.vector <- as.character(conditions.vector)
    conditions.logFCs <- getConditionLogFCs(
        exprs.normalized[, conditions.vector %in% states],
        conditions.vector[conditions.vector %in% states])
    res <- merge(res, conditions.logFCs)
    res <- res[res$max> log2.min.expr,]
    res <- res[order(res$max, decreasing=T),]        
    res <- head(res, n=top)
    res <- res[order(res$pval),]
    res <- na.omit(res)
    rownames(res) <- seq_len(nrow(res))
    return(res)
}

appendToNames <- function(x, suffix, ignore=c()) {
    require(plyr)
    rename(x, 
           sapply(
                colnames(x),
                function(n) if (n %in% ignore) n else paste0(n, suffix))) 
}

prependToNames <- function(x, prefix, ignore=c()) {
    require(plyr)
    rename(x, 
           sapply(
                names(x),
                function(n) if (n %in% ignore) n else paste0(prefix, n))) 
}


aggregateColumns <- function(x, columns, new.name, f) {
    x[[new.name]] <- apply(x[, columns], 1, f)
    x
}

combineDE <- function(de.all) {
    de.combined <- Map(
        function(n) appendToNames(
            de.all[[n]],
            paste0(".", n), ignore="ID"),
        names(de.all))
    
    de.combined <- Reduce(function(...) merge(..., all=T, by="ID") , de.combined)
    
    aggregateColumns(
        de.combined, 
        grep("^pval", colnames(de.combined)),
        "pval",
        function(x) min(x, na.rm=T))
}

getConditionLogFCs <- function(logexprs, conditions.vector) {
    condition.means <- data.frame(sapply(sort(unique(conditions.vector)), function(condition) {
        print(condition)
        condition.logexprs <- logexprs[, conditions.vector == condition, drop=F]
        apply(condition.logexprs, 1, mean)
    }))
    
    base.means <- apply(condition.means, 1, mean)
    maxs <- apply(condition.means, 1, max)
    
    res <- cbind(baseMean=base.means, max=maxs, prependToNames(condition.means - base.means, "log2FC."))
    res$ID <- names(base.means) 
    res
}
