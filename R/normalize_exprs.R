
normalize_expressions <- function(exprs, zero.rm=T, log2=T, quantile=T) {
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
        library(preprocessCore)
        exprs2 <- as.data.frame(normalize.quantiles(as.matrix(exprs), copy=T))
        colnames(exprs2) <- colnames(exprs)
        rownames(exprs2) <- rownames(exprs)
        exprs <- exprs2
    }
    return(exprs)
}