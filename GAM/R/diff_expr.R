#' @import DESeq limma
NULL

#' @export
convertPvalBiomart <- function(pval, from, to, mart) {
    map <- getBM(attributes=c(from, to), filters=from, values=pval$ID, mart=mart)
    colnames(map) <- c("from", "to")
    map$to <- as.character(map$to)
    map <- na.omit(map)
    convert.pval(pval, map$from, map$to)
}

#' @export
convertPval <- function(pval, from, to) {
    map <- data.frame(from=from, to=to, stringsAsFactors=F)
    pval.ext <- merge(map, pval, by.x = "from", by.y = "ID")
    origin.field = if ("origin" %in% colnames(pval)) "origin" else "from"
    
    res <- data.frame(
        ID=pval.ext[,"to"], 
        pval=pval.ext$pval, 
        logFC=pval.ext$logFC, 
        origin=pval.ext[,origin.field],
        stringsAsFactors=F)    
    
    get_min_pval <- function(id) {
        subset <- res[res$ID == id,]
        
        # :ToDo: more intelligent aggregate function?
        subset[which.min(subset$pval),]        
    }
    
    res <- do.call("rbind", lapply(unique(res$ID), get_min_pval))
    res <- res[order(res$pval),]    
}

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

#' @export
fixInf <- function(dm) {    
    dm[dm == -Inf] <- min(dm[dm != -Inf]) - 1
    dm[dm == Inf] <- max(dm[dm != Inf]) + 1
    dm
}

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
        res = data.frame(ID=f.top$ID, pval=f.top$adj.P.Val, logFC=-f.top$logFC, stringsAsFactors=F)        
    }
    return(res)
}
