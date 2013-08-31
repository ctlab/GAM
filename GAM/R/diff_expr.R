#!/usr/bin/env Rscript
library(limma)


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

diff.expr <- function(exprs, conditions.vector, state1, state2, top=10000, log2=F, quantile=F, use.deseq=F) {
    expression_matrix<-as.matrix(exprs)
    
    
    classes_vector <- as.character(conditions.vector)
    
    
    unique_classes_vector<-unique(classes_vector)
    unique_classes_vector
    number_samples<-ncol(expression_matrix)
    number_classes<-length(unique_classes_vector)
    print(paste('Number of samples:',number_samples,'      Number of Cell Types:',number_classes))
    
    
    ####################
    
    if (opt$deseq) {
        library(DESeq)
        counts <- expression_matrix
                
        cds <- newCountDataSet(round(counts), classes_vector)
        cds <- estimateSizeFactors(cds)
        cds <- estimateDispersions(cds)
        res <- nbinomTest(cds, opt$state1, opt$state2)
        res <- res[order(res$baseMean, decreasing=T)[1:top],]
        res <- res[order(res$pval),]
        res <- na.omit(res)
        res <- res[, c("id", "pval", "log2FoldChange")]
        colnames(res) <- c("ID", "pval", "logFC")
    } else {
        log_expression_matrix<-normalize_expressions(expression_matrix, zero.rm=T, log2=log2, quantile=quantile)
        group_names<-classes_vector
        
        names(group_names)<-group_names
        group_names<-factor(group_names)
        
        design<-model.matrix(~0+group_names)
        colnames(design)<-levels(group_names)
        
        colnames(log_expression_matrix)<-group_names
        
        fit<-lmFit(log_expression_matrix, design)
        j = which(unique_classes_vector == state1)
        if (length(j) != 1) {
            stop(c("invalied state '", opt$state1, "'", sep=""))
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
        res = data.frame(ID=f.top$ID, pval=f.top$adj.P.Val, logFC=-f.top$logFC)        
    }
    return(res)
}



main <-function() {
    library(optparse)
    
    option_list <- list(
        #    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
        #                help="Print extra output [default]"),
        #    make_option(c("-q", "--quietly"), action="store_false",
        #                dest="verbose", help="Print little output"),
        make_option(c("-e", "--expressions"),
                    dest="expressions_file",
                    help="File to read expressions from",
                    metavar="csv file"),
        make_option(c("-c", "--conditions"),
                    dest="conditions_file",
                    help="File to read conditions from",
                    metavar="csv file"),
        make_option(c("--deseq"),
                    dest="deseq",
                    action="store_true",
                    default=FALSE,
                    help="Use DESeq for differential expression (with non-normalized data)"),
        make_option(c("--log2"),
                    dest="log2",
                    action="store_true",
                    default=FALSE,
                    help="Use log2 normalization"),
        make_option(c("--quantile"),
                    dest="quantile",
                    action="store_true",
                    default=FALSE,
                    help="Use quantile normalization"),
        make_option(c("-s1", "--state-1"),
                    dest="state1",
                    help="First state for differential expression",
                    metavar="state"),
        make_option(c("-s2", "--state-2"),
                    dest="state2",
                    help="Second state for differential expression",
                    metavar="state"),
        make_option(c("-o", "--output-file"),
                    dest="output_file",
                    help="Output file",
                    metavar="file")
    )
    
    opt <- newEmptyObject()
    opt$expressions_file <- "./data_new_good/Gene.exp.tsv"
    opt$conditions_file <- "./data_new_good/Gene.conditions.csv"
    opt$log2=T
    opt$quantile=T
    opt$deseq=F
    opt$state1="MandLPSandIFNg"
    opt$state2="MandIL4"
    
    opt <- parse_args(OptionParser(option_list=option_list))
    exprs.sep=","
    
    if (grepl("tsv$", opt$expressions_file)) {
        exprs.sep <- "\t"
    }
    
    exprs<-read.csv(file=opt$expressions_file, head=TRUE, row.names=1, sep=exprs.sep)
    conditions <- read.csv(file=opt$conditions_file, head=TRUE)
    classes_vector <- conditions$condition[match(colnames(exprs), conditions$sample)]
    
    to_write <- diff.expr(exprs=exprs, conditions.vector=classes_vector,
                          state1=opt$state1, state2=opt$state2,
                          log2=opt$log2, quantile=opt$quantile,
                          use.deseq=opt$deseq)
    
    write.table(to_write, file=opt$output_file, sep="\t", row.names=F, quote=F)
    
}
