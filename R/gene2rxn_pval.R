#!/usr/bin/Rscript
# converts genes to reactions

gene2rxn_pval.main <- function() {

    #rxn2genes_map_file = "rxn2genes.csv"
    #gene_pvals_file = "./data/metabolome/MandIL4-MOandLPSandIFNg.pval"
    
    library(optparse)
    
    option_list <- list(
        #    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
        #                help="Print extra output [default]"),
        #    make_option(c("-q", "--quietly"), action="store_false",
        #                dest="verbose", help="Print little output"),
        make_option("--rxn2genes",
                    dest="rxn2genes_map_file",
                    help="File with mapping reaction to genes it depends on",
                    metavar="tsv file"),
        make_option("--gene-pvals",
                    dest="gene_pvals_file",
                    help="File with gene differential expression p-values",
                    metavar="file"),
        make_option(c("-o", "--output-file"),
                    dest="output_file",
                    help="File to output reaction p-values",
                    metavar="file")
    )
    
    opt <- parse_args(OptionParser(option_list=option_list))
    
    rxn2genes <- read.csv(opt$rxn2genes_map_file, header= TRUE, check.names=F, sep="\t")
    
    #pval<-read.table(gene_pvals_file, row.names=1)
    # somehow pval table may have duplicate entries
    pval<-read.table(opt$gene_pvals_file, header=TRUE)
    pval<-pval[!duplicated(pval[,1]),]
    
    # convert refseq to entrez
    if (T) {
        gene_id_type="refseq_mrna"
        #gene_id_type="mgi_symbol"
        library(biomaRt)
        mart = useMart("ensembl", dataset="mmusculus_gene_ensembl")
        map <- getBM(attributes=c("entrezgene", gene_id_type), filters=gene_id_type, values=pval$ID, mart=mart)
        head(map)
        head(pval)
        pval_ext <- merge(map, pval, by.x=gene_id_type, by.y="ID")
        head(pval_ext)
        
        pval_ext <- na.omit(pval_ext)
        pval <- pval_ext
        rm(pval_ext)
        pval$ID <- pval$entrezgene
    }
    
    rxn_pvals_ext <- merge(rxn2genes, pval, by.x = "gene", by.y = "ID")
    rxn_pvals_ext <- as.data.frame(cbind(
        ID=as.character(rxn_pvals_ext$rxn),
        pval=rxn_pvals_ext$pval,
        logFC=rxn_pvals_ext$logFC))
    
    get_rxn_pval <- function(rxn) {
        subset <- rxn_pvals_ext[rxn_pvals_ext$ID == rxn,]
        
        # :ToDo: more intelligent aggregate function?
        res <- subset[which.min(subset$pval),]
        return(res)
    }
    
    rxn_pvals <- do.call("rbind", lapply(levels(rxn_pvals_ext$ID), get_rxn_pval))
    rxn_pvals <- rxn_pvals[order(rxn_pvals$pval),]
    
    write.table(rxn_pvals, opt$output_file, quote=F, col.names=T, row.names=F, sep="\t")
}