context("Differential expression")


genExpressions <- function(conditions.vector, distribution.parameters) {
    res <- rep(0, length(conditions.vector))
    for (name in names(distribution.parameters)) {
        mean <- distribution.parameters[[name]][1]
        sd <- distribution.parameters[[name]][2]
        res[conditions.vector == name] <- 
            rnorm(sum(conditions.vector == name), mean, sd)
    }
    return(res)
}

test_that("diffExpr works with limma", {
    if (require(mouseMacrophages) && require(limma)) {
        data(mmpData)
        met.exprs <- exprs(mmpMetSet)
        set.seed(42)
        n <- 100
        met.exprs <- met.exprs[sample(rownames(met.exprs), n), ]
        state1 <- "MandLPSandIFNg"
        state2 <- "MandIL4"
        met.exprs <- rbind(
            met.exprs,
            HMDBx1=genExpressions(pData(mmpMetSet)$condition, 
                       list(
                           MandIL4=c(15, 0.2), 
                           MandLPSandIFNg=c(13, 0.2))),
            HMDBx2=genExpressions(pData(mmpMetSet)$condition, 
                       list(
                           MandIL4=c(8, 0.3), 
                           MandLPSandIFNg=c(8.2, 0.3))))
        met.de <- diffExpr(
            exprs=met.exprs, conditions.vector=pData(mmpMetSet)$condition,
            state1=state1, state2=state2,
            use.deseq=FALSE, top=Inf)
        expect_true(met.de[met.de$ID == "HMDBx1",]$pval < 1e-5)
        expect_true(met.de[met.de$ID == "HMDBx2",]$pval > 1e-5)
    }
})


test_that("diffExpr works with DESeq", {
    if (require(mouseMacrophages) && require(DESeq)) {
        data(mmpData)
        gene.exprs <- exprs(mmpGeneSet)
        set.seed(42)
        n <- 1000
        gene.exprs <- gene.exprs[sample(rownames(gene.exprs), n), ]
        state1 <- "MandLPSandIFNg"
        state2 <- "MandIL4"
        gene.exprs <- rbind(
            gene.exprs,
            NMx_1=genExpressions(pData(mmpGeneSet)$condition, 
                       list(
                           MandIL4=c(1000, 200), 
                           MandLPSandIFNg=c(100, 50))),
            
            NMx_2=genExpressions(pData(mmpGeneSet)$condition, 
                       list(
                           MandIL4=c(200, 100), 
                           MandLPSandIFNg=c(200, 100))))
        gene.de <- diffExpr(
            exprs=gene.exprs, conditions.vector=pData(mmpGeneSet)$condition,
            state1=state1, state2=state2,
            use.deseq=TRUE, top=Inf, min.expr=5)
        expect_true(gene.de[gene.de$ID == "NMx_1",]$pval < 1e-5)
        expect_true(gene.de[gene.de$ID == "NMx_2",]$pval > 1e-5)
    }
})

test_that("convertPval works", {
    data(met.id.map)
    met.de.hmdb <- data.frame(
        ID=c("HMDB14289", "HMDB01919", "HMDB02092"),
        pval=c(     1e-5,       1e-20,       1e-10),
        log2FC=c(    0.5,        -4.2,         2.7), 
        stringsAsFactors=FALSE)
    met.de.kegg <- convertPval(met.de.hmdb, met.id.map$HMDB, met.id.map$KEGG)
    expect_equal(met.de.kegg$pval[met.de.kegg$ID == "C00025"], 1e-20)
    expect_equal(met.de.kegg$origin[met.de.kegg$ID == "C00490"], "HMDB02092")
})

test_that("convertPvalBiomart works", {
    if (require(biomaRt) && require(RCurl)) {
        # For biomaRt to work with Squid proxy, see http://comments.gmane.org/gmane.science.biology.informatics.conductor/39218
        options(RCurlOptions = list(http.version=HTTP_VERSION_1_0))
    
        gene.de.refseq<- data.frame(
            ID=c("NM_010927", "NM_011198"),
            pval=c(     1e-6,       1e-11),
            log2FC=c(    0.9,         1.6), 
            stringsAsFactors=FALSE)
        
        mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")    
        gene.de.entrez <- convertPvalBiomart(
            gene.de.refseq, 
            "refseq_mrna", "entrezgene", mart)
        
        expect_equal(gene.de.entrez$pval[gene.de.entrez$ID == "18126"], 1e-6)
        expect_equal(gene.de.entrez$origin[gene.de.entrez$ID == "19225"],
                     "NM_011198")
        
    }
})

