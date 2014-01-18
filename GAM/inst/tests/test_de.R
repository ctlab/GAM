context("Differential expression")

test_that("diffExpr works", { # :ToDo:
})

test_that("convertPval works", { # :ToDo:
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

test_that("convertPvalBiomart works", { # :ToDo:
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

