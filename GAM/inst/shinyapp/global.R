library(GAM)

data(kegg.mouse.network)
data(kegg.human.network)
networks <- list(
    "Mouse musculus"=kegg.mouse.network,
    "Homo sapiens"=kegg.human.network)

gene.ids <- c("RefSeq mRNA accession"="RefSeq", 
                 "Entrez gene"="Entrez",
                 "Gene symbol"="name")
met.ids <- c("HMDB"="HMDB", "KEGG"="KEGG")