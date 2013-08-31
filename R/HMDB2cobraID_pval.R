#!/usr/bin/env Rscript

# convert genes to reactions

### make matrix that would map HMDB IDs one-to-many cobraIDs based on the output dictionary table (cobraIDs vertically and HMDBs horisontally)


#map<-mat.or.vec(dim(expM)[1],dim(expM)[1])

library(optparse)
option_list <- list(
#    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
#                help="Print extra output [default]"),
#    make_option(c("-q", "--quietly"), action="store_false",
#                dest="verbose", help="Print little output"),
    make_option(c("-i", "--pvals"),
                dest="pvals_file",
                help="File to read p-values from",
                metavar="file"),
    make_option(c("-o", "--output-file"),
                dest="output_file",
                help="Output file",
                metavar="file"),
    make_option("--hmdb2cobra-dict",
                dest="hmdb2cobra_dict",
                help="Matrixmap between rxns and genes that is outputed by Cobra",
                metavar="csv file")
)

opt <- parse_args(OptionParser(option_list=option_list))

outputfile=opt$output_file
outputfile_noncomp=paste("CobraID","noncomp",basename(opt$pvals_file),sep="_")



#read the hmdb to cobraID
map<-read.table(opt$hmdb2cobra_dict, header=T, colClasses="character")
pval<-read.table(opt$pvals_file, header=T)
rownames(pval) <- pval$ID
#pval<-pvals[!duplicated(pvals[,1]),]  # making sure there is no
#rownames(pval)<-pval[,1]

new_pval <- pval[map$HMDB,]
new_pval$ID <- map$Cobra
new_pval <- na.omit(new_pval)

write.table(new_pval, outputfile, quote=FALSE, row.names=F,col.names=T,sep="\t")
