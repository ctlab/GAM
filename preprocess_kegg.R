enz2gene <- read.table("./networks//kegg/enz2gene.tsv", header=T, colClasses="character")
rxn2enz <- read.table("./networks//kegg/rxn2enz.tsv", header=T, colClasses="character")
net <- read.table("./networks//kegg/net.sif", header=T, colClasses="character")
rxn2name <- read.table("./networks//kegg/rxn2name.tsv", header=T, colClasses="character")
met2name <- read.table("./networks//kegg/met2name.tsv", header=T, colClasses="character")


#rxns2keep <- read.table("./networks//kegg/rxns2keep.lst", colClasses="character")[,1]

mets2mask <- read.table("./networks//kegg/mets2mask.lst", colClasses="character")[,1]


rxn2gene = merge(rxn2enz, enz2gene)[, c("rxn", "gene")]
rxns2keep <- unique(rxn2gene$rxn)
write.table(rxns2keep, "./networks//kegg/rxns2keep.lst", row.names=F, col.names=F, quote=F)




get_rxn2met <- function(net) {
    rxn2met.x <- data.frame(rxn=net$rxn, met=net$met.x, stringsAsFactors=F)
    rxn2met.y <- data.frame(rxn=net$rxn, met=net$met.y, stringsAsFactors=F)
    return(unique(rbind(rxn2met.x, rxn2met.y)))
}

get_met.degree <- function(net) {
    rxn2met <- get_rxn2met(net)
    met.degree <- table(rxn2met$met)
    met.degree <- met.degree[order(met.degree, decreasing=T)]
    met.degree <- as.data.frame(met.degree)
    met.degree$name <- met2name$name[match(rownames(met.degree), met2name$met)]
    return(met.degree)
}

met.degree <- get_met.degree(net)
write.csv(met.degree, "met.degree.1.csv")

net <- net[net$rxn %in% rxns2keep,]

met.degree <- get_met.degree(net)
write.csv(met.degree, "met.degree.2.csv")

net <- net[!(net$met.x %in% mets2mask) & !(net$met.y %in% mets2mask),]
met.degree <- get_met.degree(net)
write.csv(met.degree, "met.degree.3.csv")

rxn2met <- get_rxn2met(net)
met2rxn = rxn2met[,c("met", "rxn")]

rxn2rxn = unique(merge(rxn2met, met2rxn, by="met")[,c("rxn.x", "rxn.y")])

write.table(rxn2rxn, "./networks//kegg/rxn2rxn.masked.tsv", row.names=F, quote=F, col.names=T)



net.sq.1 = cbind(u=net$met.x, e="m2m", v=net$met.y)
net.sq.2 = cbind(u=rxn2met$rxn, e="r2m", v=rxn2met$met)
net.sq.3 = cbind(u=rxn2rxn$rxn.x, e="r2r", v=rxn2rxn$rxn.y)
net.sq <- rbind(net.sq.1, net.sq.2, net.sq.3)

write.table(net.sq, "./networks//kegg/net.sq.sif", row.names=F, quote=F, col.names=F)

long_names.met <- met2name$name
names(long_names.met) <- met2name$met
long_names.rxn <- rxn2name$name
names(long_names.rxn) <- rxn2name$rxn

long_names <- data.frame(longName = na.omit(c(long_names.met, long_names.rxn)), stringsAsFactors=F)
write.table(long_names, "./networks//kegg/net.sq_LongName.NA", row.names=T, quote=T, col.names=T, sep=" = ")
# remove quotes from first line!

short_names <- long_names
colnames(short_names) <- "shortName"
short_names$shortName <- gsub(";.*$", "", short_names$shortName)
write.table(short_names, "./networks//kegg/net.sq_ShortName.NA", row.names=T, quote=T, col.names=T, sep=" = ")
# remove quotes from first line!
