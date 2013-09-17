networkFromSif <- function(file.basename) {
    print("Loading network...")
    network <- newEmptyObject()
    class(network) <- "network"
    network$edges <- read.table(paste(file.basename, "sif", sep="."), header=FALSE, colClasses=c("character", "factor", "character"), comment.char="#", col.names=c("u", "e", "v"))
    
    for (f in list.files(dirname(file.basename), paste(basename(file.basename), "_\\w+.NA", sep=""), full.names=T)) {
        meta_type.cur <- sub(paste("^", file.basename, "_(\\w+).NA", sep=""), "\\1", f)        
        print(paste("Loading attribute", meta_type.cur))
        meta_temp <- unique(read.table(f, skip=1, comment.char="#", colClasses="character"))
        meta.cur <- data.frame(meta_temp[,1], meta_temp[,3], stringsAsFactors=F)
        colnames(meta.cur) <- c("v", meta_type.cur)
        if (is.null(network$meta)) {
            network$meta <- meta.cur
        } else {
            network$meta <- merge(network$meta, meta.cur, by="v", all=T)
        }        
    }
    
    return(network)
}

networkToSif <- function(network, file.basename) {
    write.table(network$edges[,c("u", "e", "v")], paste(file.basename, "sif", sep="."),
                col.names=F, row.names=F, quote=F, sep="\t")
    for (m in colnames(network$meta)) {
        if (m != "v") {
            x <- network$meta[,m, drop=F]
            rownames(x) <- network$meta$v
            x <- na.omit(x)
            f <- paste(file.basename, "_", colnames(x), ".NA", sep="")
            write.table(x, f, row.names=T, quote=T, col.names=T, sep=" = ")
        }
    }
}

networkToTab <- function(edges, file.name) {
    write.table(edges[,c("u", "v")], file.name, col.names=F, row.names=F, quote=F, sep="\t")
}