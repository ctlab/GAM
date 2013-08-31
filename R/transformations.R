MyProduct<-function(expM,expN){
	mypr<-expM %*% expN
	mypr<-mypr-diag(diag(mypr))#*IdentityM
	return(mypr)
}

sif2matrix<-function(filename){
	sif_data<-read.table(filename, header=FALSE, colClasses=c("character", "factor", "character"), comment.char="#")
	nodes<-unique(c(sif_data[,1],sif_data[,3]))
	dimen<-length(nodes)
	expM<-mat.or.vec(dimen,dimen)
	rownames(expM)<-nodes
	colnames(expM)<-nodes

    I <- match(sif_data$V1, nodes)
    J <- match(sif_data$V3, nodes)

    expM[ cbind(I, J) ] <- 1
    expM[ cbind(J, I) ] <- 1

	return(expM)
}

#subset_names<-c("Gene_2","Gene_4","Gene_5","Gene_6")

matrix_subset<-function(expM,subset_names){
	sub_expM<-subset(expM,rownames(expM) %in% subset_names,colnames(expM) %in% subset_names)
	return(sub_expM)
}

rel_env = new.env()
rel_env[["1"]] = "act"
rel_env[["-1"]] = "inh"

matrix2frame <- function(expM) {
    nonzeros<-which(expM!=0,arr.ind=TRUE)
    
    rn <- rownames(expM)[nonzeros[,1]]
    cn <- colnames(expM)[nonzeros[,2]]
    rel <- unlist(mget(as.character(expM[nonzeros]), rel_env), use.names=F)
    
    res = data.frame(rn, rel, cn)
}

#use matrix with rownames and colnames
#outputfile="test.sif"
matrix2sif<-function(expM,outputfile){
		res <- matrix2frame(expM)
        write.table(res, outputfile, col.names=F, row.names=F, quote=F, sep="\t")
}

matrix2tab<-function(expM,outputfile){
        res <- matrix2frame(expM)
        write.table(res[,c("rn", "cn")], outputfile, col.names=F, row.names=F, quote=F, sep="\t")
}


