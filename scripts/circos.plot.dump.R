circos.dump.line.plot <- function(chr, pos, data,file){
  stopifnot(length(chr)==length(pos) || length(pos)==length(data))
  output <- data.frame(chr=paste('mm',tolower(chr),sep=""),left=pos,right=pos,data=as.vector(data))
  write.table(output,row.names=F,col.names=FALSE,sep='\t',file=file,qmethod='escape',quote=FALSE)
}


name <- "129xPWK"
parental.color <- "pink"
marental.color <- "red"

name <- "PWKx129"
parental.color <- "red"
marental.color <- "pink"


t<-read.table(paste("result/cuffsequence/",name,"/",name,".info",sep=""),stringsAsFactors=F)
N <- length(t[,1])
cuffcmp<-read.table(paste("result/cuffsequence/",name,"/cuffcmp.tracking",sep=""),stringsAsFactors=F)
meta <- do.call(rbind,strsplit(cuffcmp[,5],"\\|"))
gene.index <- do.call(rbind,strsplit(cuffcmp[,3],"\\|"))
gene.index <- gene.index[,1]
transcript.index <- do.call(rbind,strsplit(meta[,1],":"))
transcript.index <- transcript.index[,2]
rownames(cuffcmp)=meta[,2]


total <- list()
s129 <- list()
spwk <- list()
pos <- list()
chr <- list()
length <- list()
for (i in 1:N){##length(t[,1])){
  ##  print(i)
  plot.all.data <-
    read.table(paste(t[i,2],t[i,1],"landscape.plot.info",sep=""))
  counts.data <- matrix(NA,nrow=1,ncol=3)
  if (file.exists(paste(t[i,2],t[i,1],"landscape.plot.info.meta",sep="")))
    counts.data <-
      read.table(paste(t[i,2],t[i,1],"landscape.plot.info.meta",sep=""))
  else {
    counts.data <- matrix(NA,nrow=1,ncol=3)
    counts.data[1,1] <- sum(plot.all.data[,3])
    counts.data[1,2] <- sum(plot.all.data[,4])
    counts.data[1,3] <- sum(plot.all.data[,5])

  }
  pos[[i]] <- plot.all.data[1,9]
  chr[[i]] <- sub("chr","",plot.all.data[1,14])
  if (name=='129xPWK'){
    s129[[i]] <- counts.data[1,2]
    spwk[[i]] <- counts.data[1,3]
  }
  else {
    spwk[[i]] <- counts.data[1,2]
    s129[[i]] <- counts.data[1,3]
  }
  total[[i]] <- sum(counts.data)
  length[[i]] <- dim(plot.all.data)[1]
}

dir.create("circos/")
circos.dump.line.plot(unlist(chr),unlist(pos),log(unlist(s129)+1),file=paste("circos/",name,".129.data",sep=""))
circos.dump.line.plot(unlist(chr),unlist(pos),log(unlist(spwk)+1),file=paste("circos/",name,".pwk.data",sep=""))
circos.dump.line.plot(unlist(chr),unlist(pos),log(unlist(total)+1),file=paste("circos/",name,".total.data",sep=""))

print(sum(unlist(total))/sum(unlist(length)))
