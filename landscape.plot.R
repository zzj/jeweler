plot.landscape <- function(plot.data,parental.color,marental.color,cuffcmp,name){
  total=plot.data[,3]+plot.data[,4]+plot.data[,5]
  plot(plot.data[,2,],total,col='black',type="n",lwd=3,ylab="num of reads",
       xlab='pos')
  abline(v=which(plot.data[,6]==1),col="green",lwd=2)
  abline(v=which(plot.data[,7]==1),col="purple",lwd=2)
  lines(plot.data[,2,],total,col='black',type="l",lwd=3,ylab="num of reads",
        xlab='pos')
  lines(plot.data[,2,],plot.data[,3],col='gray',type="l",lwd=3)
  lines(plot.data[,2],plot.data[,4],col=parental.color,type="l",lwd=3)
  lines(plot.data[,2],plot.data[,5],col=marental.color,type="l",lwd=3)
  print(plot.data[which(plot.data[,6]==1),3])
  cuff <- cuffcmp[as.character(levels(plot.data[1,1])),3:4]
  title(paste(name,cuff[1,]))
}

name <- "PWKx129"
parental.color <- "red"
marental.color <- "pink"

name <- "129xPWK"
parental.color <- "pink"
marental.color <- "red"

t<-read.table(paste("result/cuffsequence/",name,"/",name,".info",sep=""),stringsAsFactors=F)
cuffcmp<-read.table(paste("result/cuffsequence/",name,"/cuffcmp.tracking",sep=""),stringsAsFactors=F)
meta <- do.call(rbind,strsplit(cuffcmp[,5],"\\|"))
gene.index <- do.call(rbind,strsplit(cuffcmp[,3],"\\|"))
gene.index <- gene.index[,1]
transcript.index <- do.call(rbind,strsplit(meta[,1],":"))
transcript.index <- transcript.index[,2]
rownames(cuffcmp)=meta[,2]
pdffolder <- paste("result/cuffsequence_parent/",name,"_plot/",sep="")

selected=(c('Kcnk9','Eif2c2','Osbpl5','Copg2'))
dir.create(pdffolder, showWarnings = FALSE, recursive = TRUE)
exon.num=mat.or.vec(length(t[,1]),1)
snp.num=mat.or.vec(length(t[,1]),1)
transcript.num=mat.or.vec(length(t[,1]),1)

exon.total <- 0
exon.has.snp <- 0
transcript.has.snp <- 0
snp.total <- 0
exon.snp.num <- NULL
for (i in 1:length(t[,1])){
  if (!file.exists(paste(t[i,2],t[i,1],"landscape.plot.info",sep=""))){
    break
  }
  plot.all.data <-
    read.table(paste(t[i,2],t[i,1],"landscape.plot.info",sep=""))
  transcripts=levels(plot.all.data[,1])
  transcript.num[i]=length(transcripts)
  if (length(transcripts)>1){
    exon.num[i] <- -1
    snp.num[i] <--1
    next
  }

  b<-hist(which(plot.all.data[,6]==1),breaks=c(0,which(plot.all.data[,7]==1)),plot=F)
  exon.has.snp <- exon.has.snp+sum(b$counts!=0)
  exon.snp.num <- c(exon.snp.num,b$counts)
  transcript.has.snp <-transcript.has.snp+sum(sum(b$counts!=0)!=0)
  exon.total <- exon.total+sum(plot.all.data[,7])
  snp.total <- snp.total+sum(plot.all.data[,6])
  exon.num[i] <- sum(plot.all.data[,7])
  snp.num[i] <- sum(plot.all.data[,6])
}

pdf(paste(name,"exon","pdf",sep="."))
exon.num <- exon.num[exon.num>=0]
hist(exon.num,breaks=min(exon.num):max(exon.num),main=name)
dev.off()
pdf(paste(name,"snp","pdf",sep="."))
snp.num <- snp.num[exon.num>=0]
hist(snp.num,breaks=min(snp.num):max(snp.num),main=name)
dev.off()
pdf(paste(name,"transcript","pdf",sep="."))
transcript.num <- transcript.num[exon.num>=0]
hist(transcript.num,breaks=min(transcript.num):max(transcript.num),main=name)
dev.off()
pdf(paste(name,"exon.snp","pdf",sep="."))
exon.snp.num <- exon.snp.num[exon.num>=0]
hist(exon.snp.num,breaks=min(exon.snp.num):max(exon.snp.num),main=name)
dev.off()



stop()
for (i in 1:length(t[,1])){
  if (sum(t[i,1]==transcript.index)==0) next
  ttt <- (gene.index[t[i,1]==transcript.index])
  if (!(ttt[1] %in% selected)){
    next
  }
  print(i)
  plot.all.data <-
    read.table(paste(t[i,2],t[i,1],"landscape.plot.info",sep=""))
  ##working on the genes that only have one transcript
  transcripts=levels(plot.all.data[,1])
  for (j in transcripts){
    plot.data <- plot.all.data[plot.all.data[,1]==as.character(j),]

    pdf(paste(pdffolder,ttt[1],".",j,".landscape.plot",sep=""),width=20)
    plot.landscape(plot.data,parental.color,marental.color,cuffcmp,name)
    dev.off()
  }
}

