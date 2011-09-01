t<-read.table("result/cuffsequence/129xPWK/129xPWK.info",stringsAsFactors=F)
cuffcmp<-read.table("result/cuffsequence/129xPWK/cuffcmp.tracking",stringsAsFactors=F)
meta <- do.call(rbind,strsplit(cuffcmp[,5],"\\|"))
rownames(cuffcmp)=meta[,2]
pdffolder <- "result/cuffsequence/129xPWK_plot/"

dir.create(pdffolder, showWarnings = FALSE, recursive = TRUE)
for (i in 1:length(t[,1])){
  print(i)
  plot.data <-
    read.table(paste(t[i,2],t[i,1],"landscape.plot.info",sep=""))
  ##working on the genes that only have one transcript
  transcripts=levels(plot.data[,1])
  if (length(transcripts)>1)
    next
  total=plot.data[,3]+plot.data[,4]+plot.data[,5]
  if (sum(total)/100<1000) next
  pdf(paste(pdffolder,t[i,1],".landscape.plot",sep=""),width=20)
  plot(plot.data[,2,],total,col='black',type="n",lwd=3,ylab="num of reads",
       xlab='pos')
  abline(v=which(plot.data[,6]==1),col="green",lwd=2)
  abline(v=which(plot.data[,7]==1),col="purple",lwd=2)
  lines(plot.data[,2,],total,col='black',type="l",lwd=3,ylab="num of reads",
        xlab='pos')
  lines(plot.data[,2,],plot.data[,3],col='gray',type="l",lwd=3)
  lines(plot.data[,2],plot.data[,4],col='red',type="l",lwd=3)
  lines(plot.data[,2],plot.data[,5],col='blue',type="l",lwd=3)
  print(plot.data[which(plot.data[,6]==1),3])
  title(cuffcmp[as.character(levels(plot.data[1,1])),3:4])
  dev.off()

}

