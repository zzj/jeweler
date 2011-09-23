plot.landscape <- function(plot.data,parental.color,marental.color,cuffcmp,name){
  total=plot.data[,3]+plot.data[,4]+plot.data[,5]
  plot(plot.data[,2,],total,col='black',type="n",lwd=3,ylab="num of reads",
       xlab='pos')
  abline(v=which(plot.data[,6]==1),col="green",lwd=2)
  abline(v=which(plot.data[,7]==1),col="purple",lwd=2)
                                        # abline(v=which(plot.data[,8]>10), col="chocolat",lwd=2)
  lines(plot.data[,2,],total,col='black',type="l",lwd=3,ylab="num of reads",
        xlab='pos')
  lines(plot.data[,2,],plot.data[,3],col='gray',type="l",lwd=3)
  lines(plot.data[,2],plot.data[,4],col=parental.color,type="l",lwd=3)
  lines(plot.data[,2],plot.data[,5],col=marental.color,type="l",lwd=3)
  lines(plot.data[,2],plot.data[,8],col='brown4',type="l",lwd=3)
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
##pdffolder <- paste("result/cuffsequence_parent/",name,"_plot/",sep="")
pdffolder <- paste("result/cuffsequence_no_novelSNP/",name,"_plot/",sep="")


dir.create(pdffolder, showWarnings = FALSE, recursive = TRUE)


for (i in 1:length(t[,1])){
  ##  print(i)
  plot.all.data <-
    read.table(paste(t[i,2],t[i,1],"landscape.plot.info",sep=""))
  ##working on the genes that only have one transcript
  transcripts=levels(plot.all.data[,1])
  for (j in transcripts){
    satisfy <- TRUE
    plot.data <- plot.all.data[plot.all.data[,1]==as.character(j),]
    for(k in 1:length(plot.data[,1]))
      {

        if(plot.data[k,8]>20)
          {
            print(i)
            print(j)
            satisfy <- FALSE
            break
          }
      }

    if(!satisfy)
      next

    cuff <- cuffcmp[as.character(levels(plot.data[1,1])),3]
    pdf(paste(pdffolder,t[i,1],".",j,".",cuff,".",ceiling(sum(plot.data[,3:5])/100),".landscape.plot.pdf",sep=""),width=20)
                                        #pdf(paste(pdffolder,ttt[1],".",j,".landscape.plot",sep=""),width=20)
    plot.landscape(plot.data,parental.color,marental.color,cuffcmp,name)
    dev.off()
  }
}
