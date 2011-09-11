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
  lines(plot.data[,2],plot.data[,8],col='yellow',type="l",lwd=3)
  print(plot.data[which(plot.data[,6]==1),3])
  cuff <- cuffcmp[as.character(levels(plot.data[1,1])),3:4]
  title(paste(name,cuff[1,]))
}

name <- "PWKx129"
parental.color <- "blue"
marental.color <- "pink"

#name <- "129xPWK"
#parental.color <- "pink"
#marental.color <- "red"

t<-read.table(paste("result/cuffsequence/",name,"/",name,".info",sep=""),stringsAsFactors=F)
cuffcmp<-read.table(paste("result/cuffsequence/",name,"/cuffcmp.tracking",sep=""),stringsAsFactors=F)
pdffolder <- paste("result/cuffsequence_parent/",name,"_plot/",sep="")

dir.create(pdffolder, showWarnings = FALSE, recursive = TRUE)


for (i in 1:length(t[,1])){
  print(i)
  plot.all.data <-
    read.table(paste(t[i,2],t[i,1],"landscape.plot.info",sep=""))
  ##working on the genes that only have one transcript
  transcripts=levels(plot.all.data[,1])
  for (j in transcripts){
    plot.data <- plot.all.data[plot.all.data[,1]==as.character(j),]

    pdf(paste(pdffolder,t[i,1],".",j,".landscape.plot",sep=""),width=20)
    #pdf(paste(pdffolder,ttt[1],".",j,".landscape.plot",sep=""),width=20)
    plot.landscape(plot.data,parental.color,marental.color,cuffcmp,name)
    dev.off()
  }
}


#t<-read.table("result/cuffsequence/PWKx129/PWKx129.info",stringsAsFactors=F)
#cuffcmp<-read.table("result/cuffsequence/PWKx129/cuffcmp.tracking",stringsAsFactors=F)
#meta <- do.call(rbind,strsplit(cuffcmp[,5],"\\|"))
#rownames(cuffcmp)=meta[,2]
#pdffolder <- "result/cuffsequence/PWKx129_plot/"

#dir.create(pdffolder, showWarnings = FALSE, recursive = TRUE)
#for (i in 1:length(t[,1])){
#for (i in 1:length(t[,1])){
  #print(i)
  #plot.data <-
  #  read.table(paste(t[i,2],t[i,1],"landscape.plot.info",sep=""))
  ##working on the genes that only have one transcript
  #transcripts=levels(plot.data[,1])
  #if (length(transcripts)>1)
  #  next
  #total=plot.data[,3]+plot.data[,4]+plot.data[,5]
  #if (sum(total)/100<1000) next
  #pdf(paste(pdffolder,t[i,1],".landscape.plot",sep=""),width=20)
  #plot(plot.data[,2,],total,col='black',type="n",lwd=3,ylab="num of reads",
  #     xlab='pos')
  #abline(v=which(plot.data[,6]==1),col="green",lwd=2)
  #abline(v=which(plot.data[,7]==1),col="purple",lwd=2)
  #lines(plot.data[,2,],total,col='black',type="l",lwd=3,ylab="num of reads",
  #      xlab='pos')
  #lines(plot.data[,2,],plot.data[,3],col='gray',type="l",lwd=3)
  #lines(plot.data[,2],plot.data[,4],col='pink',type="l",lwd=3)
  #lines(plot.data[,2],plot.data[,5],col='blue',type="l",lwd=3)
  #lines(plot.data[,2],plot.data[,8],col='yellow',type="l",lwd=3)
  #print(plot.data[which(plot.data[,6]==1),3])
  #if(sum(plot.data[,8])>0)
  #  print(i)
  #title(cuffcmp[as.character(levels(plot.data[1,1])),3:4])
  #dev.off()

#}



#stop()
 # v=which(plot.data[,6]==1)
 # if (length(v)==0)
 # next
 # if (length(v) > 0)
 # {
 #     temp = plot.data[v,4]-plot.data[v,5]
 #     npositive <- 0
 #     nnegative <- 0
#      print(temp)
 #     for(j in 1:length(temp))
#{
#if(temp[j]>0)
#npositive <- 1
#if(temp[j]<0)
#nnegative <- 1
#if (npositive*nnegative==1)
#break
#}
 # }
#if (length(v)==1)
#{
#next
#}
#else
#{
#if(npositive*nnegative==1)
#next
#}


