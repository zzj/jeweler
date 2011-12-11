
plot.2d <- function(xlab,ylab,xdata,ydata){
  X <- rowMeans(xdata)
  Y <- rowMeans(ydata)
  show <- intersect(which((X!=0)),which(Y!=0))
  if (length(show))
    plot(X[show],Y[show],xlab=xlab,ylab=ylab,pch=18)
  else{
    plot(X,Y,xlab=xlab,ylab=ylab,pch=18)
  }
  abline(a=0,b=1,col='red')
}


tracking=c('result/inbred_list/jeweler/FF_0684_F_alignedto_CAST/FF_0684_F_alignedto_CAST.info',
  'result/inbred_list/jeweler/FF_0727_F_alignedto_CAST/FF_0727_F_alignedto_CAST.info',
  'result/inbred_list/jeweler/FF_0758_M_alignedto_CAST/FF_0758_M_alignedto_CAST.info',
  'result/inbred_list/jeweler/FF_6136_F_alignedto_CAST/FF_6136_F_alignedto_CAST.info',
  'result/inbred_list/jeweler/GG_1236_M_alignedto_PWK/GG_1236_M_alignedto_PWK.info',
  'result/inbred_list/jeweler/GG_1241_F_alignedto_PWK/GG_1241_F_alignedto_PWK.info',
  'result/inbred_list/jeweler/GG_1260_F_alignedto_PWK/GG_1260_F_alignedto_PWK.info',
  'result/inbred_list/jeweler/GG_6208_F_alignedto_PWK/GG_6208_F_alignedto_PWK.info',
  'result/inbred_list/jeweler/HH_1345_M_alignedto_WSB/HH_1345_M_alignedto_WSB.info',
  'result/inbred_list/jeweler/HH_1359_F_alignedto_WSB/HH_1359_F_alignedto_WSB.info',
  'result/inbred_list/jeweler/HH_1361_F_alignedto_WSB/HH_1361_F_alignedto_WSB.info',
  'result/inbred_list/jeweler/HH_1379_F_alignedto_WSB/HH_1379_F_alignedto_WSB.info')

tracking.names=c('FF_0684_F_alignedto_CAST',
  'FF_0727_F_alignedto_CAST',
  'FF_0758_M_alignedto_CAST',
  'FF_6136_F_alignedto_CAST',
  'GG_1236_M_alignedto_PWK',
  'GG_1241_F_alignedto_PWK',
  'GG_1260_F_alignedto_PWK',
  'GG_6208_F_alignedto_PWK',
  'HH_1345_M_alignedto_WSB',
  'HH_1359_F_alignedto_WSB',
  'HH_1361_F_alignedto_WSB',
  'HH_1379_F_alignedto_WSB')


tracking.selected=c(1,2,4,6,7,8,10,11,12)

tracking.selected <- 1:12
tracking.selected <- c(1,2,6,8,10,11)
tracking.num=length(tracking.selected)
tracking.info=list()
id <- 1
for (i in 1:12){
  a <- read.table(tracking[i],stringsAsFactors=F)
  rownames(a)=a[,1]
  tracking.info[[id]]=a
  id <- id+1
}

cuffcompare.result=c('result/inbred_list/cuffcompare/cuffcompare.tracking')
cuffcompare.temp.info=read.table(cuffcompare.result,stringsAsFactors=F)

gene.list=unique(cuffcompare.temp.info[,3])
load('result.Rdata')
## cuffcompare.info<-matrix(data=NA,nrow=length(gene.list),ncol=tracking.num+1,dimnames=list(gene.list,c('pattern',paste('tracking',1:12))))
## for (i in 1:dim(cuffcompare.temp.info)[1]){
##   if (is.na(cuffcompare.info[cuffcompare.temp.info[i,3],'pattern']) || cuffcompare.temp.info[i,4]=='.'){
##     cuffcompare.info[cuffcompare.temp.info[i,3],'pattern']=cuffcompare.temp.info[i,4]
##     # skip pattern coloumn
##     id <- 2
##     for (j in 1:12){
##       ##trascript name
##       if ((cuffcompare.temp.info[i,4+j])!='-')
##         cuffcompare.info[cuffcompare.temp.info[i,3],id] <- strsplit((cuffcompare.temp.info[i,4+j]),'\\|')[[1]][2] 
##       id <- id+1
##     }
##   }
## }



info <- list()
prefix <- 'result/inbred_list/new_combined/'
corr.result <- list()
mean.result <- list()
counts.result <- list()
for (i in 1:dim(cuffcompare.info)[1]){
  print(i)
  id <- 1
  info <- list()
  meta <- list()
  for (j in 2:dim(cuffcompare.info)[2]){
    if (!is.na(cuffcompare.info[i,j])){
      gene.id=paste(strsplit(cuffcompare.info[i,j],'\\.')[[1]][1:2],collapse='.')
      transcript.id=cuffcompare.info[i,j]
      if (file.exists(paste(tracking.info[[j-1]][gene.id,2],'/',transcript.id,'.landscape.plot.info',sep=""))){
        info[[id]]<-read.table(paste(tracking.info[[j-1]][gene.id,2],'/',transcript.id,'.landscape.plot.info',sep=""),header=T)
        meta[[id]]<-read.table(paste(tracking.info[[j-1]][gene.id,2],'/',gene.id,'.landscape.plot.meta',sep=""),header=F)
      }
      else{
        info[[id]]=NA
        meta[[id]]=NA
      }
    }
    else {
        meta[[id]]=NA
        info[[id]]=NA
    }
    id <- id+1
  }

  locations <- NULL
  for (j in 1:12){
    if (!is.list(info[[j]])) next
    locations <- union(locations,info[[j]]$location)
  }
  locations <- sort(locations)
  if (is.null(locations)) next
  tracking.plot <- matrix(data=0,nrow=length(locations),ncol=12)
  tracking.exon.jump <- as.vector(matrix(data=0,nrow=length(locations),ncol=1))

  meanv <- rep(0,12)
  counts <- rep(0,12)
  for (j in 1:12){
    if (!is.list(info[[j]])) next
    index <- which((locations %in% info[[j]]$location))
    tracking.plot[index,j]=info[[j]]$unknown+info[[j]]$paternal+info[[j]]$maternal
    meanv[j]=sum(tracking.plot[index,j])/length(index)
    counts[j]=sum(meta[[j]][1,2]+meta[[j]][1,3]+meta[[j]][1,4])
    if (j %in% tracking.selected)
      tracking.exon.jump[index]=tracking.exon.jump[index] + info[[j]]$exon_jump
  }
  mean.result[[i]] <- meanv
  counts.result[[i]] <- counts
  ## pdf(paste(prefix,gene.list[i],".pdf",sep=""))
  ## par(mfrow=c(tracking.num,1),mar=c(0.1,0.1,0.1,0.1))
  ## ymax=max(tracking.plot[,tracking.selected])
  ## for (j in tracking.selected){
  ##   plot((locations),tracking.plot[,j],ylim=c(0,ymax),col='blue',xaxt='n', ann=FALSE,yaxt='n',type='h',axes=FALSE)
  ## }
  ## par(mfrow=c(1,1),mar=c(6,6,6,6))
  ## plot.2d('CAST','PWK',tracking.plot[,c(1,2)],tracking.plot[,c(6,8)])
  ## plot.2d('CAST','WSB',tracking.plot[,c(1,2)],tracking.plot[,c(10,11)])
  ## plot.2d('PWK','WSB',tracking.plot[,c(6,8)],tracking.plot[,c(10,11)])
  ## dev.off()
  corr.result[[i]] <- cor(tracking.plot)
}
save(counts.result,mean.result,cuffcompare.info,gene.list,corr.result,file='result.Rdata')
## select either the only one or the mixed one 'pattern code ' is dot (.)




