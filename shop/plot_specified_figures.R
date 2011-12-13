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
  'result/inbred_list/jeweler/HH_1379_F_alignedto_WSB/HH_1379_F_alignedto_WSB.info'
)


tracking.selected <- c(1,2,10,11,6,8)
tracking.num=length(tracking.selected)
 tracking.info=list()
 id <- 1
 for (i in 1:12){
   a <- read.table(tracking[i],stringsAsFactors=F)
   rownames(a)=a[,1]
   tracking.info[[id]]=a
   id <- id+1
 }

selected <- c('Sv2c|ENSMUST00000161263','Ddx19a|ENSMUST00000040416','Cldn12|ENSMUST00000115445','Vegfa|ENSMUST00000133605')
ymaxs <- c(0,600,50,50,600)
info <- list()
prefix <- 'result/inbred_list/selected/'
n <- 1
for ( i in (which( gene.list %in% selected))){
    print(i)
    print(n)
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
  png(paste(prefix,gene.list[i],".pnp",sep=""))
  par(mfrow=c(tracking.num,1),mar=c(0.1,0.1,0.1,0.1),oma=c(1,1,1,1))

  ymax=100                                 #max(tracking.plot)               
  for (j in tracking.selected){
    print(ymax)
    plot((-locations),tracking.plot[,j],ylim=c(0,ymax),col='blue',xaxt='n', ann=FALSE,yaxt='n',type='h',axes=FALSE,yaxs='i')
    if (sum(tracking.plot[,j]>=ymax)){
      points((-locations)[tracking.plot[,j]>=ymax],rep(ymax,sum(tracking.plot[,j]>=ymax)),col='red',pch=19)
    }
    abline(h=0,col='blue')
  }
  
    ##plot.2d('CAST','PWK',tracking.plot[,c(1,2)],tracking.plot[,c(6,8)])
  ##plot.2d('CAST','WSB',tracking.plot[,c(1,2)],tracking.plot[,c(10,11)])
  ##plot.2d('PWK','WSB',tracking.plot[,c(6,8)],tracking.plot[,c(10,11)])
  dev.off()
  corr.result[[i]] <- cor(tracking.plot)
  n <- n+1
}
