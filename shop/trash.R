
tracking.selected <- c(1,2,6,8,10,11)
t <- c()

for (i in 1:length(corr.result)){
  x <- corr.result[[i]]
  if (is.null(x)) {
    t <- c(t,1)
    next
  }
  x=x[tracking.selected,tracking.selected]
  if (sum(is.na(x))>0){
    x[is.na(x)]=0
    if (x[1,2]==0 && x[3,4]==0 && x[5,6]==0) {
      t <- c(t,1 )
      next    
    }
    if (sum(is.na(cuffcompare.info[i,tracking.selected+1]))==6){
      t <- c(t,1 )
      next    
      
    }
    if (x[1,2]==0 ){
      if (is.na(cuffcompare.info[i,tracking.selected[1]+1] ) && is.na(cuffcompare.info[i,tracking.selected[2]+1] )) x[1,2]=1
      else {
        t <- c(t,1 )
        next
      }
    }
    if (x[3,4]==0 ){
      if (is.na(cuffcompare.info[i,tracking.selected[3]+1] ) && is.na(cuffcompare.info[i,tracking.selected[4]+1] )) x[3,4]=1
      else {
        t <- c(t,1 )
        next
      }
    }
    if (x[5,6]==0 ){
      if (is.na(cuffcompare.info[i,tracking.selected[5]+1] ) && is.na(cuffcompare.info[i,tracking.selected[6]+1] )) x[5,6]=1
      else {
        t <- c(t,1 )
        next
      }
    }
  }
  s<-t.test(x=c(x[1,3:6],x[2,3:6],x[3,c(1:2,5:6)],x[4,c(1:2,5:6)],x[5,1:4],x[6,1:4]),y=(c(x[1,2],x[3,4],x[5,6])),alternative="less")
  t <- c(t,(s$p.value))
}

b<-lapply(corr.result,function(x){
  ret=x[tracking.selected,tracking.selected]
  ret
})
t <- 0
for (i in b){
  if (sum(is.na(i))>0) next
  t <- t+i
}



##dump html

p <- sort.int(unlist(t),index.return=T)
for ( i in 1:length(p$ix)){
  if (grepl('^Gm',gene.list[p$ix[i]]))
    next
  write(paste("<a href=http://csbio-desktop107.cs.unc.edu/rna_seq/combined/",gene.list[p$ix[i]],".pdf>",gene.list[p$ix[i]],".pdf,","</a>  ",p$x[i],"<br> ",sep=""), append=T,sep="\n",file='diff.html')
}


## avarage distance
S <- corr.result[[1]]
id <- 0
i <- 1
for (x in corr.result){
  i <- i+1
  if (is.null(x)) next
  if (sum(is.na(x))==0){
    S <- S+x*sum(mean.result[[i]])/12
    id <- id+sum(mean.result[[i]])/12
  }
  
}
S <- S-corr.result[[1]]
S <- S/id
tracking.selected <- c(1,2,6,8,10,11)

S <- S[tracking.selected,tracking.selected]
S <- 1-S
rownames(S) <- c('CAST_1','CAST_2','PWK_1','PWK_2','WSB_1','WSB_2')
colnames(S) <- c('CAST_1','CAST_2','PWK_1','PWK_2','WSB_1','WSB_2')
tr <- fastme.bal(S)
tr <- read.tree(text="(((PWK_2:0.2737770449,PWK_1:0.2673013275):0.2218329502,(WSB_2:0.2762423346,WSB_1:0.2786957348):0.2160633933):0.1357675731,(CAST_2:0.2885666223,CAST_1:0.2910210608):0.1);")

pdf('rna_seq_tree.pdf',width=3.4,height=3.4)
plot(tr)
dev.off()

###calculate counts

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

  counts <- rep(0,12)
  for (j in 1:12){
    if (!is.list(info[[j]])) next
    counts[j]=sum(meta[[j]][1,2]+meta[[j]][1,3]+meta[[j]][1,4])
  }  counts.result[[i]] <- counts 
}


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


tracking=c(
  'result/inbred_list/jeweler/FF_0684_F_alignedto_CAST/FF_0684_F_alignedto_CAST.info',
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




