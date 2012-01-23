gene.pileup.plot <- function(data.folder, landscape.file, mismatcher.file, title.text){

  plot.info <- read.table(landscape.file,stringsAsFactors=F)
  plot.data <- list()
  locations <- NULL
  num.transcripts=dim(plot.info)[1]
  for (i in 1:num.transcripts){
    plot.data[[i]] <- read.table(paste(data.folder,plot.info[i,1],".landscape.plot.info",sep=""),
                                 stringsAsFactors=F,header=T)
    locations <- union(locations,plot.data[[i]]$location)
  }

  mismatcher <- read.table((mismatcher.file),header=T)
  locations <- union(locations,mismatcher$location)
  locations <- sort(locations)
  tracking.unknown <- matrix(data=0,nrow=length(locations),ncol=num.transcripts)
  tracking.total <- matrix(data=0,nrow=length(locations),ncol=num.transcripts)
  tracking.paternal <- matrix(data=0,nrow=length(locations),ncol=num.transcripts)
  tracking.maternal <- matrix(data=0,nrow=length(locations),ncol=num.transcripts)
  tracking.exon.jump <- as.vector(matrix(data=0,nrow=length(locations),ncol=1))
  tracking.snp <- as.vector(matrix(data=0,nrow=length(locations),ncol=1))
  tracking.mismatches <- vector(mode='numeric',length=length(locations))
  for (j in 1:num.transcripts){
    index <- which((locations %in% plot.data[[j]]$location))

    tracking.total[index,j]=plot.data[[j]]$unknown+plot.data[[j]]$paternal+plot.data[[j]]$maternal
    tracking.maternal[index,j]=plot.data[[j]]$maternal
    tracking.paternal[index,j]=plot.data[[j]]$paternal
    tracking.unknown[index,j]=plot.data[[j]]$unknown
    tracking.exon.jump[index]=tracking.exon.jump[index] + plot.data[[j]]$exon_jump
    tracking.snp[index]=tracking.snp[index] + plot.data[[j]]$is_snp
  }
  
  if (!is.null(mismatcher) && length(index)>10){
    index <-which (locations %in% mismatcher$location)
    tracking.mismatches[index] = mismatcher$mismatches
  }
  par(mfrow=c(num.transcripts,1),mar=c(3,3,3,3))
  ymax=max(tracking.total)
  xrange=1:length(locations)

  for (j in 1:num.transcripts){
    plot(xrange,tracking.total[,j],col='black',type="n",lwd=3,ylab="num of reads",
         xlab='pos')
    abline(v=which(tracking.snp!=0),col="green",lwd=2)
    abline(v=which(tracking.exon.jump!=0),col="purple",lwd=2)
    lines(xrange,tracking.total[,j],col='black',type="l",lwd=3)
    lines(xrange,tracking.unknown[,j],col='gray',type="l",lwd=3)
    lines(xrange,tracking.paternal[,j],col='lightblue',type="l",lwd=3)
    lines(xrange,tracking.maternal[,j],col='pink',type="l",lwd=3)
    lines(xrange,tracking.mismatches,col='brown',type="h",lwd=3)
    title(title.text)
  }

}

transcript.pileup.plot <- function(locations, plot.data, mismatcher = NULL, title.text){
  tracking.unknown <- vector(mode='numeric',length=length(locations))
  tracking.total <- vector(mode='numeric',length=length(locations))
  tracking.paternal <- vector(mode='numeric',length=length(locations))
  tracking.maternal <- vector(mode='numeric',length=length(locations))
  tracking.exon.jump <- vector(mode='numeric',length=length(locations))
  tracking.snp <- vector(mode='numeric',length=length(locations))
  ## total coverage, not only the aligned reads to the current transcript
  tracking.coverage <- vector(mode='numeric',length=length(locations))
  tracking.mismatches <- vector(mode='numeric',length=length(locations))
  if (!is.na(plot.data) ) {
    index <- which((locations %in% plot.data$location))
    tracking.total[index]=plot.data$unknown+plot.data$paternal+plot.data$maternal
    tracking.maternal[index]=plot.data$maternal
    tracking.paternal[index]=plot.data$paternal
    tracking.unknown[index]=plot.data$unknown
    tracking.exon.jump[index]=tracking.exon.jump[index] + plot.data$exon_jump
    tracking.snp[index]=tracking.snp[index] + plot.data$is_snp
    if (!is.null(mismatcher) && length(index)>10){
      index <-which (locations %in% mismatcher$location)
      tracking.coverage[index] = mismatcher$coverage
      tracking.mismatches[index] = mismatcher$mismatches
    }
  }
  ymax=max(tracking.total)
  xrange=1:length(locations)
  plot(xrange,tracking.total,col='black',type="n",lwd=3,ylab="num of reads",
       xlab='pos')
  abline(v=which(tracking.snp!=0),col="green",lwd=2)
  abline(v=which(tracking.exon.jump!=0),col="purple",lwd=2)
  lines(xrange,tracking.total,col='black',type="l",lwd=3)
  lines(xrange,tracking.unknown,col='gray',type="l",lwd=3)
  lines(xrange,tracking.paternal,col='blue',type="l",lwd=3)
  lines(xrange,tracking.maternal,col='red',type="l",lwd=3)
  lines(xrange,tracking.mismatches,col='brown',type="h",lwd=3)
  title(title.text)
}
