info.file <- "result/merged_list/jeweler//FG_0168_F_merged/FG_0168_F_merged.info"
data.folder <- "result/merged_list/jeweler//FG_0168_F_merged/"
output.main.folder <- "result/merged_list/allele_specific/"
output.folder <- paste(output.main.folder,gsub(".info","",basename(info.file)),"/",sep="")
dir.create(output.folder,recursive=T)

landscape.plot <- function(data.folder, landscape.file, result){

  plot.info <- read.table(landscape.file,stringsAsFactors=F)
  plot.data <- list()
  locations <- NULL
  num.transcripts=dim(plot.info)[1]
  for (i in 1:num.transcripts){
    plot.data[[i]] <- read.table(paste(data.folder,plot.info[i,1],".landscape.plot.info",sep=""),
                                 stringsAsFactors=F,header=T)
    locations <- union(locations,plot.data[[i]]$location)
  }
  locations <- sort(locations)
  tracking.unknown <- matrix(data=0,nrow=length(locations),ncol=num.transcripts)
  tracking.total <- matrix(data=0,nrow=length(locations),ncol=num.transcripts)
  tracking.paternal <- matrix(data=0,nrow=length(locations),ncol=num.transcripts)
  tracking.maternal <- matrix(data=0,nrow=length(locations),ncol=num.transcripts)
  tracking.exon.jump <- as.vector(matrix(data=0,nrow=length(locations),ncol=1))
  tracking.snp <- as.vector(matrix(data=0,nrow=length(locations),ncol=1))
  for (j in 1:num.transcripts){
    index <- which((locations %in% plot.data[[j]]$location))

    tracking.total[index,j]=plot.data[[j]]$unknown+plot.data[[j]]$paternal+plot.data[[j]]$maternal
    tracking.maternal[index,j]=plot.data[[j]]$maternal
    tracking.paternal[index,j]=plot.data[[j]]$paternal
    tracking.unknown[index,j]=plot.data[[j]]$unknown
    tracking.exon.jump[index]=tracking.exon.jump[index] + plot.data[[j]]$exon_jump
    tracking.snp[index]=tracking.snp[index] + plot.data[[j]]$is_snp
  }
  
  pdf(result,width=20)
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
    title(result)
  }
  dev.off()
}

info <- read.table(info.file,stringsAsFactors=F)
for ( i in 1:dim(info)[1]){
  path.file=paste(a[i,3],".allele.specific.path",sep="")
  graph.file=paste(a[i,3],".allele.specific.graph",sep="")
  graph.info.file=paste(a[i,3],".allele.specific.info",sep="")
  landscape.file=paste(a[i,3],".landscape.plot.meta",sep="")

  if (file.exists(graph.info.file)){
    if (file.info(path.file)$size!=0){
      is.allele.specific=read.table(graph.info.file,stringsAsFactors=F)
      if (is.allele.specific[1]=='Yes'){
        system(paste('dot -Tpng ',graph.file, " -o ",output.folder,a[i,1],".png", sep=""))
        landscape.plot(paste(data.folder,a[i,1],"/",sep=""),
                       landscape.file, paste(output.folder,a[i,1],".plot.pdf",sep=""))
      }
    }
  }
}

  

