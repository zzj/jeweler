gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": \n", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
    header=FALSE, comment.char="#", nrows = nrows,
    colClasses=c("character", "character", "character",
      "integer",
      "integer",
      "character", "character", "character",
      "character"))
  colnames(gff) = c("seqname", "source", "feature", "start",
            "end",
            "score", "strand", "frame",
            "attributes")

  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}


info<-read.table('result/cuffsequence/129xPWK/129xPWK.info',stringsAsFactors=F)
colnames(info) <-
  c('geneid','folder','gtffile','fa.seq','ma.seq','read.seq','fa.map','ma.map')


for (i in 1:dim(info)[1]){

  gtfinfo <- gffRead(info$gtffile[i])
  fa.map<-read.table(info$fa.map[i],stringsAsFactors=F,skip=5,comment.char='*')
  fa.map <- fa.map[order(fa.map[,16]),]
  print(dim(fa.map)[1])
  if (dim(fa.map)[1]>10000)
    next
  ma.map<-read.table(info$ma.map[i],stringsAsFactors=F,skip=5,comment.char='*')
  ma.map <- ma.map[order(ma.map[,16]),]

  fa.color=mat.or.vec(1,length(fa.map))
  ma.color=mat.or.vec(1,length(fa.map))
  ma.channel=mat.or.vec(1,dim(ma.map)[1])
  ma.channel.end=c(0)
  ma.channel.num=1
  fa.channel=mat.or.vec(1,dim(fa.map)[1])
  fa.channel.end=c(0)
  fa.channel.num=1
  for (k in 1:dim(fa.map)[1]){
    if (sum(ma.map[,11]==fa.map[k,11] & ma.map[,10]==fa.map[k,10] & ma.map[,14]==fa.map[k,14])){
      q=which(ma.map[,11]==fa.map[k,11] & ma.map[,10]==fa.map[k,10]& ma.map[,14]==fa.map[k,14])

      if (length(q)>1) {
        print(ma.map[q,10])
        print(ma.map[q,11])
        print('warning')
        next
      }
    }
    else {
      next
    }
    if (sum(fa.map[k,3:8])<=2 && sum(ma.map[q,3:8])<=2){
      ## same number of mismatch
      ## might be no error (mismatch=0) or read error (mismatch>0)
      if (fa.map[k,2]==ma.map[q,2]){
        fa.color[k]='gray'
        ma.color[q]='gray'
      }
      else  if (fa.map[k,2]<ma.map[q,2]){
        fa.color[k]='blue'
        ma.color[q]='white'
      }
      else{
        fa.color[k]='white'
        ma.color[q]='red'
      }
      if (sum(fa.channel.end<fa.map[k,16])>0){
        selected=min(which(fa.channel.end<fa.map[k,16]))
        fa.channel[k]=selected
        fa.channel.end[selected]=fa.map[k,17]
      }
      else {
        fa.channel.num=fa.channel.num+1
        fa.channel[k]=fa.channel.num
        fa.channel.end=c(fa.channel.end,fa.map[k,17])
      }
      if (sum(ma.channel.end<ma.map[q,16])>0){
        selected=min(which(ma.channel.end<ma.map[q,16]))
        ma.channel[q]=selected
        ma.channel.end[selected]=ma.map[q,17]
      }
      else {
        ma.channel.num=ma.channel.num+1
        ma.channel[q]=ma.channel.num
        ma.channel.end=c(ma.channel.end,ma.map[q,17])
      }
    }
  }

  xleft=min(gtfinfo$start)
  xright=max(gtfinfo$end)
  shifted.start <- 0
  shifted.end <- 0


  pdf(paste(info$folder[i],info$geneid[i],".seq.plot.pdf",sep=""))
  transcript.id=0
  for(k in 1:dim(gtfinfo)[1]){
    if (gtfinfo$feature[k]=='transcript'){
      plot(c(xleft,xright),c(0,1))
      transcript.id=transcript.id+1
      title(main=paste("transcript ", transcript.id))
      shifted.start <- 0
      shifted.end <- 0
    }
    if (gtfinfo$feature[k]=='exon'){
      segments(gtfinfo$start[k],0.5,gtfinfo$end[k],0.5)
      shifted.start <- shifted.end
      shifted.end <- shifted.end+gtfinfo$end[k]-gtfinfo$start[k]+1;
      selected <- which(fa.map[,16]<=shifted.end &
                        fa.map[,17]>=shifted.start)
      for (j in selected){

        if (fa.channel[j]==0)
          next
        fa.y=fa.channel[j]/fa.channel.num*0.5+0.5
        segments(gtfinfo$start[k]+max(0,fa.map[j,16]-shifted.start),fa.y,gtfinfo$end[k]+min(0,fa.map[j,17]-shifted.end),fa.y,col=fa.color[j])
      }
      selected <- which(ma.map[,16]<=shifted.end &
                        ma.map[,17]>=shifted.start)
      for (j in selected){
        if (ma.channel[j]==0)
          next
        ma.y=0.5-ma.channel[j]/ma.channel.num*0.5
        segments(gtfinfo$start[k]+max(0,ma.map[j,16]-shifted.start),ma.y,gtfinfo$end[k]+min(0,ma.map[j,17]-shifted.end),ma.y,col=ma.color[j])
      }
    }
  }
  transcript.id=0
  exon.id=1
  for(k in 1:dim(gtfinfo)[1]){
    if (gtfinfo$feature[k]=='transcript'){
      shifted.start <- 0
      shifted.end <- 0
      transcript.id=transcript.id+1
      exon.id=1
    }
    if (gtfinfo$feature[k]=='exon'){
      plot(c(gtfinfo$start[k],gtfinfo$end[k]),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
      title(main=paste("transcript ", transcript.id,"exon",exon.id))
      exon.id=exon.id+1
      segments(gtfinfo$start[k],0.5,gtfinfo$end[k],0.5)
      shifted.start <- shifted.end
      shifted.end <- shifted.end+gtfinfo$end[k]-gtfinfo$start[k]+1;
      selected <- which(fa.map[,16]<=shifted.end &
                        fa.map[,17]>=shifted.start)
      for (j in selected){

        if (fa.channel[j]==0)
          next
        fa.y=fa.channel[j]/fa.channel.num*0.5+0.5
        segments(gtfinfo$start[k]+max(0,fa.map[j,16]-shifted.start),fa.y,gtfinfo$end[k]+min(0,fa.map[j,17]-shifted.end),fa.y,col=fa.color[j])
      }
      selected <- which(ma.map[,16]<=shifted.end &
                        ma.map[,17]>=shifted.start)
      for (j in selected){
        if (ma.channel[j]==0)
          next
        ma.y=0.5-ma.channel[j]/ma.channel.num*0.5
        segments(gtfinfo$start[k]+max(0,ma.map[j,16]-shifted.start),ma.y,gtfinfo$end[k]+min(0,ma.map[j,17]-shifted.end),ma.y,col=ma.color[j])
      }
    }
  }

  dev.off()

}
