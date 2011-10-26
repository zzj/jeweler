circos.dump.line.plot <- function(chr, pos, data,file){
  stopifnot(length(chr)==length(pos) || length(pos)==length(data))
  output <- data.frame(chr=paste('mm',tolower(chr),sep=""),left=pos,right=pos,data=as.vector(data))
  write.table(output,row.names=F,col.names=FALSE,sep='\t',file=file,qmethod='escape',quote=FALSE)
}


getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  a<-sapply(s, function(atts) {
    a = strsplit(atts, split = " ", fixed = TRUE)
    m = match(field, sapply(a, "[", 2))
    if (!is.na(m)) {
      rv = a[[m]][3]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
  sub("\\\"","",sub("\\\"","",a))
}


gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
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

  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "),
      "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}


gfffile <- "../data/ensembl/Mus_musculus.NCBIM37.63.chr.gtf"
gff <- gffRead(gfffile)
gff$transcript.id=getAttributeField(gff$attributes,"transcript_id")


gene.list2circos.list <- function(gff,genefile,outputfile){
  a<-read.table(genefile,stringsAsFactors=F)
  b <- sapply((strsplit(a[,1],split='\\|')),"[",2)
  chr <- list()
  pos <- list()
  ones <- list()
  p <- 1
  for (i in b){
    t <- which(gff$transcript.id==i)
    chr[[p]] <- sub("chr","",gff$seqname[t[1]])
    pos[[p]] <- gff$start[t[1]]
    ones[[p]] <- 1
    p <- p+1
  }
  circos.dump.line.plot(unlist(chr),unlist(pos),unlist(ones),file=outputfile)
}
genefile <- '../result/differential_expression/relaxed_129'
outputfile <- paste("circos/",'relaxed_129',sep="")
gene.list2circos.list(gff,genefile,outputfile)
genefile <- '../result/differential_expression/relaxed_PWK'
outputfile <- paste("circos/",'relaxed_PWK',sep="")
gene.list2circos.list(gff,genefile,outputfile)
genefile <- '../result/differential_expression/relaxed_maternal'
outputfile <- paste("circos/",'relaxed_maternal',sep="")
gene.list2circos.list(gff,genefile,outputfile)
genefile <- '../result/differential_expression/relaxed_paternal'
outputfile <- paste("circos/",'relaxed_paternal',sep="")
gene.list2circos.list(gff,genefile,outputfile)
