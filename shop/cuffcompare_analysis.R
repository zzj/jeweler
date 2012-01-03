source('pileup.plot.R')

CuffcomparePlotter <- setRefClass("CuffcomparePlotter",
                                  fields = list(
                                    num.samples = "numeric",
                                    num.transcripts = "numeric",
                                    start.column = "numeric",
                                    gene.list = "vector",
                                    result.info = "data.frame",
                                    result.file = "character",
                                    result.folder="character",
                                    is.from.start = "logical",
                                    cuffcompare.info = "matrix",
                                    jeweler.result.info = "list",
                                    jeweler.result.list = "vector"
                                    )
                                  )

CuffcomparePlotter$methods(
initialize = function( result.file, jeweler.result.file, result.folder)
 {
   is.from.start <<- T

   result.info.file <- paste(result.file,".initialize.Rdata",sep="")
   if (is.from.start) {
     result.info <<- read.table(result.file, stringsAsFactors = F)
     result.info.temp <- result.info
     save(result.info.temp, file= result.info.file)
   }
   else {
     load(result.info.file)
     result.info <<- result.info.temp
   }
   gene.list <<- unique( result.info[,3] )
   start.column <<- 4
   num.samples <<- dim(result.info)[2]-4
   num.transcripts <<- length(gene.list)
   result.folder <<- result.folder

   jeweler.result.list <<-
     c(as.vector(read.table(jeweler.result.file, stringsAsFactors = F))[[1]])

   stopifnot(num.samples != length(jeweler.result.file))
   jeweler.result.info.file <- paste(result.file,".initialize.jeweler.Rdata",sep="")
   if (is.from.start) {
     for ( i in 1:num.samples){
       a <- read.table(file=paste("../",jeweler.result.list[[i]],sep=""),
                       stringsAsFactors=F)
       rownames(a)=a[,1]
       jeweler.result.info[[i]]<<-a
     }
     jeweler.result.info.temp <- jeweler.result.info
     save(jeweler.result.info.temp, file=jeweler.result.info.file)
   }
   else {
     load(jeweler.result.info.file)
     jeweler.result.info <<- jeweler.result.info.temp 
   }
 }
)
CuffcomparePlotter$methods(
get.cuffcompare.info = function ( ) {
  cuffcompare.info.file <- paste(result.file, ".get.cuffcompare.info.Rdata", sep = "")
  if (is.from.start) {
    cuffcompare.info.temp <- matrix(data=NA, nrow=length(gene.list), ncol=num.samples+1,
                               dimnames=list(
                                 gene.list,
                                 c('pattern',paste('tracking',1:num.samples))
                                 )
                               )
    for (i in 1:dim(result.info)[1]) {
      print(i)
      if (i>10) break
      transcript.name <- result.info[i,3]
      transcript.pattern <- result.info[i,4]
      if (is.na(cuffcompare.info.temp[transcript.name,'pattern']) ||
          transcript.pattern=='.'){
        cuffcompare.info.temp[transcript.name,'pattern'] <-transcript.pattern
        ## skip pattern coloumn
        id <- 2
        for (j in 1:num.samples){
          ##trascript name
          if ((result.info[i,start.column + j]) != '-'){
            cuffcompare.info.temp[transcript.name,id] <- 
              strsplit((result.info[i,start.column+j]),'\\|')[[1]][2]
          }
          id <- id+1
        }
      }
    }
    cuffcompare.info <<- cuffcompare.info.temp
    save(cuffcompare.info.temp, file = cuffcompare.info.file)
  }
  else {
    load(cuffcompare.info.file)
    cuffcompare.info <<- cuffcompare.info.temp
  }
}
)


CuffcomparePlotter$methods(
plot = function () {
  get.cuffcompare.info()
  for ( i in 2:dim(cuffcompare.info)[1]){
    info <- list()
    meta <- list()
    if (i>3) break
    print(i)
    for ( j in 2:dim(cuffcompare.info)[2]) {
      if ( !is.na(cuffcompare.info[i,j]) ) {
        transcript.id <- cuffcompare.info[i,j]
        ## get gene id from ttranscript id
        gene.id <- paste(strsplit(transcript.id,'\\.')[[1]][1:2],collapse='.')
        plot.info.file <- paste('../',jeweler.result.info[[j-1]][gene.id,2], '/',
                                transcript.id, '.landscape.plot.info',sep="")
        plot.meta.file <- paste('../',jeweler.result.info[[j-1]][gene.id,2], '/',
                                transcript.id, '.landscape.plot.meta',sep="")
        if (file.exists(plot.info.file)) {
          info[[j-1]] <- read.table(plot.info.file,header=T)
          meta[[j-1]] <- read.table(plot.info.file,header=F)
        }
        else {
          info[[j-1]] <- NA
          meta[[j-1]] <- NA
        }
      }
      else {
        info[[j-1]] <- NA
        meta[[j-1]] <- NA
      }
    }

    locations <- NULL
    for ( j in 1:num.samples){
      if (!is.list(info[[j]])) next
      locations <- union(locations,info[[j]]$location)
    }
    locations <- sort(locations)
    if (is.null(locations)) next
    pdf ( paste(result.folder,"/", rownames(cuffcompare.info)[i],".pdf", sep=""),
         height =20)
    par(mfrow=c(4,1))
    for (j in 1:num.samples)
      transcript.pileup.plot(locations, plot.data=info[[j]], jeweler.result.list[j])
    dev.off()


  }


}
)


