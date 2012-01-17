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
                                    joint.plot.result.folder = "character",
                                    active.transcript.result.folder = "character",
                                    is.from.start = "logical",
                                    cuffcompare.info = "matrix",
                                    jeweler.result.info = "list",
                                    jeweler.result.list = "vector"
                                    )
                                  )

CuffcomparePlotter$methods(
initialize = function( result.file, jeweler.result.file, result.folder)
 {

   print("Initializing CuffcomparePlotter")
   is.from.start <<- F

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
   joint.plot.result.folder <<- paste(result.folder, "/joint_plot/",sep="")
   active.transcript.result.folder <<- paste(result.folder, "/active/",sep="")
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
  print("Get CuffcompareInfo")
  cuffcompare.info.file <- paste(result.file, ".get.cuffcompare.info.Rdata", sep = "")
  if (is.from.start) {
    cuffcompare.info.temp <- matrix(data=NA, nrow=length(gene.list), ncol=num.samples+1,
                               dimnames=list(
                                 gene.list,
                                 c('pattern',paste('tracking',1:num.samples))
                                 )
                               )
    for (i in 1:dim(result.info)[1]) {
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
filter.transcript = function() {
  get.cuffcompare.info()
  gene.prefix <- vector()
  transcript.prefix <- vector()
  gene.names <- vector()
  ## stupid R, the performance is bad, need to hand working on it for append function
  cache.size = 10000
  temp.t.vector <- vector(mode = 'character', length = cache.size)
  temp.g.vector <- vector(mode = 'character', length = cache.size)
  temp.names.vector <- vector(mode = 'character', length = cache.size)
  idx <- 1

  for ( i in 2:dim(cuffcompare.info)[1]){
    print(i)
    for ( j in 2:dim(cuffcompare.info)[2]) {
      if ( !is.na(cuffcompare.info[i,j]) ) {
        transcript.id <- cuffcompare.info[i,j]
        gene.id <- paste(strsplit(transcript.id,'\\.')[[1]][1:2],collapse='.')
        data.folder <- jeweler.result.info[[j-1]][gene.id,2]
        plot.info.file <- paste('../', data.folder, '/',
                                transcript.id, '.landscape.plot.info',sep="")
        temp.t <- paste('../', data.folder, '/', transcript.id,sep="")
        temp.g <- paste('../', data.folder, '/', gene.id,sep="")
        temp.names <- rownames(cuffcompare.info)[i]
        if (file.exists(plot.info.file)) {
          temp.t.vector[idx]=temp.t
          temp.g.vector[idx]=temp.g
          temp.names.vector[idx]=temp.names
          idx <- idx + 1
          if (idx == cache.size + 1){
            idx=1
            gene.prefix <- c(gene.prefix, temp.g.vector)
            transcript.prefix <- c(transcript.prefix, temp.t.vector)
            gene.names <- c(gene.names,temp.names.vector)
          }
        }
      }
    }
    if (i %% 1000 == 0) 
      save(gene.names, gene.prefix,transcript.prefix,
           file = paste(active.transcript.result.folder,"Active.Rdata",sep=""))
  }
  if (idx > 1) {
    gene.prefix <- c(gene.prefix, temp.g.vector[1:(idx-1)])
    transcript.prefix <- c(transcript.prefix, temp.t.vector[1:(idx-1)])
    gene.names <- c(gene.names,temp.names.vector[1:(idx-1)])
  }
  save(gene.names, gene.prefix,transcript.prefix,
       file = paste(active.transcript.result.folder,"Active.Rdata",sep=""))
}
)

CuffcomparePlotter$methods(
plot = function () {
  get.cuffcompare.info()
  for ( i in 29370:dim(cuffcompare.info)[1]){
    info <- list()
    meta <- list()
    mismatcher <- list()
    print(i)
    print(rownames(cuffcompare.info)[i])
    for ( j in 2:dim(cuffcompare.info)[2]) {
      if ( !is.na(cuffcompare.info[i,j]) ) {
        transcript.id <- cuffcompare.info[i,j]
        gene.id <- paste(strsplit(transcript.id,'\\.')[[1]][1:2],collapse='.')
        data.folder <- jeweler.result.info[[j-1]][gene.id,2]
        plot.info.file <- paste('../', data.folder, '/',
                                transcript.id, '.landscape.plot.info',sep="")
        plot.meta.file <- paste('../', data.folder, '/',
                                gene.id, '.landscape.plot.meta',sep="")
        plot.mismatcher.file <- paste('../', data.folder, '/',
                                gene.id, '.mismacher',sep="")
        if (file.exists(plot.info.file)) {
          info[[j-1]] <- read.table(plot.info.file,header=T)
          meta[[j-1]] <- read.table(plot.meta.file,header=F)
          ## skip the heaer row due to the bug in mismatcher.cpp
          mismatcher[[j-1]] <- read.table(plot.mismatcher.file, header=T)
        }
        else {
          info[[j-1]] <- NA
          meta[[j-1]] <- NA
          mismatcher[[j-1]] <- NA
        }
      }
      else {
        info[[j-1]] <- NA
        meta[[j-1]] <- NA
        mismatcher[[j-1]] <- NA
      }
    }

    locations <- NULL
    for ( j in 1:num.samples){
      if (!is.list(info[[j]])) next
      locations <- union(locations,info[[j]]$location)
      locations <- union(locations,mismatcher$location)
    }
    locations <- sort(locations)
    if (is.null(locations)) next
    pdf(paste(joint.plot.result.folder,"/", rownames(cuffcompare.info)[i],".pdf",
              sep=""),
        height =20)
    par(mfrow=c(4,1))
    for (j in 1:num.samples){
      transcript.pileup.plot(locations = locations,
                             plot.data = info[[j]],
                             mismatcher = mismatcher[[j]],
                             title.text = jeweler.result.list[j])
    }
    dev.off()
  }
}
)


