library('glmnet')
library('e1071')
library('class')
library('LiblineaR')
ActiveAnalyzer <- setRefClass("ActiveAnalyzer",
                              fields = list(
                                result.folder = "character",
                                gene.prefix = "vector",
                                transcript.prefix = "vector",
                                gene.names = "vector",
                                data = "matrix",
                                feature.names = "vector",
                                selected.features = "vector",
                                target.name = "character",
                                train.results = "list"
                                )
                              )

ActiveAnalyzer$methods(
initialize = function ( result.folder.temp){
  result.folder <<- result.folder.temp
  feature.names <<- c('label', 'length', 'num.exons', 'num.reads',
                     'multiple.propotion', 'num.mismatches','error.rate',
                     'num.consistent.mismatches', 'num.ag.mismatches')

  selected.features <<-c( 'num.exons', 'num.reads',
                     'multiple.propotion', 'error.rate',
                     'num.consistent.mismatches')

  target.name <<- "label"

}
)

ActiveAnalyzer$methods(
feature.extraction = function() {
  load(paste(result.folder,"/Active.Rdata",sep=""))
  gene.prefix <<- gene.prefix
  transcript.prefix <<- transcript.prefix
  gene.names <<- gene.names

  l = length(gene.names)

  length.vector=vector(mode = 'numeric', length = l)
  ne.vector = vector(mode = 'numeric', length = l)
  nr.vector = vector(mode = 'numeric', length = l)
  mp.vector = vector(mode = 'numeric', length = l)
  nm.vector = vector(mode = 'numeric', length = l)
  er.vector = vector(mode = 'numeric', length = l)
  ncm.vector = vector(mode = 'numeric', length = l)
  nag.vector = vector(mode = 'numeric', length = l)
  
  for ( i in 1:length(gene.prefix) ) {
    if (i %% 1000 == 0) {
      print(i)
      print(length(gene.prefix))

    }
    plot.info.file = paste(transcript.prefix[[i]], '.landscape.plot.info', sep = "")
    mismatcher.file = paste(gene.prefix[[i]], '.mismacher', sep="")
    meta.file = paste(gene.prefix[[i]],'.landscape.plot.meta', sep = "")
    mamf.meta.file = paste(gene.prefix[[i]], '.mamf.meta', sep = "")
    if ( ! (file.exists(plot.info.file) && file.exists(mismatcher.file) &&
            file.exists(meta.file) && file.exists(mamf.meta.file) ) ){
      next
    }
        
    plot.info <- read.table(plot.info.file, header = T , stringsAsFactors = F)
    mismatcher <- read.table(mismatcher.file, header =T , stringsAsFactors = F)
    meta <- read.table(meta.file, stringsAsFactors = F)

    mamf.meta <- read.table(mamf.meta.file, stringsAsFactors = F)
    length.vector[i] <- length(plot.info$exon_jump)
    ne.vector[i] <- sum(plot.info$exon_jump != 0)
    ## TODO: in jeweler, add column name, and here use column name to access the value
    ## meta[1,2] total reads
    ## meta[1,1] single Reads
    nr.vector[i] <- mamf.meta[1,1]
    mp.vector[i] <- (mamf.meta[1,2]-mamf.meta[1,1])/mamf.meta[1,2]
    nm.vector[i] <- sum(mismatcher$mismatches)
    er.vector[i] <- sum(mismatcher$mismatches)/sum(mismatcher$coverage)
    ## select position where read coverage > 8
    idx <- which(mismatcher$coverage > 8)
    ## at least larger than 0.6 to avoid missing SNPs
    selected <- which((mismatcher$mismatches[idx] / mismatcher$coverage[idx]) > 0.6)
    ncm.vector[i] <- sum(length(selected))
    if (length(selected) > 0){
      selected <- idx[selected]
      nag.vector[i] <- sum(mismatcher$Maternal[selected] == 'A' &&
                           mismatcher$Paternal[selected] == 'A' &&
                           (mismatcher$G[selected]/mismatcher$coverage[selected]>0.6)
                           )
    }
    else {
      nag.vector[i] <- 0
    }
    if (i %% 1000 == 0 ){
      print('saving current result')
      data <<- matrix(data=NA, nrow = length(gene.names), ncol = length(feature.names),
                      dimnames = list(
                        gene.names,
                        feature.names
                        )
                      )
      
      data[,'label'] <<- as.numeric(grepl('^Gm[0-9]+\\|',gene.names))
      data[,'length'] <<- length.vector
      data[,'num.exons'] <<- ne.vector
      data[,'num.reads'] <<- nr.vector
      data[,'multiple.propotion'] <<- mp.vector
      data[,'num.mismatches'] <<- nm.vector
      data[,'error.rate'] <<- er.vector
      data[,'num.consistent.mismatches'] <<- ncm.vector
      data[,'num.ag.mismatches'] <<- nag.vector
      data.temp = data
      save(data.temp, file = paste(result.folder, 'ActiveAnalyzer.data.Rdata', sep = ""))
    }
  }

  data <<- matrix(data=NA, nrow = length(gene.names), ncol = length(feature.names),
                      dimnames = list(
                        gene.names,
                        feature.names
                        )
                      )

  data[,'label'] <<- as.numeric(grepl('^Gm[0-9]+\\|',gene.names))
  data[,'length'] <<- length.vector
  data[,'num.exons'] <<- ne.vector
  data[,'num.reads'] <<- nr.vector
  data[,'multiple.propotion'] <<- mp.vector
  data[,'num.mismatches'] <<- nm.vector
  data[,'error.rate'] <<- er.vector
  data[,'num.consistent.mismatches'] <<- ncm.vector
  data[,'num.ag.mismatches'] <<- nag.vector
  data.temp <- data
  save(data.temp, file = paste(result.folder, 'ActiveAnalyzer.data.Rdata',sep = ""))
}
)

ActiveAnalyzer$methods(
analyze = function() {
  load(paste( result.folder, 'ActiveAnalyzer.data.Rdata', sep = ""))
  data <<- data.temp
  x <- as.matrix(data[,selected.features])
  y <- as.numeric(grepl('^Gm[0-9]+\\|',rownames(x))) ##as.vector(data[,target.name])
  rownames(x) <- 1:dim(x)[1]
  selected <- 1:18000
  x[is.na(x)]=0
  y <- y[selected]
  x <- x[selected,]
  print(paste('base line is', 1-sum(y)/length(y)))
  idx <- 1
  train.results.temp <- list()
  models <- c('libsvm') ##,'svm','lasso', 'elasticnet', 'ridge')
  for (model in models) {
    print(paste('working on model', model))
    if (model %in% c('lasso', 'elasticnet', 'ridge') ) {
      alpha <- 1
      if (model=='lasso')
        alpha <- 1
      else if (model=='elasticnet'){
        alpha <- 0.5
      }
      else if (model=='ridge'){
        alpha <- 0
      }
      result <- cv.glmnet(x,y,alpha=alpha, family='binomial', type.measure='auc')
    }
    else if (model == 'svm') {
      result <- svm(x,as.factor(y),cross=10)
    }
    else if (model == 'libsvm'){
      result <- LiblineaR(scale(x),as.factor(y),cross=10)
      print(result)
    }
    else if (model == 'knn') {
      result <- knn.cv(x,as.factor(y))
    }
    train.results.temp [[idx]] <- result
    idx <- idx + 1
  }
  train.results <<- train.results.temp
  save(models, train.results.temp,
       file = paste(result.folder, 'ActiveAnalyzer.models.Rdata', sep = ""))
}
)

ActiveAnalyzer$methods(
visualize = function() {
  load(paste( result.folder, 'ActiveAnalyzer.data.Rdata', sep = ""))
  data <<- data.temp
  feature.names <<- c(  'num.exons',
                      'multiple.propotion',
                     'num.consistent.mismatches')
  count <- 18000
  selected <- sample(1:count)[1:2000]
  plot.data <- (data[selected,feature.names])
  str(plot.data)
  pdf(paste(result.folder, "pairs.plot.pdf",sep=""), width = 7, height = 7)
  pairs(plot.data, pch = 23, bg = c('red', 'blue')[data[selected,'label']+1])
  dev.off()
}
)                       
