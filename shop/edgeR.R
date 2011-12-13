library(edgeR)
##edge R
plotSmear <- function (object, pair = NULL, de.tags = NULL, xlab = "logConc",
          ylab = "logFC", pch = 19, cex = 0.2, smearWidth = 0.5, panel.first = grid(),
          smooth.scatter = FALSE, lowess = FALSE, ...)
  {
    if (!(class(object) %in% c("DGEList", "DGELRT", "DGEExact")))
      stop("Currently only supports DGEList/DGELRT/DGEExact objects as the object argument.")
    if (is(object, "DGEList") & is.null(object$samples$group))
      stop("Cannot produce a smear plot if no experimental groups are defined. Here, d$samples$groups is NULL.\n")
    if (is(object, "DGEList")) {
      levs.group <- levels(object$samples$group)
      if (length(levs.group) == 1)
        stop("Cannot produce an MA-plot with only one group. The one group defined is: ",
             levs.group)
      if (is.null(pair))
        pair <- levs.group[1:2]
      if (!all(pair %in% levs.group))
        stop("At least one element of given pair is not a group.\n Groups are: ",
             paste(levs.group, collapse = " "), "\n")
      stopifnot(length(pair) == 2)
      cols1 <- pair[1] == object$samples$group
      cols2 <- pair[2] == object$samples$group
      lib.size <- object$samples$lib.size * object$samples$norm.factors
      x <- rowMeans(object$counts[, cols1, drop = FALSE]/expandAsMatrix(lib.size[cols1],
                                             dim(object$counts[, cols1, drop = FALSE])))
      y <- rowMeans(object$counts[, cols2, drop = FALSE]/expandAsMatrix(lib.size[cols2],
                                             dim(object$counts[, cols2, drop = FALSE])))
      i <- match(de.tags, rownames(object$counts))
      i <- i[!is.na(i)]
      maPlot(x, y, xlab = xlab, ylab = ylab, pch = pch, cex = cex,
             smearWidth = smearWidth, de.tags = i, panel.first = panel.first,
             smooth.scatter = smooth.scatter, lowess = lowess,
             ...)
    }
    else {
      if (is.null(object$table$logFC))
        stop("table$logFC slot in DGELRT object is NULL. We cannot produce an MA (smear) plot if more than one coefficient from the GLM is being tested in the likelihood ratio test as this results in more one logFC value per gene---one for each coefficient.\n")
      i <- match(de.tags, rownames(object$table))
      i <- i[!is.na(i)]
      maPlot(x = NULL, y = NULL, logAbundance = object$table$logConc,
             logFC = object$table$logFC, xlab = xlab, ylab = ylab,
             pch = pch, cex = cex, smearWidth = smearWidth, de.tags = i,
             panel.first = panel.first, smooth.scatter = smooth.scatter,
             lowess = lowess, ...)
    }
  }
ctable<-do.call(rbind,counts.result)

rownames(ctable)<-paste('Tag.',1:dim(ctable)[1])

plot.edgeR <- function(folder,name,selected,group){
  s <- c()
  data <- ctable[,selected]

  for (i in 1:dim(data)[1]){
    if (sum(data[i,group==1])==0) next
    if (sum(data[i,group==2])==0) next
    if (sum(data[i,]==0)>0) next
    s <- c(s,i)
  }
  data <- data[s,]
  y <- DGEList(counts=data,group=group)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  et <- exactTest(y)
  q <- p.adjust(et$table$p.value,'BH')
  pdf(paste(folder,name,'.pdf',sep=""),width=4,height=4)
  print(name)
  print(dim(et$table))
  print(sum(et$table$logFC[et$table$p.value<0.05]>0))
  print(sum(et$table$logFC[et$table$p.value<0.05]<0))
  
  if (name!='mixed'){
    if (name=='CASTvsWSB'){
      plotSmear(y,de.tags=rownames(et$table)[et$table$p.value<0.05],cex=0.5,ylim=c(6,-6),ylab='logFC')
    }
    else
      plotSmear(y,de.tags=rownames(et$table)[et$table$p.value<0.05],cex=0.5,ylim=c(6,-6),ylab='')
  }
  else
    plotSmear(y,cex=0.5,ylim=c(6,-6),ylab='')
  dev.off()
}

result.folder <- 'result/inbred_list/global/'
plot.edgeR(result.folder,'CASTvsPWK',c(1,2,6,8),c(1,1,2,2))
plot.edgeR(result.folder,'CASTvsWSB',c(1,2,10,11),c(1,1,2,2))
plot.edgeR(result.folder,'PWKvsWSB',c(6,8,10,11),c(1,1,2,2))
plot.edgeR(result.folder,'mixed',c(1,6,10,2,8,11),c(1,1,1,2,2,2))

## plot.edgeR(result.folder,'CASTvsPWK_with_zero_reads',c(1,2,6,8),c(1,1,2,2))
## plot.edgeR(result.folder,'CASTvsWSB_with_zero_reads',c(1,2,10,11),c(1,1,2,2))
## plot.edgeR(result.folder,'PWKvsWSB_with_zero_reads',c(6,8,10,11),c(1,1,2,2))
## plot.edgeR(result.folder,'mixed_with_zero_reads',c(1,6,10,2,8,11),c(1,1,1,2,2,2))
