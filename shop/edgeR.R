library(edgeR)
##edge R

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
  print(topTags(et))
  q <- p.adjust(et$table$p.value,'BH')
  pdf(paste(folder,name,'.pdf',sep=""),width=4,height=4)
  if (name!='mixed')
    plotSmear(y,de.tags=rownames(et$table)[et$table$p.value<0.05],cex=0.5)
  else
    plotSmear(y,cex=0.5)
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
