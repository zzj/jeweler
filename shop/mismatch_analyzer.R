
pnames <- read.table("../result/merged_list/mismatch_analyzer/output",stringsAsFactors=F)
for(p in pnames[[1]]){
  print(p)
  name <- p
  result.folder <- paste("../result/merged_list/mismatch_analyzer/",name,"/",sep="")
  

  if (!file.exists(paste(result.folder,'result.error.first', sep = "")))
    next
  print('herer')
  first <- read.table(paste(result.folder,'result.error.first', sep = ""),stringsAsFactors = F, header = T)
  last <- read.table(paste(result.folder,'result.error.last', sep = ""), stringsAsFactors = F, header = T )
  first <- first[-1,]
  last <- last[-1,]
  pdf(paste(name,".pdf",sep=""))
  yr <- c((first$num.error),(last$num.error), (first$num.quality),(last$num.quality))
  ymin <- min(yr)
  ymax <- max(yr)
  x <- as.numeric(rownames(first))
  
  y1 <- (first$num.quality)
  y2 <- (last$num.quality)
  ## par(mar=c(5,4,4,5)+.1)
  ## plot(x, y1, col='blue', type = 'l', xlab = 'Quality Score',ylab='Num of Calls',log = 'y')
  ## lines(x, y1, col='blue', type = 'o',pch=0)
  ## lines(x, y2, col='red', type = 'o',pch = 0)
  ## par(new=TRUE)
  y1 <- (first$num.error)
  y2 <- (last$num.error)
  y3 <- 10^(-x/10)
  plot(x, y2, col='blue', type = 'n', yaxt = 'n', xaxt = 'n', xlab="", ylab="" ,log='y')
  lines(x, y1, col='blue', type = 'o')
  lines(x, y2, col='red', type = 'o')
  lines(x, y3, col='black', type = 'o')
  axis(4)
  mtext("error.rate",side=4,line=3)
  legend("topleft",col=c("red","blue","black"),lty=c(1,1,1),pch=c(1,1,1),legend=c("error rate after correction","error rate before correction","expected error rate"))
  
  dev.off()
}
