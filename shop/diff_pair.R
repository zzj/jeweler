
tracking.selected <- c(1,2,6,8,10,11)
t <- c()

global.diff.plot <- function(result.folder, name, mean.result, t1, t2){
  A <- c()
  B <- c()
  
  for (mv in mean.result){
    if (sum(mv[c(t1,t2)]==0)) next
    A <- c(A,mean(mv[t1]))
    B <- c(B,mean(mv[t2]))
  }
  A[A==0]=min(A[A>0])
  B[A==0]=min(B[B>0])
  X <- log2(A+B)
  Y <- log2(A/B)
  pdf(paste(result.folder,name,".pdf",sep=""),width=4,height=4)
  plot(X,Y,xlab='log(A+B)',ylab='log(A/B)',main=name,ylim=c(-3,3 ),pch=16,cex=0.5)
  dev.off()
}


result.folder <- 'result/inbred_list/global/'
global.diff.plot(result.folder,'CASTvsPWK',mean.result,c(2,4),c(6,8))
global.diff.plot(result.folder,'CASTvsWSB',mean.result,c(2,4),c(10,11))
global.diff.plot(result.folder,'PWKvsWSB',mean.result,c(6,8),c(10,11))
global.diff.plot(result.folder,'mixed',mean.result,c(1,6,10),c(2,8,11))

