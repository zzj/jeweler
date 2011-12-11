
tracking.selected <- c(1,2,6,8,10,11)
diffx <- function(choose){
  t <- c()
  for (i in 1:length(corr.result)){
    x <- corr.result[[i]]
    if (is.null(x)) {
      t <- c(t,1)
      next
    }
    x=x[tracking.selected,tracking.selected]
    if (sum(is.na(x))>0){
      x[is.na(x)]=0
      if (x[1,2]==0 && x[3,4]==0 && x[5,6]==0) {
        t <- c(t,1 )
        next    
      }
      if (sum(is.na(cuffcompare.info[i,tracking.selected+1]))==6){
        t <- c(t,1 )
        next    
      }
      if (x[1,2]==0 ){
        if (is.na(cuffcompare.info[i,tracking.selected[1]+1] ) && is.na(cuffcompare.info[i,tracking.selected[2]+1] )) x[1,2]=1
        else {
          t <- c(t,1 )
          next
        }
      }
      if (x[3,4]==0 ){
        if (is.na(cuffcompare.info[i,tracking.selected[3]+1] ) && is.na(cuffcompare.info[i,tracking.selected[4]+1] )) x[3,4]=1
        else {
          t <- c(t,1 )
          next
        }
      }
      if (x[5,6]==0 ){
        if (is.na(cuffcompare.info[i,tracking.selected[5]+1] ) && is.na(cuffcompare.info[i,tracking.selected[6]+1] )) x[5,6]=1
        else {
          t <- c(t,1 )
          next
        }
      }
    }
    X <- c()
    Y <- c()
    if (sum(choose==1)>0){
      X <- c(X,x[1,3:6],x[2,3:6])
      Y <- c(Y,x[1,2])
    }
    if (sum(choose==2)>0){
      X <- c(X,x[3,c(1:2,5:6)],x[4,c(1:2,5:6)])
      Y <- c(Y,x[3,4])
    }
    if (sum(choose==3)>0){
      X <- c(X,x[5,1:4],x[6,1:4])
      Y <- c(Y,x[5,6])
    }

    if (sum(Y)==2) {
      t <- c(t,1)
      next
    }
    s<-t.test(x=X,y=Y,alternative="less")
    t <- c(t,(s$p.value))
  }
  t
}
t1 <- diffx(c(1,2))
t2 <- diffx(c(1,3))
t3 <- diffx(c(2,3))
transcript.names=rownames(cuffcompare.info)
fg <- transcript.names[t1<0.01]
fg.name <- unlist(lapply(strsplit(fg,'\\|'),function(x){return (x[1])}))
fg.list <- unlist(lapply(strsplit(fg,'\\|'),function(x){return (x[2])}))
fg.list <- fg.list[grep('^Gm',fg.name,invert=T)]
fh <- transcript.names[t2<0.01]
fh.name <- unlist(lapply(strsplit(fh,'\\|'),function(x){return (x[1])}))
fh.list <- unlist(lapply(strsplit(fh,'\\|'),function(x){return (x[2])}))
fh.list <- fh.list[grep('^Gm',fh.name,invert=T)]
gh <- transcript.names[t3<0.01]
gh.name <- unlist(lapply(strsplit(gh,'\\|'),function(x){return (x[1])}))
gh.list <- unlist(lapply(strsplit(gh,'\\|'),function(x){return (x[2])}))
gh.list <- gh.list[grep('^Gm',gh.name,invert=T)]
transcript.names <- transcript.names[grep('^Gm',transcript.names,invert=T)]
transcript.list <- unlist(lapply(strsplit(transcript.names,'\\|'),function(x){return (x[2])}))
save(transcript.list,fg.list,fh.list,gh.list,file="difflist.Rdata")
