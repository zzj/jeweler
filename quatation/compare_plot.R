plot.comparision <- function(data, name, filename, title, so, founder.names,
                             label.position, break.position, offset.position){
  pdf(paste(filename, ".pdf", sep=""),width=8)
  y1 <-  data[so,paste('before.', name, sep="")]
  y2 <-  data[so,paste('after.', name, sep="")]
  ratio <-  data[,paste('ratio.', name, sep="")]
  x <- 1:length(y1) + offset.position
  l = max(c(y1,y2)) - min(c(y1,y2))

  plot(x, y1, xlab = "samples", ylab ="number",
       main = paste(title, "mean(change) =", sprintf("%.2f%%", mean(ratio)*100)) ,
       pch = 15, col ="blue", ylim=c(min(c(y1,y2)) -l *0.1, max(c(y1,y2)) + l *0.15),
       xaxt="n")
  points(x,y2, pch = 20, col ="red")
  for (i in 1:length(x)){
    arrows(x[i],y1[i],x[i],y2[i], length= 0.1)
  }
  for ( i in 1:(length(break.position)-1)){
    abline(v = break.position[i], lty= 4)
  }
  axis(1, at=label.position, labels=founder.names, las=1, cex.axis=0.75)
  print(founder.names)
  legend('topleft', c('Original Pipeline','Gene Scissors Pipeline'), pch = c(15,20), col= c("blue","red"),bg = "white")
  dev.off()
}
##log <- read.table('app.stats', header= T, stringsAsFactors =F)
##plot.comparision(log,"alignments", "Number of total alignments")

order.founder = c("FF", "HH", "GG", "FG","FG", "GF", "GF", "FH", "FH", "HF", "HF", "GH","GH" ,"HG", "HG")
order.founder.alias = c("CC", "PP", "WW", "CP","CP", "PC", "PC", "CW", "CW", "WC", "WC", "PW","PW" ,"WP", "WP")
order.gender = c("F", "F", "F", "F","M","F","M","F","M","F","M","F","M","F","M")
order.names = c()


cuff <- read.table('new.cuff.stats', header=T, stringsAsFactors=F)
order.sample <- c()
for (i in 1:dim(cuff)[1]){
  for (j in 1:length(order.founder)){
    if ( substr(cuff[i,1],1,2) == order.founder[j] &&  substr(cuff[i,1],9,9) == order.gender[j]){
      order.sample <- c(order.sample,j)
      break
    }
  }
}
plotorder = sort.int(order.sample,index.return=T)
sample.order = plotorder$ix
label.position = c()
break.position = c()
offset.position = plotorder$x
for (j in 1:length(order.founder)){
  label.position = c(label.position, mean(which(plotorder$x==j)) + j)
  break.position = c(break.position, max(which(plotorder$x==j)) + j +1)
  order.names=c(order.names,paste(order.founder.alias[j], "\n", order.gender[j],"", sep= ""))
}
print(break.position)
print(label.position)
plot.comparision(cuff,"exact", "exact_new", '(b) Number of  reported transcripts that exactly or partially matched \n an annotated one\n',
                 sample.order,
                 order.names,
                 label.position,
                 break.position,
                 offset.position
                 )
plot.comparision(cuff,"pseudo", "pseudo_new", '(c) Number of annotated pseudogenes reported\n',
                 sample.order,
                 order.names,
                 label.position,
                 break.position,
                 offset.position
                 )
plot.comparision(cuff,"unknown", "unknown_new", '(d) Number of unannotated transcripts reported \n',
                 sample.order,
                 order.names,
                 label.position,
                 break.position,
                 offset.position
                 )
plot.comparision(cuff,"total", "total_new",'(a) Number of total transcripts reported\n',
                 sample.order,
                 order.names,
                 label.position,
                 break.position,
                 offset.position
                 )
