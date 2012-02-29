
plot.comparision <- function(data, name, title){
  pdf(paste(name,".pdf",sep=""))
  y1 <-  data[,paste('before.',name,sep="")]
  y2 <-  data[,paste('after.',name,sep="")]
  ratio <-  data[,paste('ratio.',name,sep="")]
  x <- 1:length(y1)
  plot(x,y1,xlab = "samples", ylab ="nums", main = paste(title, "mean(changed rate) =",mean(ratio)) , pch = 15, col ="blue", ylim=c(min(c(y1,y2))*0.9, max(c(y1,y2))*1.1) )
  points(x,y2, pch = 15, col ="red")
  for (i in x){
    arrows(i,y1[i],i,y2[i],length= 0.1)
  }
  dev.off()
}
log <- read.table('app.stats', header= T, stringsAsFactors =F)

plot.comparision(log,"alignments", "Number of tatal alignments")

cuff <- read.table('cuff.stats', header= T, stringsAsFactors =F)
plot.comparision(cuff,"exact", 'Number of transcripts extractly matched to a known one\n')
plot.comparision(cuff,"total", 'Number of total transcripts reported by cufflinks\n')
plot.comparision(cuff,"pseudo", 'Number of pseudo transcripts reported by cufflinks\n')
plot.comparision(cuff,"unknown", 'Number of unknown transcripts reported by cufflinks\n')
