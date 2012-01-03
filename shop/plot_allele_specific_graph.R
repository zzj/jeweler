source('pileup.plot.R')

info.file <- "result/merged_list/jeweler//FG_0168_F_merged/FG_0168_F_merged.info"
data.folder <- "result/merged_list/jeweler//FG_0168_F_merged/"
output.main.folder <- "result/merged_list/allele_specific/"
output.folder <- paste(output.main.folder,gsub(".info","",basename(info.file)),"/",sep="")
dir.create(output.folder,recursive=T)


info <- read.table(info.file,stringsAsFactors=F)
for ( i in 1:dim(info)[1]){
  path.file=paste(a[i,3],".allele.specific.path",sep="")
  graph.file=paste(a[i,3],".allele.specific.graph",sep="")
  graph.info.file=paste(a[i,3],".allele.specific.info",sep="")
  pileup.file=paste(a[i,3],".landscape.plot.meta",sep="")

  if (file.exists(graph.info.file)){
    if (file.info(path.file)$size!=0){
      is.allele.specific=read.table(graph.info.file,stringsAsFactors=F)
      if (is.allele.specific[1]=='Yes'){
        system(paste('dot -Tpng ',graph.file, " -o ",output.folder,a[i,1],".png", sep=""))
        result <- paste(output.folder,a[i,1],".plot.pdf",sep="")
        pdf(result,width=20)
        gene.pileup.plot(paste(data.folder,a[i,1],"/",sep=""),
                         pileup.file, title.text)
        dev.off()
      }
    }
  }
}
