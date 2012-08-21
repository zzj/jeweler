source('pileup.plot.R')

plot.bracet.figure <- function(filename, name, resultfolder){
  dir.create(resultfolder)
  fd <- file (filename, "r")
  data = readLines(fd, n = -1)
  n <- length(data)

  for (i in 1:n){
    line <- unlist(strsplit(data[i], "\t"))
    if (as.numeric(line[2]) > 0 ) {
      newresultfolder <- paste(resultfolder, line[1], "/", sep ="")
      dir.create(newresultfolder)
      k <- as.numeric(line[2])
      bracelet.easyplot(name, newresultfolder, line[1], line[1])
      for ( j in 1:k){
        t <- j * 2 + 1
        title = paste(line[t], line[t+1])
        bracelet.easyplot(name, newresultfolder, line[t], title)
      }
    }
  }
}


bracelet.easyplot <- function(name, resultfolder, id, title){
  output <- paste(resultfolder, id, ".pdf", sep = "")
  pdf(output, width = 20);
  gene.pileup.plot(paste('../result/new_combined/new_jeweler/',name,'/',id,'/',sep=""),
                   paste('../result/new_combined/new_jeweler/',name,'/',id,'/',id, '.landscape.plot.meta',sep=""),
                   paste('../result/new_combined/new_jeweler/',name,'/',id,'/',id, '.mismacher',sep=""),
                   title)
  dev.off()
}
plot.bracet.figure('../result/merged_list/bracelet/FG_0159_M_merged/result.bracelet', 'FG_0159_M_merged',
                   '../result/merged_list/bracelet/FG_0159_M_merged/figures/')

