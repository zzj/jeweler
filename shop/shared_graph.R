
## this script will be run under main folder instead of shop folder

source('shop/pileup.plot.R')


args=(commandArgs(TRUE))
if(length(args)==0){
  name = "FH_0096_M_merged"
  idx = 6
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

plot.bracet.figure <- function(filename, name, resultfolder){
  print(filename)
  fd = file(filename, "r")
  print(isOpen(fd,"r"))

  data <- readLines(fd)
  n <- length(data)
  print(n)
  if (n>0){
    
    for (i in 1:n){
      line <- unlist(strsplit(data[i], "\t"))
      print(i)
      resultfolder = line[1]
      k <- length(line)
      for ( j in 2:k){
        title = paste(line[j])
        bracelet.easyplot(name, resultfolder, unlist(strsplit(line[j],"\\|"))[1], line[j])
      }
    }
  }
}

root = 'result/new_merged_list/'
bracelet.easyplot <- function(name, resultfolder, id, title){
  output <- paste(resultfolder,"/", id, ".pdf", sep = "")
  pdf(output, width = 14);
  gene.pileup.plot(paste(root, 'jeweler/',name,'/',id,'/',sep=""),
                   paste(root, 'jeweler/',name,'/',id,'/',id, '.landscape.plot.meta',sep=""),
                   paste(root, 'jeweler/',name,'/',id,'/',id, '.mismatcher',sep=""),
                   title)
  dev.off()
}
system("whoami")
plot.bracet.figure(paste(root, '/shared_graph/',name,'/shared_graph',sep=""),
                   name,
                   paste(root, '/shared_graph/',name,'/figures/',sep = ""))

