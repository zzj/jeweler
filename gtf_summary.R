file="data/cufflinks/129xPWK_129/transcripts.gtf"

##data <- readLines(file)
##table=strsplit(data,"\t")
t <- NULL
for (i in 1:length(table)){
  print(i)
  field <- strsplit(table[[i]][9],";")
  transcript <- strsplit(field[[1]][2],"\"")
  transcriptid <- strsplit(transcript[[1]][2],"\\.")
  t <- c(t,transcriptid[[1]][3])
}

