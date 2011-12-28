t<-read.csv('result/consistent_mismatch.merged.not_SNPs.csv')
rna.editing <- list()
for (i in 1:dim(t)[1]){
  a <-
    system(paste('../fastahack/fastahack -r chr',t[i,1],':',t[i,2],' ../data/genomes/129S1_converted.fa',sep=""),intern=T)
  rna.editing[[i]] <- paste(a,"->",t[i,3],seq="")
}

write.table(table(as.factor(unlist(rna.editing))),file="result.txt")
