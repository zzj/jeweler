source('active_analyzer.R')

args=(commandArgs(TRUE))
if(length(args)==0){
  idx = 1
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
merged.alias <- unlist(read.table('../info/merged_alias',stringsAsFactors=F))
sample <- merged.alias[idx]
print(sample)

result.folder <- paste("../result/merged_list/individual_study/",sample,"/active/",sep="")

ActiveAnalyzerClass <- getRefClass("ActiveAnalyzer")
aa <- ActiveAnalyzerClass$new(result.folder)

aa$feature.extraction()
aa$analyze()
aa$visualize()
