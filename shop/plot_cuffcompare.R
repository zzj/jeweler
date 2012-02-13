source('cuffcompare_analysis.R')
result.file <- '../result/merged_list/cuffcompare/cuffcompare.tracking'
jeweler.result.file <- '../info/merged_jeweler_result'

args=(commandArgs(TRUE))
if(length(args)==0){
  idx = 6
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

CuffcomparePlotterClass <- getRefClass("CuffcomparePlotter")
## plotter <- CuffcomparePlotterClass$new(result.file = result.file,
##                                        jeweler.result.file = jeweler.result.file,
##                                        result.folder = '../result/merged_list/')
merged.alias <- unlist(read.table('../info/merged_alias',stringsAsFactors=F))

sample <- merged.alias[idx]
print(sample)
result.file <- paste('../result/merged_list/cuffcompare/',sample,'/cuffcompare.tracking', sep="")
result.folder <- paste('../result/merged_list/individual_study/',sample,sep="")
plotter <- CuffcomparePlotterClass$new(result.file = result.file,
                                       jeweler.result.input = paste('result/merged_list/jeweler/',sample,'/',sample,'.info', sep = ""),
                                       is.jeweler.result.file = F,
                                       result.folder = result.folder)

plotter$filter.transcript()
##plotter$plot()

