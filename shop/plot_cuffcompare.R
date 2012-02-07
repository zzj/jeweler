source('cuffcompare_analysis.R')
result.file <- '../result/merged_list/cuffcompare/cuffcompare.tracking'
jeweler.result.file <- '../info/merged_jeweler_result'

CuffcomparePlotterClass <- getRefClass("CuffcomparePlotter")
## plotter <- CuffcomparePlotterClass$new(result.file = result.file,
##                                        jeweler.result.file = jeweler.result.file,
##                                        result.folder = '../result/merged_list/')

sample <- c('FG_0159_M_merged')
result.file <- paste('../result/merged_list/cuffcompare/',sample,'/cuffcompare.tracking', sep="")
plotter <- CuffcomparePlotterClass$new(result.file = result.file,
                                       jeweler.result.input = paste('result/merged_list/jeweler/',sample,'/',sample,'.info', sep = ""),
                                       is.jeweler.result.file = F,
                                        result.folder = '../result/merged_list/')

plotter$filter.transcript()
##plotter$plot()
