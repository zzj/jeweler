source('active_analyzer.R')


result.folder <- "../result/merged_list/active/"

ActiveAnalyzerClass <- getRefClass("ActiveAnalyzer")
aa <- ActiveAnalyzerClass$new(result.folder)

##aa$feature.extraction()
aa$analyze()
aa$visualize()
