suppressPackageStartupMessages({
  library(ExperimentHub)
  library(benchmarkfdrData2019)
  library(SummarizedBenchmark)
  library(dplyr)
  library(ggplot2)
  library(rlang)
  library(parallel)
})
if(packageVersion("SummarizedBenchmark") == "0.99.2"){
  stop("Summarized Benchmark version is out of date, please update with 'BiocManager::install('SummarizedBenchmark')'.")
}

add_new_metrics <- function(data){
  if(class(data)=="SummarizedBenchmark"){
    dat <- tibble(pval = assay(data)[, "unadjusted"],
                  covariate = rowData(data)$ind_covariate,
                  qvalue = rowData(data)[,"qvalue"])
    
    new_design <- BenchDesign(AdaPTGMM = new_method, data = dat)
    new_res <- buildBench(new_design,truthCols = "qvalue")
    new_res <- addPerformanceMetric(new_res,
                                    evalMetric =  c("TPR", "FDR", "TNR", "FNR", "rejections"),
                                    assay = "qvalue")
    cat(paste(format(Sys.time(), "%a %b %d %X %Y")),"\n")
    if(class(new_res) != "SummarizedBenchmark"){
      print("@@@@@@@@@@@@@@ FDR Method did not work. @@@@@@@@@@@")
      stop("FDR Method did not work.")
      return("ERROR ENCOUNTERED IN METHOD")
    }
    return(new_res)
  }else{
      print("@@@@@@@@@@@@@@ Dataset is not SummarizedBenchmark. @@@@@@@@@@@")
      stop("Dataset is not SummarizedBenchmark.")
    return("DATASET NOT SUMMARIZED BENCHMARK")
  }
   
  
}
# 90 Min per file, 9 times slower than just BH
# BH took about 11 hours, thus 100 total hours or 4.2 days

#46-79
#98-171

#i = 148
all_indices <- c(46:79,98:171)
#all_indices_sub <- all_indices[c(1,9,79:84,97:108)]
update_RDS <- function(save_data_file = T){
  
  
  for(id in all_indices){
    
    title <- bfdrData$title[id]
      
      data_file <- file.path(data_folder,paste0(title,".rds"))
      output_file <- file.path(results_folder,paste0(title,".rds"))
      
      if(!file.exists(output_file)){
        
       # id <- bfdrData$ah_id[all_file_titles == title]
        cat(paste(title,id,"\n")) 
        if (!file.exists(data_file)) {
          res <- bfdrData[[id]]
          if(save_data_file){
            saveRDS(res,data_file)
          }
        }else{
          cat("File found. Reading from RDS. \n")
          res <- readRDS(data_file)
        }
        if(is.vector(res)){
          #If it is a nested vector with informative and uninformative
          if(is.vector(res[[1]])){
            new_res_inf <- mclapply(1:length(res),function(x){
              print(paste("Informative Experiment",title,x))
              add_new_metrics(res[[x]]$informative)},mc.cores = cores)
            new_res_uninf <- mclapply(1:length(res),function(x){
              print(paste("Uninformative Experiment",title,x))
              add_new_metrics(res[[x]]$informative)},mc.cores = cores)
            new_res <- list()
            for(i in 1:length(new_res_inf)){
              new_res[[i]] <- list()
              new_res[[i]]$informative <- new_res_inf[[i]]
              new_res[[i]]$uninformative <- new_res_uninf[[i]]
            }
          }else{
            new_res <- mclapply(1:length(res),function(x){
              print(paste("Experiment",title,x))
              add_new_metrics(res[[x]])},mc.cores = cores)
          }
          
        }else{
          new_res <- add_new_metrics(res)
        }
      
       saveRDS(new_res,output_file)
    }
  }
}


cores <- 10

results_folder <- "./results"
data_folder <- "./data_files"

hub <- ExperimentHub()
bfdrData <- query(hub, "benchmarkfdrData2019")

library(splines)
library(AdaPTGMM)
new_method <- BDMethod(x = AdaPTGMM::adapt_gmm,post = function(x){print(x$args)
  x$qvals},
  params = rlang::quos(pvals=pval, x= data.frame(icov=covariate), 
                       beta_formulas = paste0("splines::ns(icov, df = ", seq(2,8,2), ")"),
                       nclasses=c(2,4,6), model_type="neural",
                       alphas = seq(0.01, 0.3, 0.01)
  ))


# Set to not save datasets to save space, can change to true if desired
update_RDS(save_data_file = F)
