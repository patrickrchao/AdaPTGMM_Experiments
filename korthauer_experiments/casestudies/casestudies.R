suppressPackageStartupMessages({
  library(ExperimentHub)
  library(benchmarkfdrData2019)
  library(SummarizedBenchmark)
  library(dplyr)
  library(rlang)
  library(parallel)
  library(splines)
  library(AdaPTGMM)
})
if(packageVersion("SummarizedBenchmark") == "0.99.2"){
  stop("Summarized Benchmark version is out of date, please update with 'BiocManager::install('SummarizedBenchmark')'.")
}
# To download AdaPTGMM:
# library(devtools)
# install_github("patrickrchao/AdaPTGMM")

add_new_metrics <- function(data){
  
  hist(assay(data)[, "unadjusted"])
  if("ind_covariate" %in% colnames(rowData(data))){
    dat <- tibble(pval = assay(data)[, "unadjusted"],
                  covariate = rowData(data)$ind_covariate)
  }else if ("ind_covar_AF" %in% colnames(rowData(data))){
    dat <- tibble(pval = assay(data)[, "unadjusted"],
                  covariate = rowData(data)$ind_covar_AF)
  }else if ("ind_covar_uninf" %in% colnames(rowData(data))){
    dat <- tibble(pval = assay(data)[, "unadjusted"],
                  covariate = rowData(data)$ind_covar_uninf)
  }else if ("ind_covar_N" %in% colnames(rowData(data))){
    dat <- tibble(pval = assay(data)[, "unadjusted"],
                  covariate = rowData(data)$ind_covar_N)
  }else if ("detection" %in% colnames(rowData(data))){
    dat <- tibble(pval = assay(data)[, "unadjusted"],
                  covariate = rowData(data)$detection)
  }else if ("meanExp" %in% colnames(rowData(data))){
    dat <- tibble(pval = assay(data)[, "unadjusted"],
                  covariate = rowData(data)$meanExp)
  }else if ("rand_covar" %in% colnames(rowData(data))){
    dat <- tibble(pval = assay(data)[, "unadjusted"],
                  covariate = rowData(data)$rand_covar)
  }
  else{
    stop("Could not find covariate name")
  }
  
  # dat$SE <- rep(1,nrow(dat))
  # if("effect_size" %in% colnames(rowData(data))){
  #   dat$effect_size <-  rowData(data)$effect_size
  # }else{
  #   dat$effect_size <- qnorm(dat$pval,lower.tail = FALSE)
  # }
  
  new_design <- BenchDesign(AdaPTGMM = new_method, data = dat)
  new_res <- buildBench(new_design)
  new_res <- addPerformanceMetric(new_res,
                                  evalMetric =  c("rejections"),
                                  assay="default")
  return(new_res)
  
}

add_new_metrics_sigma_i <- function(data,cov="ind_covariate"){
  
  if(max(abs(data$pval-2*pnorm(abs(data$test_statistic),lower.tail = F)))>0.02){
    stop("p-values are not two sided")
  }
  if(cov == "ind_covariate"){
    dat <- tibble(pval = data$pval,
                  covariate = data$ind_covariate,
                  se = data$SE)
  }else if (cov == "samplesize"){
    dat <- tibble(pval = data$pval,
                  covariate = data$ind_covar_N,
                  se = data$SE)
  }else if (cov == "maf"){
    dat <- tibble(pval = data$pval,
                  covariate = data$ind_covar_AF,
                  se = data$SE)
  }
  dat <- na.omit(dat)
  new_design <- BenchDesign(AdaPTGMM = new_method_sigma_i, data = dat)
  new_res <- buildBench(new_design)
  new_res <- addPerformanceMetric(new_res,
                                  evalMetric =  c("rejections"),
                                  assay="default")
  return(new_res)
  
}

update_RDS <- function(){
  titles <- c("cbp-csaw", "h3k4me3-csaw", "human", "mouse", 
              "enigma-al-otus", "enigma-ph-otus","baxter-genus", 
              "goodrich-genus", "papa-genus", "schubert-otus",
              "enigma-al-otus-abun", "enigma-ph-otus-abun",
              "baxter-genus-mean-abun", "goodrich-genus-abun", "papa-genus-abun",
              "schubert-otus-abun",
              "human-mast-det", "human-mast-mean", "human-scdd-det",
              "human-scdd-mean", "human-wilcox-det", "human-wilcox-mean",
              "mouse-mast-det", "mouse-mast-mean", "mouse-scdd-det",
              "mouse-scdd-mean", "mouse-wilcox-det", "mouse-wilcox-mean")
  
  all_file_titles <- gsub("-benchmark", "", bfdrData$title)
  for(title in titles){

    data_file <- file.path(data_folder,paste0(title,".rds"))
    output_file <- file.path(results_folder,paste0(title,".rds"))
    
    if(!file.exists(output_file)){
      print(title)
      
       
       id <- bfdrData$ah_id[all_file_titles == title]
       cat(paste(title,id,"\n")) 
       if (!file.exists(data_file)) {
         res <- bfdrData[[id]]
         saveRDS(res,data_file)
       }else{
         cat("File found. Reading from RDS. \n")
         res <- readRDS(data_file)
       }
      
      new_res <- add_new_metrics(res)
      saveRDS(new_res,output_file)
    }
  }
  
}


update_RDS_sigma_i <- function(){
  titles <- c("brain", "mir200c","bmi")
  for(title in titles){
  
    data_file <- file.path(data_folder,paste0(title,".rds"))
    output_file <- file.path(results_folder,paste0(title,".rds"))
    if (!file.exists(data_file)) {
      if(grepl("bmi",title)){
        warning("BMI data not loaded. Please run GWAS.Rmd file to download data first.")
        return(0)
      }else if(!grepl("bmi",title)){
        # TODO: replace this
        warning("RNAseq data not loaded. Please run RNAseq-gtex.rmd and RNAseq-mir200c.Rmd files to download data first.")
        return(0)
      }

      saveRDS(res,data_file)
    }else{
      cat("File found. Reading from RDS. \n")
      res <- readRDS(data_file)
    }
    if(title == "bmi"){
      maf_file <- file.path(results_folder,paste0(title,"-maf.rds"))
      samplesize_file <- file.path(results_folder,paste0(title,"-samplesize.rds"))
      if(!file.exists(maf_file)){
        print("Running BMI maf")
        new_res <- add_new_metrics_sigma_i(res,"maf")
        saveRDS(new_res,maf_file)
      }
      if(!file.exists(samplesize_file)){
        print("Running BMI samplesize")
        new_res <- add_new_metrics_sigma_i(res,"samplesize")
        saveRDS(new_res,samplesize_file)
      }
      
    }else{
      if(!file.exists(output_file)){
        print(paste("Running",title))
        new_res <- add_new_metrics_sigma_i(res,"ind_covariate")
        saveRDS(new_res,output_file)
      }
      
    }
  }
  
 
  
}

cores <- 1

hub <- ExperimentHub()
bfdrData <- query(hub, "benchmarkfdrData2019")

new_method <- BDMethod(x = AdaPTGMM::adapt_gmm,post = function(x){print(x$args)
  x$qvals},
  params = rlang::quos(pvals=pval, x= data.frame(icov=covariate),
                       beta_formulas = paste0("splines::ns(icov, df = ", seq(2,8,2), ")"),
                       nclasses=c(2,4,6), model_type="neural",
                       alphas = seq(0.01, 0.3, 0.01)
  ))

new_method_sigma_i <- BDMethod(x = AdaPTGMM::adapt_gmm,post = function(x){print(x$args)
  x$qvals},
  params = rlang::quos(z=test_statistic, x= data.frame(x=covariate,y=se), testing="two_sided",symmetric_modeling=TRUE,
                       beta_formulas = c("splines::ns(x, df = 3 )+splines::ns(y, df = 3)","splines::ns(x, df = 5 )+splines::ns(y, df = 5)"),
                       nclasses=c(2,4,6), model_type="neural",
                       alphas = seq(0.01, 0.3, 0.01),intercept_model=F
  ))

# library(devtools)
# load_all("~/statistics_research/Mask_AdaPT/adaptMT")
# new_method <- BDMethod(x = adapt_glm,post = function(x){x$qvals},
#                        params = rlang::quos(pvals=pval, x= data.frame(icov=covariate),
#                                             pi_formulas = c(paste0("splines::ns(icov, df = ", seq(1, 9, 2), ")")),
#                                             mu_formulas = c(paste0("splines::ns(icov, df = ", seq(1, 9, 2), ")"))
#                        ))
# 
# sigma_i_formulas <- c("splines::ns(x, df = 3 )+splines::ns(y, df = 3)","splines::ns(x, df = 5 )+splines::ns(y, df = 5)")
# new_method_sigma_i <- BDMethod(x = adapt_glm,post = function(x){x$qvals},
#                        params = rlang::quos(pvals=pval,x= data.frame(x=covariate,y=se),
#                                             pi_formulas = sigma_i_formulas,
#                                             mu_formulas = sigma_i_formulas
#                        ))

results_folder <- "./glm-results-final"
data_folder <- "./data_files"
update_RDS()
update_RDS_sigma_i()

