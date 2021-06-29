
# Load packages and source benchmark FDR
library(tidyr)
library(dplyr)
library(ggplot2)
library(magrittr)
library(cowplot)
library(tibble)
library(ggthemes)
library(grid)
library(SummarizedBenchmark)
## load helper functions
for (f in list.files("./../R", "\\.(r|R)$", full.names = TRUE)) {
  source(f)
}

insertSource("./../R/PerformanceMetricsNew.R", package = "SummarizedBenchmark",
             functions = "tidyUpMetrics")

if(packageVersion("SummarizedBenchmark") != "0.99.2"){
  stop("Summarized Benchmark version is incorrect, please run 'devtools::install_github('areyesq89/SummarizedBenchmark', ref = 'fdrbenchmark')'.")
}

print(getwd())
result_summary_dir <- "./results-summary"
gmm_results_dir <- "./results"
print(as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))-1)
res_files <- list.files(gmm_results_dir, "\\.rds$", full.names = TRUE)
res_files <- res_files[!grepl("yeast",res_files) &!grepl("polyester",res_files)]
for (ia in 1:10) {
  print(paste("############## file",ia,"out of 10 ##############"))
  summary_file <- file.path(result_summary_dir,paste0("gmm_summary_alpha", ia, ".rds"))
  if (!file.exists(summary_file)) {
    
    print(res_files)
    res <- mclapply(res_files, function(x) {
      cat(paste(x,"\n"))
      zz <- readRDS(x)
      zz_i <- lapply(zz, `[[`, "informative")
      zz_u <- lapply(zz, `[[`, "uninformative")
      zz_i <- plotsim_standardize(zz_i, alpha = ia / 100L)
      zz_i <- dplyr::select(zz_i, rep, label,  performanceMetric, alpha, value)
      zz_u <- plotsim_standardize(zz_u, alpha = ia / 100L)
      zz_u <- dplyr::select(zz_u, rep, label,  performanceMetric, alpha, value)
      print(paste("succesfully_finished",x))
      dplyr::full_join(zz_i, zz_u,
                       by = c("rep", "label", "performanceMetric", "alpha"),
                       suffix = c(".info", ".uninfo"))
    },mc.cores= 10)
    
    names(res) <- gsub("\\.rds$", "", basename(res_files))
    saveRDS(res, summary_file)
  }
}
