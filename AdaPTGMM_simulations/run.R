
library(ggplot2)
library(devtools)
library(dplyr)
library(AdaPTGMM)
#install_github("patrickrchao/adaptMT")
library(adaptMT)
library(ashr)
library(fdrtool)
library(FDRreg)
library(swfdr)
library(IHW)
library(tidyr)
library(akima)
library(VGAM)
library(RColorBrewer)
library(RadaFDR)
library(parallel)
source('color_scheme.R')
source('plot_results.R')
source('generate_data.R')
source('run_method.R')
source('run_simulations.R')
source('indiv_method_func.R')

num_hypo <- 3000
num_sims <- 100
radius <- 1
testing <- "two_sided"
se <- rep(1,num_hypo)
# number of cores in parallel processing
num_cores <- 10

all_methods <- c("AdaPT","AdaPTg","AdaPTGMM",
                 "ASH",
                 "Boca Leek",
                 "BH",
                 "FDRreg-t",
                 "Storey BH",
                 "IHW","LFDR","AdaFDR"
)

dir.create("Images/", showWarnings = FALSE)
dir.create("logs/", showWarnings = FALSE)
if(testing !="two_sided"){
  chosen_methods <- all_methods[!(all_methods%in% c("FDRreg-t","ASH"))]
}else{
  chosen_methods <- all_methods
}
log <- run_simulations(chosen_methods, num_sims,num_hypo,radius,testing,se,num_cores)
plot_results(log,testing)
