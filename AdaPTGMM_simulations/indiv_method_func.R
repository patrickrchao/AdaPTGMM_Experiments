#source("R/All_q_est_functions.R")



BH_method <- function(pvals, alpha){
  n <- length(pvals)
  khat <- max(c(0, which(sort(pvals) <= alpha * (1:n) / n)))
  which(pvals <= alpha * khat / n)
}

## Storey's BH Procedure (Storey et al. 2004)
Storey_method <- function(pvals,alpha,thr){
  n <- length(pvals)
  pi0 <- min(1, mean(pvals > thr) / (1 - thr))
  pvals[pvals > thr] <-  Inf
  khat <- max(c(0, which(sort(pvals) <= alpha * (1:n) / n / pi0)))
  which(pvals <= alpha * khat / n / pi0)
}

## Barber-Candes procedure (Barber and Candes, 2015)
BC_method <- function(pvals,alpha){
  sorted_mask_pvals <- sort(pmin(pvals, 1 - pvals))
  fdphat <- sapply(sorted_mask_pvals, function(thresh){
    (1 + sum(pvals >= 1 - thresh))/max(1, sum(pvals <= thresh))
  })
  khat <- which(fdphat <= alpha)
  if (length(khat) == 0){
    return(khat)
  } else {
    khat <- max(khat)
    phat <- sorted_mask_pvals[khat]
    return(which(pvals <= phat))
  }
}

BL_method <- function(pvals,x){
  idf <- 3
  pi0 <- (swfdr::lm_pi0(p=pvals,X=x,smooth.df=idf))$pi0
  qvals <- pi0 * p.adjust(pvals, method = "BH")
  return(qvals)
}
Solve_q_ordered_simple <- function(pvals, tau, eps, ADMM_params){
  target_num <- 5000
  n <- length(pvals)
  if (n <= target_num){
    qhat <- Solve_q_ordered(pvals, tau, eps, ADMM_params)
    return(qhat)
  }
  num_reps <- ceiling(n / target_num)
  new_pvals <- sapply(1:target_num, function(i){
    ind <- pmin((i - 1) * num_reps + (1:num_reps), n)
    pvals[ind[1]]
  })
  qhat <- Solve_q_ordered(new_pvals, tau, eps, ADMM_params)
  qhat <- rep(qhat, each = num_reps)[1:n]
  return(qhat)
}

SABHA_method <- function(pvals, qhat, alpha, tau){
  # Use the original, or estimated q as input
  n <- length(pvals)
  pvals[pvals > tau] <- Inf
  khat <- max(c(0, which(sort(qhat * pvals) <= alpha * (1:n) / n)))
  which(qhat * pvals <= alpha * khat / n)
}

## Adaptive Seqstep (Lei & Fithian, 2016)
Adaptive_SeqStep_method <- function(pvals, alpha, thr1, thr2){ # Lei & Fithian 2016's method
  # thr1 & thr2 correspond to s & lambda (with s<=lambda) in their paper
  fdphat <- thr1 / (1 - thr2) * (1 + cumsum(pvals > thr2)) / pmax(1, cumsum(pvals <= thr1))
  if(any(fdphat <= alpha)){
    khat <- max(which(fdphat <= alpha))
    return(which(pvals[1:khat] <= thr1))
  }else{
    return(NULL)
  }
}

# ## Independent Filtering (with oracle threshold)
# if (!requireNamespace("genefilter", quietly = TRUE)){
#   source("https://bioconductor.org/biocLite.R")
#   biocLite("genefilter")
# }
# library("genefilter")
#
# IF_oracle <- function(pvals, x, alpha){
#   theta_list <- seq(0, 1, 0.05)
#   R_list <- sapply(theta_list, function(theta){
#     R1 <- filtered_R(filter = x, test = pvals,
#                      theta = theta, method = "BH",
#                      alpha = alpha)
#     R2 <- filtered_R(filter = rev(x), test = rev(pvals),
#                      theta = theta, method = "BH",
#                      alpha = alpha)
#   })
#   return(max(R_list))
# }

## IHW (with oracle nbins)
ihw_oracle <- function(pvals, x, alpha){
  ihw(pvals, x, alpha, nbins = 15)
}

ash_method <- function(z,se){
  qvals <-   ashr::get_qvalue(ashr::ash(betahat=z, sebetahat=se))

}
lfdr_method <- function(pvals,x){
  qvals <- clfdr_hickswrapper(
    unadj_p = pvals, groups = IHW::groups_by_filter(x, 15))
  return(qvals)
}

fdrreg_method <- function(z,x,nulltype){

  if(!is.null(ncol(x))){
    qvals <- FDRreg::FDRreg(z=z,covars=as.matrix(x),
                            nulltype=nulltype,nmc=5000)$FDR
  }else{
    qvals <- FDRreg::FDRreg(z=z,covars=model.matrix( ~  splines::bs(x, df = 3) - 1),
                            nulltype=nulltype,nmc=5000)$FDR
  }
  return(qvals)
}
## adapted from IHWpaper::clfdr()
clfdr_hickswrapper <- function(unadj_p, groups, lfdr_estimation="fdrtool") {

  # Exclude this method if there are fewer than 200 tests within a grouping
  # (fdrtool is applied separately to each group, and throws a warning in such case)
  min_groups <- min(table(groups))
  if(min_groups< 200){
    stop(paste0("Not enough tests to apply this method. Require at least 200, only found ",min_groups,"."))
  }

  ## estimate local fdr within each stratum first
  lfdr_res <- lfdr_fit(unadj_p, groups, lfdr_estimation=lfdr_estimation)
  lfdrs <- lfdr_res$lfdr

  ## now use the rejection rule described in Cai's paper

  ## Remark:
  ## When sorting lfdrs, we break ties by pvalues so that in the end within each stratum
  ## we get monotonic adjusted p-values as a function of the p-values
  ## This is mainly needed for grenander based lfdrs, with most other
  ## lfdr estimation methods lfdr ties are not a problem usually

  o <- order(lfdrs, unadj_p)
  lfdrs_sorted <- lfdrs[o]
  fdr_estimate <- cumsum(lfdrs_sorted)/(1:length(unadj_p))
  adj_p <- rev(cummin(rev(fdr_estimate)))
  adj_p <- adj_p[order(o)]
  return(adj_p)
}


## helper function for cai
#' @importFrom fdrtool fdrtool
lfdr_fit <- function(unadj_p, group, lfdr_estimation="fdrtool"){

  pvals_list <- split(unadj_p, group)

  if (lfdr_estimation == "fdrtool"){
    lfdr_fun <- function(pv) fdrtool::fdrtool(pv, statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr

  }

  lfdr_list <- lapply(pvals_list, lfdr_fun)
  lfdrs <- unsplit(lfdr_list, group)

  fit_obj <- data.frame(pvalue=unadj_p, lfdr=lfdrs, group=group)
  fit_obj
}
