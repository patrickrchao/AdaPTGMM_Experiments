

calc_tpr <- function(rejections, h,alphas){

  if(length(rejections) == 0){
    return(rep(0,length(alphas)))
  }

  total_nonnull <- sum(h)
  tpr <- unlist(lapply(rejections,function(x)
    sum(h[x])/total_nonnull))
  return(tpr)

}

calc_fdr <- function(rejections, h,alphas){
  if(length(rejections) == 0){
    return(rep(0,length(alphas)))
  }

  fdr <- unlist(lapply(rejections,function(x)
    sum(1-h[x])/max(1,length(x))))
  return(fdr)
}


run_method <- function(x,z,pvals,se,intervals,method,alphas,h,testing){

  tau <- 0.5; eps <- 0.1 # parameters for SABHA
  thr <- 0.5 # parameter for Storey-BH
  thr1 <- 0.1; thr2 <- 0.5 # parameters for adaptive SeqStep
  ADMM_params <- c(10^2, 10^3, 2, 5000, 1e-3) # alpha_ADMM,beta,eta,max_iters,converge_thr
  n <- length(pvals)

  num_alpha <- length(alphas)
  max_alpha <- max(alphas)

  formulas <- "splines::ns(x, df = 5)"
  full_formulas <- paste("splines::ns(x, df = ",c(2,4,6,8)," )")
  gam_formulas <- c("s(x,k=2)","s(x,k=3)")
  if(is.data.frame(x)){
    if(ncol(x)>1){
      formulas <- paste(colnames(x),collapse=" + ")

    }
  }
  out <- list(rejs=rep(list(),length(alphas)),qvals=list())

  if(method=="IHW"){
    # if (!require(IHW)) install.packages('IHW')
    output <-  ihw(pvalues=pvals,covariates=x,alpha=0.1)
    for(ind in 1:length(alphas)){
      alpha <- alphas[ind]
      out$rejs[[ind]] <- which(adj_pvalues(output) <= alpha)
    }

    out$qvals <- adj_pvalues(output)

  }else if(method=="IHW Oracle"){
    output <-  ihw_oracle(pvals,x,alpha)
    for(ind in 1:length(alphas)){
      alpha <- alphas[ind]
      out$rejs[[ind]] <- which(adj_pvalues(output) <= alpha)
    }
    out$qvals <- adj_pvalues(output)
  }else if(method=="BH"){
    for(ind in 1:length(alphas)){
      alpha <- alphas[ind]
      output <-  BH_method(pvals, alpha)
      out$rejs[[ind]] <- output
    }
    out$qvals <- p.adjust(pvals,method="BH")
  }else if(method=="Storey BH"){
    for(ind in 1:length(alphas)){
      alpha <- alphas[ind]
      output <-  Storey_method(pvals, alpha,thr)
      out$rejs[[ind]] <- output
    }
  }else if(method=="Barber Candes"){
    for(ind in 1:length(alphas)){
      alpha <- alphas[ind]
      output <-  BC_method(pvals, alpha)
      out$rejs[[ind]] <- output
    }
  }else if(method=="Boca Leek"){
    qvals <- BL_method(pvals,x)
    for(ind in 1:length(alphas)){
      alpha <- alphas[ind]
      out$rejs[[ind]] <- which(qvals<alpha)
    }
    out$qvals <- qvals

  }else if(method=="SABHA (step)" | method == "SABHA (ordered)"){

    ## gather results
    NumRej <- matrix(0, nrow = 13, ncol = num_alpha)

    ## methods: 1 SeqStep, 2 HingeExp, 3 ForwardStop, 4 Adaptive SeqStep, 5 BH, 6 Storey-BH, 7 Barber-Candes, 8 SABHA (step), 9 SABHA (ordered), 10 IHW, 11 IHW (oracle), 12 IF (oracle), 13 AdaPT

    if(method =="SABHA (step)"){
      qhat_step <- Solve_q_step(pvals, tau, eps)
      for(ind in 1:length(alphas)){
        alpha <- alphas[ind]
        out$rejs[[ind]] <- SABHA_method(pvals, qhat_step, alpha, tau)
      }
    }else{
      qhat_ordered <- Solve_q_ordered_simple(pvals, tau, eps, ADMM_params)
      for(ind in 1:length(alphas)){
        alpha <- alphas[ind]
        out$rejs[[ind]] <- SABHA_method(pvals, qhat_ordered, alpha, tau)
      }
    }
  }else if(method == "LFDR"){
    #if (!require(fdrtool)) install.packages('fdrtool')

    qvals <- lfdr_method(pvals,x)
    for(ind in 1:length(alphas)){
      alpha <- alphas[ind]
      out$rejs[[ind]] <- which(qvals<alpha)
    }
    out$qvals <- qvals
  }else if(method == "FDRreg-t"){
    #if (!require(fdrtool)) install.packages('fdrtool')

    qvals <- fdrreg_method(z,x,nulltype="theoretical")
    for(ind in 1:length(alphas)){
      alpha <- alphas[ind]
      out$rejs[[ind]] <- which(qvals<alpha)
    }
    out$qvals <- qvals
  }else if(method == "FDRreg-e"){
    #if (!require(fdrtool)) install.packages('fdrtool')

    qvals <- fdrreg_method(z,x,nulltype="empirical")
    for(ind in 1:length(alphas)){
      alpha <- alphas[ind]
      out$rejs[[ind]] <- which(qvals<alpha)
    }
    out$qvals <- qvals
  }else if(method== "ASH"){
    # if (!require(ashr)) install.packages('ashr')
    if(testing != "interval"){


      qvals <- ash_method(z,se)
      for(ind in 1:length(alphas)){
        alpha <- alphas[ind]
        out$rejs[[ind]] <- which(qvals<= alpha)
      }
      out$qvals <- qvals
    }


  }else if(method == "AdaFDR"){
    for(ind in 1:length(alphas)){
      alpha <- alphas[ind]
      print(alpha)
      rejs <-  adafdr_test(pvals, as.matrix(x),alpha=alpha,fast_mode = TRUE)$decision
      out$rejs[[ind]] <- which(rejs)# which(qvals<= alpha)

    }
  }else if(method=="AdaPTGMM"){

    #devtools::install_github("patrickrchao/AdaPTGMM")
    if(!is.data.frame(x)){
      x <- data.frame(x=x)
    }
    if(testing == "interval"){
      out <- adapt_gmm(x=x,z=z,se=se,model_type="neural",intercept_model = F,
                       testing="interval",lendpoint=intervals[1],rendpoint = intervals[2],
                       beta_formulas=formulas,alphas=alphas,nclasses=c(5),nfits = 5,niter_fit = 3,
                       masking_shape = "comb")
      out_zeta <- adapt_gmm(x=x,z=z,se=se,model_type="neural",intercept_model = F,
                       testing="interval",lendpoint=intervals[1],rendpoint = intervals[2],
                       beta_formulas=formulas,alphas=alphas,nclasses=c(5),nfits = 5,niter_fit = 3,
                       masking_shape = "comb",target_alpha_level = 0.01)
    }else{
      out <- adapt_gmm(x=x,z=z,se=se,model_type="neural",intercept_model = F,testing=testing,
                       beta_formulas=formulas,alphas=alphas,nclasses=c(5),nfits = 5,niter_fit = 3)
      out_zeta <- adapt_gmm(x=x,z=z,se=se,model_type="neural",intercept_model = F,testing=testing,
                       beta_formulas=formulas,alphas=alphas,nclasses=c(5),nfits = 5,niter_fit = 3,
                       target_alpha_level=0.01)
    }
    out$rejs[1:5] <- out_zeta$rejs[1:5]
    print(out$args)
  }else if(method=="AdaPT"){
    #if (!require(adaptMT)) install_github("patrickrchao/adaptMT")
    if(!is.data.frame(x)){
      x <- data.frame(x=x)
    }
    dist <- beta_family()
    out <- adapt_glm(x = x, pvals = pvals, pi_formulas = formulas,
                     mu_formulas = formulas, dist = dist, alphas=alphas,
                     cr="AIC",
                     alpha_m=0.5,lambda=0.5,zeta=1,
    )
    print(out$masking_params)

  }else if(method=="AdaPTg"){
    #if (!require(adaptMT)) install_github("patrickrchao/adaptMT")
    dist <- beta_family()
    if(!is.data.frame(x)){
      x <- data.frame(x=x)
    }
    if( testing == "interval"){
      shape <-  "comb"
    }else{
      shape <-  "tent"
    }
    out <- adapt_glm(x = x, pvals = pvals, pi_formulas = formulas,
                     mu_formulas = formulas, dist = dist, alphas=alphas,
                     cr="AIC",masking_shape=shape)
    out_zeta <- adapt_glm(x = x, pvals = pvals, pi_formulas = formulas,
                     mu_formulas = formulas, dist = dist, alphas=alphas,
                     cr="AIC",masking_shape=shape,alpha_m=0.9/11,lambda=0.9/11,zeta=10)
    out$rejs[1:5] <- out_zeta$rejs[1:5]
    print(out$masking_params)

  }else{
    print(paste(method,"method not found."))
    return(0)
  }
  return(out)


}



create_df <- function(rejections,h,method,alphas){

  if(length(rejections)==0){
    nrejs <- rep(0,length(alphas))
    tpr <- rep(0,length(alphas))
    fdr <- rep(0,length(alphas))
  }else{
    nrejs <- sapply(rejections,length)
    tpr <- calc_tpr(rejections,h,alphas)
    fdr <- calc_fdr(rejections,h,alphas)
  }

  df <- data.frame(method,alphas,tpr,fdr,nrejs)
  colnames(df) <- c("Method","Alpha","TPR","FDR","Rej")


  return(df)
}

#
# stitch <- function(zeta_df,orig_df,title){
#   orig_df<- orig_df[,!(colnames(orig_df) %in% c("FDR_diff"))]
#   zeta_df <- zeta_df[,!(colnames(zeta_df) %in% c("FDR_diff"))]
#   output <- orig_df
#   #output <- rbind(zeta_df[zeta_df$Alpha==0.01,],output)
#   output[output$Method=="AdaPTGMM" & output$Alpha<=0.05, ] <-  zeta_df [zeta_df$Method=="AdaPTGMM" & zeta_df$Alpha<=0.05, ]
#   output[output$Method=="AdaPTg" & output$Alpha<=0.04, ] <-  zeta_df [zeta_df$Method=="AdaPTg" & zeta_df$Alpha<=0.04, ]
#   output[output$Method=="AdaPTGMM" & output$Alpha>0.05, ]$TPR <-  rowMeans(cbind(zeta_df [zeta_df$Method=="AdaPTGMM" & zeta_df$Alpha>0.05, ]$TPR,
#                                                                                  orig_df [orig_df$Method=="AdaPTGMM" & orig_df$Alpha>0.05, ]$TPR))
#   output[output$Method=="AdaPTg" & output$Alpha>0.04, ]$TPR <-  rowMeans(cbind(zeta_df [zeta_df$Method=="AdaPTg" & zeta_df$Alpha>0.04, ]$TPR,
#                                                                                orig_df [orig_df$Method=="AdaPTg" & orig_df$Alpha>0.04, ]$TPR))
#   saveRDS(output,paste0("data/combined_exp4/",title,".RDS"))
# }
#
#
# stitch2 <- function(adapt_df,orig_df,title){
#   orig_df<- orig_df[,!(colnames(orig_df) %in% c("FDR_diff"))]
#   adapt_df <- adapt_df[,!(colnames(adapt_df) %in% c("FDR_diff"))]
#   output <- orig_df
#   output <- output[output$Method != "AdaPTGMM" & output$Method!="AdaPTg",]
#   #output <- rbind(zeta_df[zeta_df$Alpha==0.01,],output)
#   #output[output$Method=="AdaPTGMM" & output$Alpha<=0.04, ] <-  zeta_df [zeta_df$Method=="AdaPTGMM" & zeta_df$Alpha<=0.04, ]
#   output <- rbind(output,adapt_df)
#   saveRDS(output,paste0("data/combined_exp4_full/",title,".RDS"))
# }
#
# # stitch(readRDS("data/exp(-4x)/zeta=10/Logistic Interval.RDS"),readRDS("data/exp(-4x)/zeta=2/Logistic Interval.RDS"),"Logistic Interval")
# # stitch(readRDS("data/exp(-4x)/zeta=10/Logistic One Sided.RDS"),readRDS("data/exp(-4x)/zeta=2/Logistic One Sided.RDS"),"Logistic One Sided")
# # stitch(readRDS("data/exp(-4x)/zeta=10/Logistic Two Sided.RDS"),readRDS("data/exp(-4x)/zeta=2/Logistic Two Sided.RDS"),"Logistic Two Sided")
# #
# # stitch2(readRDS("data/combined_exp4/Logistic One Sided.RDS"),readRDS("data/exp(-4x)/Logistic One Sided.RDS"),"Logistic One Sided")
# # stitch2(readRDS("data/combined_exp4/Logistic Two Sided.RDS"),readRDS("data/exp(-4x)/Logistic Two Sided.RDS"),"Logistic Two Sided")
# # stitch2(readRDS("data/combined_exp4/Logistic Interval.RDS"),readRDS("data/exp(-4x)/Logistic Interval.RDS"),"Logistic Interval")
#
# # stitch(readRDS("data/zeta=10/High_Dimensional Two Sided.RDS"),readRDS("data/final/High_Dimensional Two Sided.RDS"),"Logistic Two Sided")
# # stitch(readRDS("data/zeta=10/High_Dimensional Interval.RDS"),readRDS("data/final/High_Dimensional Interval.RDS"),"Logistic Interval")
# # stitch(readRDS("data/zeta=10/High_Dimensional One Sided.RDS"),readRDS("data/final/High_Dimensional One Sided.RDS"),"Logistic One Sided")
#
