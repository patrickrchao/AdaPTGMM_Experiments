run_simulations <- function(chosen_methods, num_sim, num_hypo, radius, testing, se, num_cores,alphas=seq(0.01,0.2,by=0.01)){

  full_log <- data.frame(matrix(ncol=6,nrow=0))
  colnames(full_log) <- c("Method","Alpha","TPR","FDR","Rej","Time")



  if(testing != "interval"){
    chosen_methods <- chosen_methods[chosen_methods != "AdaPTGMM Interval"]
  }else{
    intervals <- c(-radius,radius)
  }

  method_subset <- chosen_methods[!(chosen_methods%in% c("AdaFDR"))]

  full_log <- simulation_helper(num_hypo,num_sim,radius,testing,se,alphas,method_subset,num_cores,intervals)

  print("Finished methods, now running AdaFDR")
  saveRDS(full_log,paste0("logs/",gsub(" ","_",testing),"_log.RDS"))
  if("AdaFDR" %in% chosen_methods){
    #AdaFDR does not play nicely with parallelization
    full_log <- rbind(full_log,simulation_helper(num_hypo,num_sim,radius,testing,se,alphas,c("AdaFDR"),1))
  }

  saveRDS(full_log,paste0("logs/",gsub(" ","_",testing),"_log.RDS"))
  return(full_log)
}

simulation_helper <- function(num_hypo,num_sim,radius,testing,se,alphas,methods,num_cores,intervals){

  log <- mclapply(1:num_sim,function(x){
    print(paste("Iteration",x,"of",num_sim))
    full_log <- data.frame(matrix(ncol=6,nrow=0))
    data <- logistic_data(num_hypo,radius=radius,testing=testing,se=se)

    x <- data$x
    pvals <- data$pvals
    h <- data$h
    z <- data$z

    for(method in methods){
        print(paste("Current method:",method))
        out <- run_method(x=x,z=z,pvals=pvals,se=se,intervals=intervals,
                          method=method,alphas=alphas,h=h,testing=testing)

        rejs <- out$rejs

        full_log <- rbind(full_log,create_df(rejs,h,method,alphas))
      }

    return(full_log)
  },mc.cores=num_cores)

  full_log <- bind_rows(log)
  return(full_log)
}
