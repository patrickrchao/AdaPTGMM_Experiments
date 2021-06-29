#' @author Patrick Chao
combine_three_tidy_df <- function(objects,objects2,objects3, colLabels, fill, annotate, alpha,methodset){
  ## create tidy data frame where each row is a method / dataset observation
  ## of a rank 
  ranks <- data.frame()
  
  for (i in seq_along(objects)){
    if ( class(objects[[i]]) == "character") {
      x <- readRDS(objects[i])
      y <- readRDS(objects2[i])
      z <- readRDS(objects3[i])
    }else if ( class(objects[[i]]) == "SummarizedBenchmark") {
      x <- objects[[i]]
      y <- objects2[[i]]
      z <- objects3[[i]]
    }
    assayNames(x) <- "qvalue"
    x <- addDefaultMetrics( x )
    
    assayNames(y) <- "qvalue"
    y <- addDefaultMetrics( y )
    
    assayNames(z) <- "qvalue"
    z <- addDefaultMetrics( z )
    
    hasResults <- apply(!is.na( assays(x)[["qvalue"]] ), 2, sum)
    NAmethods <- names(hasResults)[hasResults == 0]
    
    
    if (fill %in% c("propMaxRejections", "meanRank")){
      pmcol <- "rejections"
    }else{
      pmcol <- fill
    }
    
    if (sum(hasResults) > 0){
      tmp <- data.frame()
      tmp2 <- data.frame()
      tmp3 <- data.frame()
      for (lev in alpha){
        df <- estimatePerformanceMetrics(x, alpha=lev)
        df$alpha <- lev
        tmp <- bind_rows(tmp, zeroRejectionsFDR_new(df, tidy=TRUE)) 
        
        df2 <- estimatePerformanceMetrics(y, alpha=lev)
        df2$alpha <- lev
        tmp2 <- bind_rows(tmp2, zeroRejectionsFDR_new(df2, tidy=TRUE)) 
        
        df3 <- estimatePerformanceMetrics(z, alpha=lev)
        df3$alpha <- lev
        tmp3 <- bind_rows(tmp3, zeroRejectionsFDR_new(df3, tidy=TRUE)) 
      }
      
      
      if (fill %in% c("TPR", "FDR", "TNR") && !(fill %in% tmp$performanceMetric))
        stop(fill, " is not found in performanceMetrics")
      
      if (!is.null(annotate)){
        if (annotate %in% c("TPR", "FDR", "TNR") &&
            !(annotate %in% tmp$performanceMetric))
          stop(annotate, " is not found in performanceMetrics")
      }
      
      tmp2 <- tmp2 %>% dplyr::rename(blabel = func.pkg)
      
      tmp3 <- tmp3 %>% dplyr::rename(blabel = func.pkg)
      
      new_names <- c("New1","New2")
      
      tmp2$blabel <- new_names[1]
      tmp3$blabel <- new_names[2]
      tmp <- bind_rows(tmp,tmp2,tmp3)
      tmp <- tmp %>%
        dplyr::filter( performanceMetric == pmcol) %>%
        dplyr::rename( method = blabel) %>%
        dplyr::filter( method %in% c(new_names,methodset))  %>%
        dplyr::filter( is.na(param.alpha) | (param.alpha == alpha)) %>%
        dplyr::filter( is.na(param.smooth.df) | (param.smooth.df == "3L")) %>%
        dplyr::filter( !method == "unadjusted") %>%
        dplyr::filter( !(method %in% NAmethods)) %>%
        dplyr::select( method, value, alpha ) %>%
        dplyr::rename( nrejects = value) %>%
        mutate( method = gsub("-df03", "", method)) %>%
        mutate( method = gsub("(-a)(.*)", "", method)) %>%
        mutate( rank = rank(nrejects) / n(),
                propMaxRej = ifelse(nrejects == min(nrejects) & 
                                      nrejects == max(nrejects) &
                                      nrejects == 0, 0, 
                                    nrejects / max(nrejects, na.rm = TRUE)),
                propPossible = nrejects / nrow(x)) %>%
        mutate(casestudy = colLabels[i])
      if ( is.list(objects)) {
        tmp <- tmp %>%
          mutate(replicate = i)
      }
      
      ranks <- rbind(ranks, tmp)
    }
  }
  return(ranks)
}


# function to replace observations of NA FDR when nrejects = 0 to 0
# useful when averaging over simulation reps so we don't penalize a method
# that almost never rejects anything but we take the mean over all the 
# times it does
#' @param tidy_df a data frame output from estimatePerformanceMetrics with the 
#' tidy = FALSE option, where 
#' each row is an observation of a performance metric for a method. Assumes that
#' method is in the blabel column, the performance metric is in the performanceMetric
#' column, and that the value of the performance metric is in the value column.
#' @param tidy logical, whether to return the metrics table in tidy format or not
#' (analagous to the same argument in tidyUpMetrics function). Default is TRUE.
#' @author Keegan Korthauer, edited by Patrick Chao
zeroRejectionsFDR_new <- function(df, tidy = TRUE){
  df$FDR <- ifelse(df$rejections == 0, 0, df$FDR)
  
  if(tidy){
    valueCols <- c("TPR", "FDR", "TNR", "FNR", "rejections")
    # PATRICK CHAO EDIT 
    #tidy_df <- gather(as.data.frame(df), keys = valueCols) %>%
    #dplyr::rename(performanceMetric = key)
    
    tidy_df <- gather(as.data.frame(df), key,value, valueCols) %>%  
      dplyr::rename(performanceMetric = key)
    #dplyr::rename(blabel = value)
    return(tidy_df)
  }else{
    return(df)
  }
}



#' Function to create a tidy data frame of performance metrics 
#' 
#' Takes as input a vector of filepaths, each one pointing to a SummarizedBenchmark
#' results object, labels for each object, as well as various annotation variables.
#' Used internally by plotMethodRanks, but can also be used on its own.
#' 
#' @param objects a character vector containing the full filepaths to 
#'  SummarizedBenchmark results objects.
#' @param colLabels a character vector containing the labels to use for the 
#'  heatmap columns (same length and order as `objects`). 
#' @param alpha numeric value  or vector indicating which alpha value(s) to use for IHW 
#'  (need to be among of the values specified in the bench object). Default 0.10.
#' @param fill character indicating the outcome variable of interest. Default
#'  is 'propMaxRejections'. For 'FDR' and 'TPR' sb objects need to contain 
#'  these performance metrics.
#' @param annotate character indicating what should be plotted as text labels
#'  (heatmap plot only). Default is proportion of total possible rejections 
#'  (propPossible). Could also be another variable (from \code{fill} choices) 
#'  or NULL (for no text label annotations).
#' @author Keegan Korthauer, edited by Patrick Chao 
new_tidy_df <- function(objects,colLabels, fill, annotate, alpha){
  ## create tidy data frame where each row is a method / dataset observation
  ## of a rank 
  ranks <- data.frame()

    

  for (i in seq_along(objects)){
    
    if ( is(objects[[i]], "character")) {
      x <- readRDS(objects[i])
    }else if ( is(objects[[i]], "SummarizedBenchmark")) {
      x <- objects[[i]]
    }
    
    assayNames(x) <- "qvalue"
    paste("file",i,"out of",length(objects))
    
    x <- addDefaultMetrics( x )
    
    hasResults <- apply(!is.na( assays(x)[["qvalue"]] ), 2, sum)
    NAmethods <- names(hasResults)[hasResults == 0]
    
    if (fill %in% c("propMaxRejections", "meanRank")){
      pmcol <- "rejections"
    }else{
      pmcol <- fill
    }
    
    if (sum(hasResults) > 0){
      tmp <- data.frame()
      
      for (lev in alpha){
        df <- estimatePerformanceMetrics(x, alpha=lev)
        df$alpha <- lev
        tmp <- bind_rows(tmp, zeroRejectionsFDR_new(df, tidy=TRUE))  
      }
      
      if (fill %in% c("TPR", "FDR", "TNR") && !(fill %in% tmp$performanceMetric))
        stop(fill, " is not found in performanceMetrics")
      
      if (!is.null(annotate)){
        if (annotate %in% c("TPR", "FDR", "TNR") &&
            !(annotate %in% tmp$performanceMetric))
          stop(annotate, " is not found in performanceMetrics")
      }
      if( !("blabel" %in% colnames(tmp))){
        tmp <- tmp %>% dplyr::rename(blabel = func.pkg)
      }
      
      tmp <- tmp %>%
        dplyr::filter( performanceMetric == pmcol) %>%
        dplyr::rename( method = blabel) 
      if("param.alpha" %in% colnames(tmp) &"param.smooth.df" %in% colnames(tmp)  ){
        tmp <- tmp%>% dplyr::filter( is.na(param.alpha) | (param.alpha == alpha)) %>%
            dplyr::filter( is.na(param.smooth.df) | (param.smooth.df == "3L")) 
      }
      tmp <- tmp %>%
        dplyr::filter( !method == "unadjusted") %>%
        dplyr::filter( !(method %in% NAmethods)) %>%
        dplyr::select( method, value, alpha ) %>%
        dplyr::rename( nrejects = value) %>%
        mutate( method = gsub("-df03", "", method)) %>%
        mutate( method = gsub("(-a)(.*)", "", method)) %>%
        mutate( rank = rank(nrejects) / n(),
                propMaxRej = ifelse(nrejects == min(nrejects) & 
                                      nrejects == max(nrejects) &
                                      nrejects == 0, 0, 
                                    nrejects / max(nrejects, na.rm = TRUE)),
                propPossible = nrejects / nrow(x)) %>%
        mutate(casestudy = colLabels[i])
      if ( is.list(objects)) {
        tmp <- tmp %>%
          mutate(replicate = i)
      }
      
      ranks <- rbind(ranks, tmp)
    }
  }
 
  return(ranks)
}

