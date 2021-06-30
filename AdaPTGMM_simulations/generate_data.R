logistic_data <- function(n=100, radius,testing="interval",se=se, x=NULL){
  if(is.null(x)){
    x <- rnorm(n)
  }
  x <- sort(x)
  H <- (runif(n) < 0.75 / (1+exp(-(6*x-9))))
  mu <- H * rlogis(n,location = 2, scale = 1/2)


  if(testing == "interval"){
    z <- rnorm(n,mean=mu,sd=se)
    pvals <- 1-pnorm(abs(z)+radius,sd=se) + pnorm(-abs(z)+radius,sd=se)
    h <- abs(mu) > radius
  }else if(testing == "two_sided"){

    z <- rnorm(n,mean=mu,sd=se)
    pvals <- 2*(pnorm(abs(z),sd=se,lower.tail = FALSE))
    h <- mu != 0
  }else{
    z <- rnorm(n,mean=mu,sd=se)
    pvals <- pnorm(z,sd=se,lower.tail = FALSE)
    h <- mu > 0
  }

  out <- data.frame(x=x,pvals=pvals,h=h,z=z,row.names = NULL)
  return(out)
}
