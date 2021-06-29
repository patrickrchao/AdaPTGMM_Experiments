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
#
# logistic_data <- function(n=100, radius,dim = 20, testing="interval",se=se, x=NULL){
#   if(is.null(x)){
#     #x <- matrix(rnorm(n*dim),ncol=dim)
#     x <- rnorm(n)
#   }
#   x <- sort(x)
#   mu <- rep(0,n)
#   #beta <- c(rep(2,10),rep(0,dim-10))
#   subset <- (runif(n) < 0.75/ (1+ exp(- 6*(x - 1.5 ))))
#   n1 <- sum(subset)
#
#   mean <- sample(c(2),n1,replace=T)
#   mu[subset] <- rlogis(n1,location = mean, scale = 1/2)
#
#
#   if(testing == "interval"){
#     z <- rnorm(n,mean=mu)
#     pvals <- 1-pnorm(abs(z)+radius) + pnorm(-abs(z)+radius)
#     h <- abs(mu) > radius
#   }else if(testing == "two_sided"){
#
#     z <- rnorm(n,mean=mu,sd=se)
#     pvals <- 2*(pnorm(abs(z),sd=se,lower.tail = FALSE))
#     h <- mu != 0
#   }else{
#     z <- rnorm(n,mean=mu,sd=se)
#     pvals <- pnorm(z,sd=se,lower.tail = FALSE)
#     h <- mu > 0
#   }
#   #x <- as.data.frame(x)
#   #colnames(x) <- paste0("x",1:dim)
#   out <- list()
#   #out$x <- x
#   out$x <- x #x%*% beta
#   out$pvals <- pvals
#   out$h <- h
#   out$z <- z
#
#
#   return(out)
# }
#
#
#
#
#
