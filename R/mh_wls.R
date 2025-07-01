#MH_WLS.R

## DISCLAIMER: THIS CODE WAS NOT WRITTEN BY ME, AND IS CURRENTLY HERE FOR TESTING PURPOSES.
## IT WILL EVENTUALLY BE REPLACED WITH MORE GENERAL CODE.

########################################################
# Bayesian Regression Models (McMahan et al., 2017)   
# Assay Accuracy Parameters (Se & Sp): KNOWN          
# This code was written in 2020                       
########################################################


############################################
# This weighted least square (WLS) Metropolis-
# Hastings (MH) algorithm works only with 
# "logit" link function.
#
wls.mh.alg <- function(b0,X,y,a,R){
  library(mnormt)
  
  res <- mt.ct(b=b0,y=y,X=X,sig=0,a=a,R=R)
  Mb0 <- res$M
  Cb0 <- res$C
  Lb0 <- res$L
  sig <- res$s
  
  b.star <- as.vector(rmnorm(1,Mb0,Cb0))
  qb.star <- dmnorm(b.star,as.vector(Mb0),Cb0)
  
  res2 <- mt.ct(b=b.star,y=y,X=X,sig=sig,a=a,R=R) 
  Mb.star <- res2$M
  Cb.star <- res2$C
  Lb.star <- res2$L
  
  qb0 <- dmnorm(b0,as.vector(Mb.star),Cb.star)
 
  ## THIS SECTION HAS BEEN CHANGED TO 
  ## FIX CRASHING CAUSED BY RBINOM FAILING
  ratio.pi <- exp(Lb.star - Lb0)
  raw_alpha <- (ratio.pi * qb0) / qb.star
  
  if (!is.finite(raw_alpha)) {
    alpha <- 0
  } else {
    alpha <- min(1, max(0, raw_alpha))
  }
  
  accept <- rbinom(1, 1, alpha)
  #########################################
  
  ##ORIGINAL:
  # ratio.pi <- exp(Lb.star - Lb0)
  # alpha <- min(1, max(0, ((ratio.pi * qb0) / (qb.star))))
  # accept <- rbinom(1, 1, alpha)
  
  if(accept==1){b1 <- b.star}
  if(accept==0){b1 <- b0}
  return(list("param"=b1,"accept"=accept))
}

############################################
# Supporting program for the MH sampling 
# function 'wls.mh.alg()' based on WLS method.
# Works only with logit link fn.
#
mt.ct <- function(b,y,X,sig,a,R){
  P <- length(b)
  xb <- X%*%b
  
  p <- plogis(xb)
  W <- p*(1-p)
  Yb <- xb*(p*(1-p))+(y-p)
  XWX <- t(X*matrix(W,ncol=P,nrow=length(W),byrow=FALSE))%*%X
  diag(XWX) <- diag(XWX) + 1e-6
  
  log_den <- log1p(exp(-abs(xb))) + pmax(xb, 0)
  
  if(is.null(R)){
    C <- solve(XWX)
    C <- (C+t(C))/2
    M <- C%*%(t(X)%*%(Yb))
    L <- sum(y * xb - log_den)
  }else{
    R.inv <- solve(R)
    C <- solve(R.inv+XWX)
    C <- (C+t(C))/2
    M <- C%*%(R.inv%*%a+t(X)%*%(Yb))
    L <- (-1/2)*((b-a)%*%R.inv%*%(b-a))+sum(y*(xb)-log_den)
  }
  return(list("M"=M, "C"=C, "L"=L, "s"=sig))
}
