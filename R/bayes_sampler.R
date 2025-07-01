#########################################
# bayes_sampler.R                       #
# Author:  Caleb Frankenberger          #
# Purpose: Runs Bayesian Gibbs sampler  #
#########################################


bayes_sampler <- function(beta_init, Z, X, N, S, g, a, R, post_git, known_acc, se_sp){
  X <- as.matrix(X)
  Yt <- rbinom(N, 1, apply(X %*% beta_init, 2, g))
  
  # Sort and extract assay info
  Z <- Z[order(Z[,5]), ]
  assay_vec <- Z[, 5]
  unique_assays <- sort(unique(assay_vec))
  num_assays <- length(unique_assays)
  
  # Construct Y-til matrix
  tmp <- Z[, -(1:5)]
  ytm <- matrix(-9L, N, S)
  for (d in 1:N) {
    idx <- which(tmp == d, arr.ind = TRUE)[, "row"]
    if (length(idx) > 0) {
      ytm[d, 1:length(idx)] <- sort(idx)
    }
  }
  Ytmat <- cbind(Yt, rowSums(ytm > 0), ytm)
  Ycol <- ncol(Ytmat)
  
  # Store necessary Z data 
  Z <- Z[ ,-(3:5)]
  Zrow <- nrow(Z)
  Zcol <- ncol(Z)
  
  accept <- rep(-9, post_git)
  beta_sv <- matrix(-9, post_git, length(beta_init))
  
  if (!known_acc) {
    se_sv <- matrix(NA, post_git, num_assays)
    sp_sv <- matrix(NA, post_git, num_assays)
  }
  
  # Begin Gibbs sampling:
  for (s in 1:post_git) {
    pvec <- g(drop(X %*% beta_init))
    U <- runif(N)
    
    # Sampling individual true statuses
    Ztil <- integer(Zrow * Zcol)
    res <- .C("sample_latent_statuses",
              as.double(pvec),
              as.integer(Ytmat),
              as.integer(Z),
              as.integer(N),
              as.double(se_sp),
              as.integer(Ycol),
              as.integer(Zrow),
              as.integer(Zcol),
              as.double(U),
              Ztil = as.integer(Ztil))
    
    Ytmat <- matrix(res[[2]], N, Ycol)
    Ztil <- matrix(res$Ztil, Zrow, Zcol) 
    
    # Sampling beta
    out <- wls.mh.alg(b0 = beta_init, X = X, y = Ytmat[, 1], a = a, R = R)
    beta_init <- out$param
    beta_sv[s, ] <- beta_init
    
    # Sampling Se & Sp:
    if (!known_acc) {
      Ztest <- Z[, 1]         # Observed pool result
      Ztrue <- Ztil[, 1]      # True pool result
      
      se <- sp <- numeric(num_assays)
      a_se0 <- b_se0 <- a_sp0 <- b_sp0 <- 1 
      
      for (asy in seq_len(num_assays)) {
        idx <- assay_vec == unique_assays[asy]
        
        Zj <- Ztest[idx]
        Ztilj <- Ztrue[idx]
        
        a_se_star <- a_se0 + sum(Zj * Ztilj)
        b_se_star <- b_se0 + sum((1 - Zj) * Ztilj)
        
        a_sp_star <- a_sp0 + sum((1 - Zj) * (1 - Ztilj))
        b_sp_star <- b_sp0 + sum(Zj * (1 - Ztilj))
        
        se[asy] <- rbeta(1, a_se_star, b_se_star)
        sp[asy] <- rbeta(1, a_sp_star, b_sp_star)
      }
      
      se_sv[s, ] <- se
      sp_sv[s, ] <- sp
      
      se_sp <- cbind(
        se[match(assay_vec, unique_assays)],
        sp[match(assay_vec, unique_assays)]
      )
    }
    
    
    # Track acceptance in the MH algorithm:
    accept[s] <- out$accept
  }
  
  output <- list(param = beta_sv, convergence = 0, accept = accept)
  if (!known_acc) {
    colnames(se_sv) <- paste0("Se_a", unique_assays)
    colnames(sp_sv) <- paste0("Sp_a", unique_assays)
    output$se <- se_sv
    output$sp <- sp_sv
  }
  return(output)
}


## END OF FILE ##