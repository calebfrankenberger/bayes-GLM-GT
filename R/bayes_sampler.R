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
  
  Ztil_out <- integer(Zrow * Zcol)
  se_out   <- integer(2L * num_assays)
  sp_out   <- integer(2L * num_assays)
  
  
  accept <- rep(-9, post_git)
  beta_sv <- matrix(-9, post_git, length(beta_init))
  
  if (!known_acc) {
    se_sv <- matrix(NA, post_git, num_assays)
    sp_sv <- matrix(NA, post_git, num_assays)
    se <- rep(settings$se_0, length.out = num_assays)
    sp <- rep(settings$sp_0, length.out = num_assays)
  }
  
  U_all <- matrix(runif(N * post_git), nrow = N, ncol = post_git)
  assay_map <- match(assay_vec, unique_assays)
  se_sp <- matrix(NA_real_, nrow = length(assay_vec), ncol = 2)
  
  # Begin Gibbs sampling:
  for (s in 1:post_git) {
    pvec <- g(drop(X %*% beta_init))
    res <- .C("sample",
              p            = as.double(pvec),
              Ytmat        = as.integer(Ytmat),
              Z_mat        = as.integer(Z),
              N            = as.integer(N),
              Ycols        = as.integer(Ycol),
              Zrows        = as.integer(Zrow),
              Zcols        = as.integer(Zcol),
              U            = as.double(U_all[, s]),
              Ztil         = as.integer(Ztil_out),
              assay_vec    = as.integer(assay_vec),
              L            = as.integer(num_assays),
              se_in        = as.double(se),
              sp_in        = as.double(sp),
              se_counts    = as.integer(se_out),
              sp_counts    = as.integer(sp_out)
    )
    Ytmat     <- matrix(res$Ytmat, N, Ycol)
    Ztil      <- matrix(res$Ztil, Zrow, Zcol)
    se_counts <- matrix(res$se_counts, nrow = 2)
    sp_counts <- matrix(res$sp_counts, nrow = 2)
    
    # Sampling beta
    out <- wls.mh.alg(b0 = beta_init, X = X, y = Ytmat[, 1], a = a, R = R)
    beta_init <- out$param
    beta_sv[s, ] <- beta_init
    
    # Sampling Se & Sp:
    if (!known_acc) {
      a_se0 <- b_se0 <- a_sp0 <- b_sp0 <- 1
      se <- sp <- numeric(num_assays)
      
      for (asy in seq_len(num_assays)) {
        tp <- se_counts[1, asy]
        fn <- se_counts[2, asy]
        tn <- sp_counts[1, asy]
        fp <- sp_counts[2, asy]
        
        se[asy] <- rbeta(1, a_se0 + tp, b_se0 + fn)
        sp[asy] <- rbeta(1, a_sp0 + tn, b_sp0 + fp)
      }
      
      se_sv[s, ] <- se 
      sp_sv[s, ] <- sp
      se_sp[, 1] <- se[assay_map]
      se_sp[, 2] <- sp[assay_map]
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