################################
# estimator.R                  #
# Author: Caleb Frankenberger  #
################################

library(parallel)

#########################
## INFER POSTERIOR     ##
#########################
infer_posterior <- function(test_data, settings) {
  X <- test_data$X
  Z <- test_data$Z
  tsts <- test_data$tsts
  
  N <- nrow(X)
  S <- length(settings$psz)
  
  Z <- Z[order(Z[, 5]), ]
  assay_vec <- Z[, 5]
  unique_assays <- sort(unique(assay_vec))
  
  if (settings$known_acc) {
    se_sp <- Z[, 3:4]
  } else {
    n_assays <- length(unique(assay_vec))
    
    se_0 <- settings$se_0
    sp_0 <- settings$sp_0
    
    if (length(se_0) == 1) se_0 <- rep(se_0, n_assays)
    if (length(sp_0) == 1) sp_0 <- rep(sp_0, n_assays)
    
    stopifnot(length(se_0) == n_assays, length(sp_0) == n_assays)
    
    se_sp_init <- cbind(se_0, sp_0)
    
    se_sp <- cbind(
      se_sp_init[match(assay_vec, unique(assay_vec)), 1],
      se_sp_init[match(assay_vec, unique(assay_vec)), 2]
    )
  }
  
  beta_init <- rep(0, length(settings$beta_true))
  res <- bayes_sampler(
    beta_init = beta_init,
    Z         = Z,
    X         = X,
    N         = N,
    S         = S,
    g         = settings$g,
    a         = NULL,
    R         = NULL,
    post_git   = settings$post_git,
    known_acc = settings$known_acc,
    se_sp     = se_sp
  )
  
  all_samp   <- res$param
  post_samp  <- res$param[-(1:settings$burn), , drop = FALSE]
  post_means <- colMeans(post_samp)
  post_sds   <- apply(post_samp, 2, sd)
  
  P <- length(settings$beta_true)
  cred_int <- vector("list", P)
  for (j in 1:P) {
    cred_int[[j]] <- quantile(post_samp[, j], probs = c(settings$alpha / 2, 1 - settings$alpha / 2))
  }
  
  output <- list(
    beta_mean = post_means,
    beta_sd   = post_sds,
    tests     = tsts,
    cred_int  = cred_int
  )
  
  if (settings$keep_raw) {
    output$beta_all <- all_samp
    if (!settings$known_acc) {
      output$se <- res$se
      output$sp <- res$sp
    }
  } else if(!settings$known_acc) {
    output$se <- res$se[-(1:settings$burn), , drop = FALSE]
    output$sp <- res$sp[-(1:settings$burn), , drop = FALSE]
  }
  
  return(output)
}


#######################
## RUN REPLICATES    ##
#######################
run_replicates <- function(test_data, settings) {
  start_time <- Sys.time()
  
  cores <- parallel::detectCores(logical=TRUE)
  cl <- parallel::makeCluster(cores)
  
  parallel::clusterEvalQ(cl, {
    dyn.load("Fortran/sampler.dll")
  })
  
  parallel::clusterExport(cl, varlist = c(
    "infer_posterior", "bayes_sampler", "hier.gt.simulation",
    "wls.mh.alg", "mt.ct", "settings"
  ), envir = environment())
  
  sim_list <- parallel::parLapplyLB(cl, test_data, function(d) {
    infer_posterior(d, settings)
  })
  
  parallel::stopCluster(cl)
  
  beta_means_mat <- t(sapply(sim_list, function(s) s$beta_mean))
  beta_sds_mat   <- t(sapply(sim_list, function(s) s$beta_sd))
  tests_vec      <- sapply(sim_list, function(s) s$tests)
  cred_ints      <- lapply(sim_list, function(s) s$cred_int)
  
  out <- list(
    beta_means = beta_means_mat,
    beta_sds   = beta_sds_mat,
    tests      = tests_vec,
    cred_ints  = cred_ints
  )
  
  if (!settings$known_acc) {
    out$se_samps <- lapply(sim_list, `[[`, "se")
    out$sp_samps <- lapply(sim_list, `[[`, "sp")
  }
  
  if (settings$keep_raw) {
    out$beta_raw <- lapply(sim_list, `[[`, "beta_all")
    if (!settings$known_acc) {
      out$se_raw   <- lapply(sim_list, `[[`, "se")
      out$sp_raw   <- lapply(sim_list, `[[`, "sp")
    }
  }
  
  end_time <- Sys.time()
  out$runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  return(compute_results(out, settings))
}


#######################
## COMPUTE RESULTS   ##
#######################
compute_results <- function(output, settings) {
  beta_means <- output$beta_means
  beta_sds   <- output$beta_sds
  tests      <- output$tests
  cred_ints  <- output$cred_ints
  runtime    <- output$runtime
  beta_true  <- settings$beta_true
  alpha      <- settings$alpha
  
  beta_est <- colMeans(beta_means)
  bias     <- beta_est - beta_true
  ssd      <- apply(beta_means, 2, sd)
  ese      <- colMeans(beta_sds)
  avg_tests <- mean(tests)
  pct_reduction <- 100 * (1 - avg_tests / settings$N)
  
  P <- length(beta_true)
  nsim <- nrow(beta_means)
  cp <- numeric(P)
  for (j in 1:P) {
    indicators <- numeric(nsim)
    for (i in 1:nsim) {
      ci_j <- cred_ints[[i]][[j]]
      indicators[i] <- is.numeric(ci_j) && length(ci_j) == 2 &&
        beta_true[j] >= ci_j[1] && beta_true[j] <= ci_j[2]
    }
    cp[j] <- mean(indicators, na.rm = TRUE)
  }
  
  param_summary <- data.frame(
    Parameter   = paste0("beta", seq_along(beta_true) - 1),
    TrueValue   = beta_true,
    EstMean     = beta_est,
    Bias        = bias,
    CP          = cp,
    SSD         = ssd,
    ESE         = ese,
    stringsAsFactors = FALSE
  )
  
  summary <- list(
    settings       = settings,
    runtime        = runtime,
    avg_tests      = avg_tests,
    pct_reduction  = pct_reduction,
    param_summary  = param_summary
  )
  

  if (!settings$known_acc) {
    if (settings$keep_raw) {
      se_post <- lapply(output$se_samps, function(m) m[-(1:settings$burn), , drop = FALSE])
      sp_post <- lapply(output$sp_samps, function(m) m[-(1:settings$burn), , drop = FALSE])
    } else {
      se_post <- output$se_samps
      sp_post <- output$sp_samps
    }
    
    se_means <- do.call(rbind, lapply(se_post, colMeans))
    sp_means <- do.call(rbind, lapply(sp_post, colMeans))
    se_sds   <- do.call(rbind, lapply(se_post, function(m) apply(m, 2, sd)))
    sp_sds   <- do.call(rbind, lapply(sp_post, function(m) apply(m, 2, sd)))
    
    se_post_mean <- colMeans(se_means)
    sp_post_mean <- colMeans(sp_means)
    se_ese       <- colMeans(se_sds)
    sp_ese       <- colMeans(sp_sds)
    
    lb <- alpha / 2
    ub <- 1 - lb
    se_cp <- sp_cp <- numeric(ncol(se_means))
    for (j in seq_len(ncol(se_means))) {
      se_cis <- t(sapply(se_post, function(m) quantile(m[, j], c(lb, ub))))
      sp_cis <- t(sapply(sp_post, function(m) quantile(m[, j], c(lb, ub))))
      se_cp[j] <- mean(settings$se_t[j] >= se_cis[, 1] & settings$se_t[j] <= se_cis[, 2])
      sp_cp[j] <- mean(settings$sp_t[j] >= sp_cis[, 1] & settings$sp_t[j] <= sp_cis[, 2])
    }
    
    assay_labels <- colnames(output$se_samps[[1]])
    assay_ids <- as.integer(gsub("Se_a", "", assay_labels))
    matched_se_t <- settings$se_t[match(assay_ids, unique(settings$assay_id))]
    matched_sp_t <- settings$sp_t[match(assay_ids, unique(settings$assay_id))]
    
    se_bias <- se_post_mean - matched_se_t
    sp_bias <- sp_post_mean - matched_sp_t
    
    acc_summary <- data.frame(
      Assay    = assay_ids,
      True_Se  = matched_se_t,
      Est_Se   = se_post_mean,
      Bias_Se  = se_bias,
      CP_Se    = se_cp,
      SSD_Se   = apply(se_means, 2, sd),
      ESE_Se   = se_ese,
      True_Sp  = matched_sp_t,
      Est_Sp   = sp_post_mean,
      Bias_Sp  = sp_bias,
      CP_Sp    = sp_cp,
      SSD_Sp   = apply(sp_means, 2, sd),
      ESE_Sp   = sp_ese
    )
    
    summary$acc_summary <- acc_summary
  }
  
  if (settings$keep_raw) {
    summary$beta_raw <- output$beta_raw
    if (!settings$known_acc) {
      summary$se_raw <- output$se_raw
      summary$sp_raw <- output$sp_raw
    }
  }
  
  return(summary)
}


## END OF FILE ##