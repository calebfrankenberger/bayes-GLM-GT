################################
# Simulator.R                  #
# Author: Caleb Frankenberger  #
################################

library(parallel)

#########################
## SIMULATE DATA       ##
#########################
simulate_data <- function(X, settings) {
  p.t <- settings$g(X %*% settings$beta_true)
  sim_out <- hier.gt.simulation(
    N       = nrow(X),
    p       = p.t,
    S       = length(settings$psz),
    psz     = settings$psz,
    Se      = settings$se_t,
    Sp      = settings$sp_t,
    assayID = settings$assay_id
  )
  
  list(X = X, Z = sim_out$gtData, tsts = sim_out$testsExp)
}


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
  
  post_samp  <- res$param[-(1:settings$burn), , drop = FALSE]
  post_means <- colMeans(post_samp)
  post_sds   <- apply(post_samp, 2, sd)
  
  P <- length(settings$beta_true)
  cred_int <- vector("list", P)
  for (j in 1:P) {
    cred_int[[j]] <- quantile(post_samp[, j], probs = c(settings$alpha / 2, 1 - settings$alpha / 2))
  }
  
  extra_out <- list()
  if (!settings$known_acc) {
    extra_out$se <- res$se[-(1:settings$burn), , drop = FALSE]
    extra_out$sp <- res$sp[-(1:settings$burn), , drop = FALSE]
  }
  
  return(c(
    list(
      beta_mean = post_means,
      beta_sd   = post_sds,
      tests     = tsts,
      cred_int  = cred_int
    ),
    extra_out
  ))
}


#######################
## RUN REPLICATES    ##
#######################
run_replicates <- function(test_data, settings) {
  start_time <- Sys.time()
  
  cores <- parallel::detectCores(logical=TRUE)
  cl <- parallel::makeCluster(cores)
  
  parallel::clusterEvalQ(cl, {
    dyn.load("latent_sampler.dll")
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
  
  end_time <- Sys.time()
  out$runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  return(out)
}


## END OF FILE ##