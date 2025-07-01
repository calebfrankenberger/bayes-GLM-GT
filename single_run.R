####################
# single_run.R     #
####################

rm(list = ls())
Sys.setenv(OMP_NUM_THREADS = "12")
options(error = recover, warn = 1)

proj_root <- "C:/Users/caleb/Desktop/Research/bayes-GLM-GT"
setwd(proj_root)

library(mnormt)
library(matrixcalc)
library(groupTesting)
library(profvis)

# Load source files
src <- function(f) source(file.path("R", f), chdir = TRUE)
src("bayes_sampler.R")
src("estimator.R")
src("mh_wls.R")
src("results.R")

# Load DLL if needed
dll_path <- file.path("Fortran", paste0("sampler", .Platform$dynlib.ext))
if (!is.loaded("sample"))
  dyn.load(dll_path)

# Settings
settings <- list(
  seed       = 123,
  N          = 5000,
  model      = "M1",
  beta_true  = NULL,
  psz        = c(5, 1),
  assay_id   = c(1, 2),
  se_t       = c(0.95, 0.98),
  sp_t       = c(0.98, 0.99),
  known_acc  = FALSE,
  se_0       = 0.9,
  sp_0       = 0.9,
  nsim       = 1,
  post_git   = 6000,
  burn       = 1000,
  alpha      = 0.05,
  g          = stats::plogis,
  keep_raw   = TRUE
)
settings$beta_true <- switch(settings$model,
                             "M1" = c(-3, 2),
                             "M2" = c(-3, 2, -1),
                             "M3" = c(-3, 2, -0.5))

# Profile the entire run
profvis({
  
  # --- Simulate data ---
  set.seed(settings$seed)
  x1 <- rnorm(settings$N)
  x2 <- rbinom(settings$N, 1, 0.5)
  
  X <- switch(settings$model,
              "M1" = cbind(1, x1),
              "M2" = cbind(1, x1, x2),
              "M3" = cbind(1, x1, x1^2))
  
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
  
  test_data <- list(X = X,
                    Z = sim_out$gtData,
                    tsts = sim_out$testsExp)
  
  # --- Inference ---
  res <- infer_posterior(test_data, settings)
  
  # --- Summarize ---
  summary_out <- compute_results(
    output = list(
      beta_means = matrix(res$beta_mean, nrow = 1),
      beta_sds   = matrix(res$beta_sd,   nrow = 1),
      cred_ints  = list(res$cred_int),
      tests      = res$tests,
      beta_raw   = list(res$beta_all), 
      se_samps   = if (!settings$known_acc) list(res$se),
      sp_samps   = if (!settings$known_acc) list(res$sp)
    ),
    settings = settings
  )
  
  print_results(summary_out, settings, digits = 2)
  
  
}) # END profvis