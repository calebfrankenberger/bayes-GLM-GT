################################
# run.R                        #
# Author: Caleb Frankenberger  #
################################


#TODO LIST:
# rbeta in Fortran (test in R first)
# Test SeSp starting values + beta starting values
# (0.99, 0.95, 0.75, 0.50, 0.25)
# 5 beta values as well (0 as one)
# Any number of assays (L is the notation they used in the paper)
# As much of the sampler in Fortran as we can feasibly do

rm(list=ls())
setwd("C:/Users/caleb/Desktop/Research/bayes-GLM-GT")

library(mnormt)
library(matrixcalc)
library(groupTesting)
library(parallel)

source("R/bayes_sampler.R")
source("R/estimator.R")
source("R/mh_wls.R")
source("R/results.R")


###############
##  SETTINGS ##
###############
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
  nsim       = 10,
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


##########################
##  SIMULATE TEST DATA  ##
##########################
set.seed(settings$seed)
test_data <- lapply(1:settings$nsim, function(i) {
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
  
  list(X = X, Z = sim_out$gtData, tsts = sim_out$testsExp)
})


##########################
##  RUNNING SIMULATION  ##
##########################
res <- run_replicates(test_data, settings)

print_results(res, settings, digits=2)
#print_settings(settings)

plot_trace(res, settings, replicate=2, parameter=1, type="se")
plot_post_hist(res, settings, replicate=2, parameter=1, type="se")

#save_results(res, settings)

## END OF FILE ##