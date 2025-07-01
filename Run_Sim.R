################################
# Run_Sim.R                    #
# Author: Caleb Frankenberger  #
################################


#TODO LIST:
# Se0 and Sp0 separately
# Script for individual run (with debugging)
# Traceplot + Histogram capability - plot convergence
# rbeta in Fortran (test in R first)
# Test SeSp starting values + beta starting values
# (0.99, 0.95, 0.75, 0.50, 0.25)
# 5 beta values as well (0 as one)
# ESE for SeSp
# Any number of assays (L is the notation they used in the paper)
# As much of the sampler in Fortran as we can feasibly do

rm(list=ls())
setwd("C:/Users/caleb/Desktop/Research/bayes-GLM-GT")

library(mnormt)
library(matrixcalc)
library(groupTesting)
library(parallel)

source("Bayes_Sampler.R")
source("MH_WLS.R")
source("Results.R")
source("Simulator.R")


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
  nsim       = 3,
  post_git   = 6000,
  burn       = 1000,
  alpha      = 0.05,
  g          = function(t) exp(t) / (1 + exp(t))
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
  simulate_data(X, settings)
})


##########################
##  RUNNING SIMULATION  ##
##########################
res <- run_replicates(test_data, settings, keep_raw=TRUE)
summary <- summarize_results(res, settings) 

print_results(summary, settings, digits=2)
#print_settings(settings)

plot_trace(res, settings, replicate=2, parameter=1, type="se")
plot_post_hist(res, settings, replicate=2, parameter=1, type="se")

#save_results(summary, settings)

## END OF FILE ##