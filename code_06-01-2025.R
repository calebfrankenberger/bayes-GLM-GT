#####################################
# Author: Caleb Frankenberger       #
# Date: 06-01-2025                  #
#####################################

rm(list=ls())
setwd("C:/Users/caleb/Desktop/Research/06-01-2025")
library(groupTesting)

#####################################
# DATA SIMULATION
#####################################
set.seed(123)
N       <- 1000
S       <- 2
psz     <- c(5, 1)
Se      <- c(.95, .95)
Sp      <- c(.98, .98)
assayID <- c(1, 1)

param.t <- c(-3, 2)
x1 <- rnorm(N, mean=0, sd=1)
X  <- cbind(1, x1)
colnames(X) <- c("Intercept", "Predictor 1")

h <- function(t) exp(t)/(1 + exp(t))
pReg <- h(X %*% param.t)

gtOut <- hier.gt.simulation(N, pReg, S, psz, Se, Sp, assayID)$gtData

#####################################
# SETTING UP Z_mat AND Yt_mat
#####################################

Z_mat <- as.matrix(gtOut[, -(3:5)])  
Zrows <- nrow(Z_mat)
Zcols <- ncol(Z_mat)

Memb <- as.matrix(gtOut[, -(1:5)])  
maxAssign <- max(table(Memb[Memb > 0]))
ytm <- matrix(0L, nrow = N, ncol = maxAssign)
poolIDs <- seq_len(Zrows)

for (i in seq_len(N)) {
  inWhich <- poolIDs[ apply(Memb, 1, function(r) i %in% r) ]
  if (length(inWhich)) {
    ytm[i, seq_along(inWhich)] <- sort(inWhich)
  }
}

Yt_initial <- rep(0L, N)
numPools    <- as.integer(rowSums(ytm > 0))   
Yt_mat      <- cbind(Yt_initial, numPools, ytm)
Ycols       <- ncol(Yt_mat)

#####################################
# BEGIN SAMPLING
#####################################
G        <- 13000
burn     <- 3000
beta.save <- matrix(-999, nrow = G, ncol = length(param.t))

# Tracks # of positive results over every iteration
diag <- numeric(G)

dyn.load("gibbs_bayes.dll")

bta   <- c(-3, 2)
for (g in seq_len(G)) {
  U <- stats::runif(N)
  p_vec  <- h(X %*% bta)
  
  SeSp <- as.matrix(gtOut[, 3:4])
  
  res <- .C("gibbs_bayes",
            as.double(p_vec),
            as.integer(Yt_mat),
            as.integer(Z_mat),
            as.integer(N),
            as.double(SeSp),
            as.integer(Ycols),
            as.integer(Zrows),
            as.integer(Zcols),
            as.double(U)
  )[[2]]
  
  Yt_mat <- matrix(res, nrow = N, ncol = Ycols)
  diag[g] <- sum(Yt_mat[, 1])
  
  # Assume bta is known for now
  bta <- c(-3, 2)
  beta.save[g, ] <- bta
}

dyn.unload("gibbs_bayes.dll")

#####################################
# RESULTS
#####################################
beta.save

# Plot the diag tracker
ymin <- min(diag)
ymax <- max(diag)
plot(diag, type = "l", col = "blue",
     main = "Number of Positive Individuals",
     xlab = "Iteration", ylab = "Num. Positives")

rect(xleft = 0, xright = burn, ybottom = ymin, ytop = ymax,
     col = rgb(1, 0, 0, 0.5), border = NA)
