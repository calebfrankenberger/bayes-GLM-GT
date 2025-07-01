################################
# Results.R                    #
# Author: Caleb Frankenberger  #
################################


#######################
## SUMMARIZE RESULTS ##
#######################
# Eventually this should output all the RAW summary data
summarize_results <- function(results, settings) {
  beta_means <- results$beta_means
  beta_sds   <- results$beta_sds
  tests      <- results$tests
  cred_ints  <- results$cred_ints
  runtime    <- results$runtime
  beta_true  <- settings$beta_true
  se_t       <- settings$se_t
  sp_t       <- settings$sp_t
  alpha      <- settings$alpha
  
  beta_est <- colMeans(beta_means)
  bias     <- beta_est - beta_true
  ssd      <- apply(beta_means, 2, sd)
  ese      <- colMeans(beta_sds)
  avg_tests <- mean(tests)
  
  P <- length(beta_true)
  nsim <- nrow(beta_means)
  cp <- numeric(P)
  for (j in 1:P) {
    indicators <- numeric(nsim)
    for (i in 1:nsim) {
      ci_j <- cred_ints[[i]][[j]]
      if (is.numeric(ci_j) && length(ci_j) == 2 && all(!is.na(ci_j))) {
        indicators[i] <- beta_true[j] >= ci_j[1] && beta_true[j] <= ci_j[2]
      } else {
        indicators[i] <- NA
      }
    }
    cp[j] <- mean(indicators, na.rm = TRUE)
  }
  
  df_par <- data.frame(
    Parameter   = paste0("beta", seq_along(beta_true) - 1),
    TrueValue   = beta_true,
    EstMean     = beta_est,
    Bias        = bias,
    CP          = cp,
    SSD         = ssd,
    ESE         = ese,
    stringsAsFactors = FALSE
  )
  
  output <- list(
    param_summary = df_par,
    avg_tests     = avg_tests,
    runtime       = runtime
  )
  
  if (!is.null(results$se_samps) && !is.null(results$sp_samps)) {
    se_means <- do.call(rbind, lapply(results$se_samps, colMeans))
    sp_means <- do.call(rbind, lapply(results$sp_samps, colMeans))
    se_sds   <- do.call(rbind, lapply(results$se_samps, function(m) apply(m, 2, sd)))
    sp_sds   <- do.call(rbind, lapply(results$sp_samps, function(m) apply(m, 2, sd)))
    
    se_post_mean <- colMeans(se_means)
    sp_post_mean <- colMeans(sp_means)
    se_ese       <- colMeans(se_sds)
    sp_ese       <- colMeans(sp_sds)
    
    lb <- alpha / 2
    ub <- 1 - lb
    se_cp <- sp_cp <- numeric(ncol(se_means))
    
    for (j in seq_len(ncol(se_means))) {
      se_cis <- t(sapply(results$se_samps, function(m) quantile(m[, j], c(lb, ub))))
      sp_cis <- t(sapply(results$sp_samps, function(m) quantile(m[, j], c(lb, ub))))
      se_cp[j] <- mean(se_t[j] >= se_cis[, 1] & se_t[j] <= se_cis[, 2])
      sp_cp[j] <- mean(sp_t[j] >= sp_cis[, 1] & sp_t[j] <= sp_cis[, 2])
    }
    
    assay_labels <- colnames(results$se_samps[[1]])
    assay_ids <- as.integer(gsub("Se_a", "", assay_labels))
    
    matched_se_t <- se_t[match(assay_ids, unique(settings$assay_id))]
    matched_sp_t <- sp_t[match(assay_ids, unique(settings$assay_id))]
    
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
    
    output$Se_summary  <- se_post_mean
    output$Sp_summary  <- sp_post_mean
    output$acc_summary <- acc_summary
  }
  
  return(output)
}


#######################
## PRINT RESULTS     ##
#######################
print_results <- function(summary, settings, digits = 4) {
  known_acc   <- settings$known_acc
  runtime     <- summary$runtime
  avg_tests   <- summary$avg_tests
  N           <- settings$N
  model_name  <- settings$model
  
  mins <- floor(runtime / 60)
  secs <- round(runtime %% 60, 1)
  pct_reduction <- round(100 * (1 - avg_tests / N), 1)
  
  param_df <- summary$param_summary
  
  # Append Se/Sp rows if needed
  if (!known_acc && !is.null(summary$acc_summary)) {
    acc <- summary$acc_summary
    acc_rows <- data.frame(
      Parameter   = c(paste0("Se", acc$Assay - 1), paste0("Sp", acc$Assay - 1)),
      TrueValue   = c(acc$True_Se, acc$True_Sp),
      EstMean     = c(acc$Est_Se, acc$Est_Sp),
      Bias        = c(acc$Bias_Se, acc$Bias_Sp),
      CP          = c(acc$CP_Se, acc$CP_Sp),
      SSD         = c(acc$SSD_Se, acc$SSD_Sp),
      ESE         = c(acc$ESE_Se, acc$ESE_Sp),
      stringsAsFactors = FALSE
    )
    param_df <- rbind(param_df[, names(acc_rows)], acc_rows)
  }
  
  # Compute derived columns
  EstMean <- round(param_df$EstMean, digits)
  Bias    <- round(param_df$Bias, digits)
  SSD     <- round(param_df$SSD, digits)
  ESE     <- round(param_df$ESE, digits)
  CP      <- round(param_df$CP, 2)
  
  # Labels like "beta0 (-3)", "Se0 (0.95)", etc.
  labels <- paste0(param_df$Parameter, " (", param_df$TrueValue, ")")
  
  # Compute column width
  est_strs  <- formatC(EstMean, digits = digits, format = "f")
  bias_strs <- mapply(function(b, cp) {
    if (is.na(cp)) formatC(b, digits = digits, format = "f")
    else sprintf(paste0("%.", digits, "f (%0.2f)"), b, cp)
  }, Bias, CP)
  ssd_strs <- mapply(function(s, e) {
    if (is.na(e)) formatC(s, digits = digits, format = "f")
    else sprintf(paste0("%.", digits, "f (%0.2f)"), s, e)
  }, SSD, ESE)
  
  col_width <- max(
    nchar(labels),
    nchar(est_strs),
    nchar(bias_strs),
    nchar(ssd_strs)
  )
  
  label_col_width <- 14  # width for row labels like "Bias (CP95)"
  
  # Header
  cat("=====================================\n")
  cat("Model:", model_name, "\n")
  cat("=====================================\n\n")
  cat(sprintf("Total runtime: %d min %.1f sec\n", mins, secs))
  cat(sprintf("Average # of tests: %.1f (%.1f%% reduction)\n", avg_tests, pct_reduction))
  cat("=====================================\n\n")
  
  # Column header
  header <- sprintf("%-*s", label_col_width, "Parameter")
  for (label in labels) {
    header <- paste0(header, " | ", formatC(label, width = col_width, flag = "-"))
  }
  cat(header, "\n")
  cat(strrep("-", nchar(header)), "\n")
  
  # Helper to print each row
  print_row <- function(name, values) {
    line <- sprintf("%-*s", label_col_width, name)
    for (val in values) {
      line <- paste0(line, " | ", formatC(val, width = col_width, flag = "-"))
    }
    cat(line, "\n")
  }
  
  print_row("Est.", est_strs)
  print_row("Bias (CP95)", bias_strs)
  print_row("SSD (ESE)", ssd_strs)
  
  cat(strrep("=", nchar(header)), "\n")
}


#######################
## PRINT SETTINGS    ##
#######################
print_settings <- function(settings) {
  cat("\nReproducibility Settings:\n")
  cat("==========================\n")
  
  fmt <- "%-20s : %s\n"  
  
  cat(sprintf(fmt, "Seed",               settings$seed))
  cat(sprintf(fmt, "Sample size",    settings$N))
  cat(sprintf(fmt, "Simulations", settings$nsim))
  cat(sprintf(fmt, "Model",              settings$model))
  cat(sprintf(fmt, "Pool sizes",         paste(settings$psz, collapse = ", ")))
  cat(sprintf(fmt, "Assay IDs",          paste(settings$assay_id, collapse = ", ")))
  cat(sprintf(fmt, "Known Accuracy",     ifelse(settings$known_acc, "Yes", "No")))
  
  cat(sprintf(fmt, "True Se",            paste(settings$se_t, collapse = ", ")))
  cat(sprintf(fmt, "True Sp",            paste(settings$sp_t, collapse = ", ")))
  
  if (!settings$known_acc) {
    cat(sprintf(fmt, "Initial Se",       paste(settings$se_0, collapse = ", ")))
    cat(sprintf(fmt, "Initial Sp",       paste(settings$sp_0, collapse = ", ")))
  }
  
  cat(sprintf(fmt, "Gibbs samples",      settings$post_git))
  cat(sprintf(fmt, "Burn-in",            settings$burn))
  cat(sprintf(fmt, "Credible level",     sprintf("%.1f%%", 100 * (1 - settings$alpha))))
  
  cat("==========================\n")
}


## END OF FILE ##