################################
# results.R                    #
# Author: Caleb Frankenberger  #
################################


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
  
  if (!known_acc && !is.null(summary$acc_summary)) {
    acc <- summary$acc_summary
    acc_rows <- data.frame(
      Parameter   = c(paste0("Se", acc$Assay - 1), paste0("Sp", acc$Assay)),
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
  
  EstMean <- round(param_df$EstMean, digits)
  Bias    <- round(param_df$Bias, digits)
  SSD     <- round(param_df$SSD, digits)
  ESE     <- round(param_df$ESE, digits)
  CP      <- round(param_df$CP, 2)
  
  labels <- paste0(param_df$Parameter, " (", param_df$TrueValue, ")")
  
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
  
  label_col_width <- 14 
  
  cat("=====================================\n")
  cat("Model:", model_name, "\n")
  cat("=====================================\n\n")
  cat(sprintf("Total runtime: %d min %.1f sec\n", mins, secs))
  cat(sprintf("Average # of tests: %.1f (%.1f%% reduction)\n", avg_tests, pct_reduction))
  cat("=====================================\n\n")
  
  header <- sprintf("%-*s", label_col_width, "Parameter")
  for (label in labels) {
    header <- paste0(header, " | ", formatC(label, width = col_width, flag = "-"))
  }
  cat(header, "\n")
  cat(strrep("-", nchar(header)), "\n")
  
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


#######################
## SAVE RESULTS      ##
#######################
save_results <- function(summary, settings, filename = NULL) {
  dir.create('results', showWarnings = FALSE, recursive = TRUE)
  if (is.null(filename)) {
    timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
    filename <- paste0(settings$model, '_', timestamp)
  }
  out_file <- file.path('results', paste0(filename, '.csv'))
  
  combined <- summary$param_summary
  
  if (!is.null(summary$acc_summary)) {
    acc <- summary$acc_summary
    se_df <- data.frame(
      Parameter = paste0('Se', acc$Assay),
      TrueValue = acc$True_Se,
      EstMean   = acc$Est_Se,
      Bias      = acc$Bias_Se,
      CP        = acc$CP_Se,
      SSD       = acc$SSD_Se,
      ESE       = acc$ESE_Se,
      stringsAsFactors = FALSE
    )
    sp_df <- data.frame(
      Parameter = paste0('Sp', acc$Assay),
      TrueValue = acc$True_Sp,
      EstMean   = acc$Est_Sp,
      Bias      = acc$Bias_Sp,
      CP        = acc$CP_Sp,
      SSD       = acc$SSD_Sp,
      ESE       = acc$ESE_Sp,
      stringsAsFactors = FALSE
    )
    combined <- rbind(combined, se_df, sp_df)
  }
  

  meta_df <- data.frame(
    Key = c('Runtime_Sec', 'Avg_Tests', 'Seed', 'Sample_Size',
            'Simulations', 'Pool_Sizes', 'Assay_IDs', 'Known_Accuracy',
            'True_Se', 'True_Sp', 'Initial_Se', 'Initial_Sp',
            'Gibbs_Samples', 'Burn_in', 'Credible_Level'),
    Value = c(
      summary$runtime,
      summary$avg_tests,
      settings$seed,
      settings$N,
      settings$nsim,
      paste(settings$psz, collapse = ';'),
      paste(settings$assay_id, collapse = ';'),
      ifelse(settings$known_acc, 'Yes', 'No'),
      paste(settings$se_t, collapse = ';'),
      paste(settings$sp_t, collapse = ';'),
      if (!settings$known_acc) paste(settings$se_0, collapse = ';') else NA,
      if (!settings$known_acc) paste(settings$sp_0, collapse = ';') else NA,
      settings$post_git,
      settings$burn,
      sprintf('%.1f%%', 100 * (1 - settings$alpha))
    ),
    stringsAsFactors = FALSE
  )
  
  write.table(meta_df,
              file = out_file,
              sep = ',',
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)
  
  write.table('', file = out_file, append = TRUE,
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  header_line <- paste0(names(combined), collapse = ',')
  write(header_line, file = out_file, append = TRUE)
  write.table(combined,
              file = out_file,
              sep = ',',
              row.names = FALSE,
              col.names = FALSE,
              append = TRUE,
              quote = TRUE)
}


########################
##  DIAGNOSTIC PLOTS  ##
########################
plot_trace <- function(results, settings, replicate = 1, parameter = 1, type = c("beta","se","sp")) {
  if (!settings$keep_raw) {
    stop("Trace and posterior plots require settings$keep_raw = TRUE.")
  }
  
  type <- match.arg(type)
  raw  <- switch(type,
                 beta = results$beta_raw,
                 se   = results$se_raw,
                 sp   = results$sp_raw)
  chain <- raw[[replicate]][, parameter]
  burn  <- settings$burn
  
  param_label <- switch(type,
                        beta = paste0("beta_", parameter - 1),
                        se   = paste0("Se_",    parameter),
                        sp   = paste0("Sp_",    parameter))
  
  oldpar <- par(mar = c(4,4,2,1))
  plot(chain,
       type = "l",
       col  = "steelblue",
       lwd  = 1.5,
       xlab = "Iteration",
       ylab = "Value",
       main = sprintf("Traceplot: %s (Replication %d)",
                      param_label, replicate))
  abline(v = burn, col = "red", lty = 2, lwd = 1.5)
  par(oldpar)
}

plot_post_hist <- function(results, settings, replicate = 1, parameter = 1, type = c("beta","se","sp"), breaks = 30) {
  if (!settings$keep_raw) {
    stop("Trace and posterior plots require settings$keep_raw = TRUE.")
  }
  
  type <- match.arg(type)
  raw  <- switch(type,
                 beta = results$beta_raw,
                 se   = results$se_raw,
                 sp   = results$sp_raw)
  chain <- raw[[replicate]][-(1:settings$burn), parameter]
  
  oldpar <- par(mar = c(4,4,2,1))
  hist(chain, breaks = breaks, freq = FALSE,
       col    = "steelblue", border = "black",
       xlab   = "Value",
       ylab   = "Density",
       main   = sprintf("Histogram: %s (Replication %d)",
                        switch(type,
                               beta = paste0("beta_", parameter-1),
                               se   = paste0("Se_",   parameter),
                               sp   = paste0("Sp_",   parameter)),
                        replicate))
  par(oldpar)
}


## END OF FILE ##