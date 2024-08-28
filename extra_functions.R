# warning, any column named time in regressor_df will be removed
glm_spec_paper <- function(mrs_data, regressor_df) {
  
  # needs to be a FD operation
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  mrs_mat <- Re(mrs_data2mat(mrs_data))
  
  # drop the time column if present
  regressor_df<- regressor_df[, !names(regressor_df) %in% c("time"),
                              drop = FALSE]
  
  lm_res_list <- vector("list", ncol(mrs_mat))
  for (n in 1:ncol(mrs_mat)) {
    lm_res_list[[n]] <- summary(lm(mrs_mat[, n] ~ ., regressor_df))
  }
  
  # extract stats
  get_glm_stat <- function(x, name) x$coefficients[-1, name, drop = FALSE]
  
  beta_weight <- as.data.frame(t(as.data.frame(sapply(lm_res_list, get_glm_stat,
                                               "Estimate", simplify = FALSE))))
  p_value <- as.data.frame(t(as.data.frame(sapply(lm_res_list, get_glm_stat,
                                           "Pr(>|t|)", simplify = FALSE))))
  
  row.names(beta_weight) <- NULL
  row.names(p_value)     <- NULL
  
  ppm_sc      <-  ppm(mrs_data)
  beta_weight <-  cbind(ppm = ppm_sc, beta_weight)
  p_value_log <- -log10(p_value)
  p_value     <-  cbind(ppm = ppm_sc, p_value)
  p_value_log <-  cbind(ppm = ppm_sc, p_value_log)
  
  p_value_log_mrs <- mat2mrs_data(t(p_value_log[, -1]), fs = fs(mrs_data),
                                  ft = mrs_data$ft, ref = mrs_data$ref,
                                  nuc = mrs_data$nuc, fd = TRUE)
  
  p_value_mrs     <- mat2mrs_data(t(p_value[, -1]), fs = fs(mrs_data),
                                  ft = mrs_data$ft, ref = mrs_data$ref,
                                  nuc = mrs_data$nuc, fd = TRUE)
  
  beta_weight_mrs <- mat2mrs_data(t(beta_weight[, -1]), fs = fs(mrs_data),
                                  ft = mrs_data$ft, ref = mrs_data$ref,
                                  nuc = mrs_data$nuc, fd = TRUE)
  
  return(list(beta_weight = beta_weight, p_value = p_value,
              p_value_log = p_value_log, p_value_log_mrs = p_value_log_mrs,
              p_value_mrs = p_value_mrs, beta_weight_mrs = beta_weight_mrs,
              lm_res_list = lm_res_list))
}


preproc_fmrs_single_paper <- function(path, label = NULL, output_dir = NULL) {
  
  # TODO combine coils if needed, make the noise region a parameter
  # TODO deal with GE style data with wref included in the same file
  
  mrs_data <- read_mrs(path)
  
  if (is.null(output_dir)) {
    output_dir <- getwd()
  } else {
    if (!dir.exists(output_dir)) dir.create(output_dir)
  }
  
  if (is.null(label)) {
    label <- basename(path)
    label <- tools::file_path_sans_ext(label)
    label <- tools::file_path_sans_ext(label)
  }
  
  mrs_rats  <- rats(mrs_data, xlim = c(4, 1.9), zero_freq_shift_t0 = TRUE,
                    ret_corr_only = FALSE)
  
  # perform simple baseline offset corrected based on a noisy spectral region
  mrs_rats$corrected <- bc_constant(mrs_rats$corrected, xlim = c(-0.5, -2.5))
  
  mean_mrs  <- mean_dyns(mrs_rats$corrected)
    
  # frequency and phase correct the mean spectrum
  ref <- sim_resonances(acq_paras = mrs_data, freq = c(2.01, 3.03, 3.22),
                        amp = 1, lw = 4, lg = 0)
  
  res <- rats(mean_mrs, ref, xlim = c(4, 1.9), p_deg = 4, ret_corr_only = FALSE)
  
  # apply mean spectrum phase and shift to the single shots
  mrs_proc <- phase(mrs_rats$corrected, -as.numeric(res$phases))
  mrs_proc <- shift(mrs_proc, -as.numeric(res$shifts), units = "hz")
  
  mrs_uncorr <-  phase(mrs_data,   -as.numeric(res$phases))
  mrs_uncorr <-  shift(mrs_uncorr,
                       -as.numeric(res$shifts) - mean(mrs_rats$shifts),
                       units = "hz")
  
  mean_uncorr <- mean_dyns(mrs_uncorr)
  
  snr         <- as.numeric(calc_spec_snr(mrs_proc))
  
  # single shot SNR
  ss_median_snr <- median(snr)
  
  dyn_peak_info <- peak_info(mrs_proc, xlim = c(1.8, 2.2))
  lw_ppm        <- as.numeric(dyn_peak_info$fwhm_ppm)
  tnaa_height   <- as.numeric(dyn_peak_info$height)
  
  if (anyNA(lw_ppm)) {
    lw_ppm_smo <- rep(NA, length(lw_ppm))
  } else {
    lw_ppm_smo <- smooth.spline(lw_ppm, spar = 0.8)$y
  }
  
  diag_table <- data.frame(dynamics = 1:length(snr),
                           shifts_hz = as.numeric(mrs_rats$shifts),
                           phases = as.numeric(mrs_rats$phases),
                           snr = snr, lw_ppm = lw_ppm,
                           lw_ppm_smo = lw_ppm_smo, tnaa_height = tnaa_height)
  
  # scale data to the tCr peak
  amp <- spec_op(zf(res$corrected), xlim = c(2.9, 3.1), operator = "max-min")
  amp <- as.numeric(amp)
  mrs_proc      <- scale_mrs_amp(mrs_proc, 1 / amp)
  mrs_uncorr    <- scale_mrs_amp(mrs_uncorr, 1 / amp)
  res$corrected <- scale_mrs_amp(res$corrected, 1 / amp)
  mean_uncorr   <- scale_mrs_amp(mean_uncorr, 1 / amp)
  
  # mean spec SNR and LW
  mean_corr_spec_snr   <- calc_spec_snr(res$corrected)
  corr_peak_info       <- peak_info(res$corrected, xlim = c(1.8, 2.2))
  mean_corr_spec_lw    <- corr_peak_info$fwhm_ppm 
  mean_uncorr_spec_snr <- calc_spec_snr(mean_uncorr)
  uncorr_peak_info     <- peak_info(bc_constant(mean_uncorr,
                                                xlim = c(-0.5, -2.5)),
                                    xlim = c(1.8, 2.2))
  mean_uncorr_spec_lw  <- uncorr_peak_info$fwhm_ppm
  
  # measure the dynamic fluctuation range
  mrs_proc_smoothed <- smooth_dyns(crop_spec(lb(mrs_proc, 5)), 10)
  mrs_mean_sub      <- sub_mean_dyns(mrs_proc_smoothed)
  mrs_mean_sub_bc   <- bc_poly(mrs_mean_sub, 2)
  dfr               <- diff(range(Re(mrs_data2mat(mrs_mean_sub))))
  
  summary_diags <- c(mean_corr_spec_snr = mean_corr_spec_snr,
                     mean_corr_spec_lw = mean_corr_spec_lw,
                     mean_uncorr_spec_snr = mean_uncorr_spec_snr,
                     mean_uncorr_spec_lw = mean_uncorr_spec_lw,
                     ss_median_spec_snr = ss_median_snr,
                     lw_ppm_smo_range = diff(range(lw_ppm_smo)),
                     shift_hz_range = diff(range(diag_table$shifts_hz)),
                     dfr = dfr)
  
  res <- list(corrected = mrs_proc, uncorrected = mrs_uncorr,
              mean_corr = res$corrected, mean_uncorr = mean_uncorr,
              diag_table = diag_table, summary_diags = summary_diags,
              mrs_mean_sub = mrs_mean_sub, mrs_mean_sub_bc = mrs_mean_sub_bc)
  
  rmarkdown::render("single_subject_fmrs_qa.Rmd",
                    params = list(data = res, label = label),
                    output_file = file.path(output_dir, label))
  
  return(res)
}

preproc_fmrs_dataset_paper <- function(paths, labels = NULL,
                                       output_dir = "fmrs_analysis") {
  
  if (is.null(labels)) {
    labels <- basename(paths)
    labels <- tools::file_path_sans_ext(labels)
    labels <- tools::file_path_sans_ext(labels)
  }
  
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # detect non unique labels and quit
  if (any(table(labels) > 1)) stop("Labels are non-unique.")
  
  tot_num <- length(paths)
  preproc_res_list <- vector(mode = "list", length = tot_num)
  for (n in 1:tot_num) {
    cat(c("Processing ", n, " of ", tot_num, " : ",
          labels[n],"\n"), sep = "")
    preproc_res_list[[n]] <- preproc_fmrs_single_paper(paths[n], labels[n],
                                                       output_dir)
  }
  
  preproc_summary <- data.frame(t(sapply(preproc_res_list,
                                         (\(x) x$summary_diags))))
  
  preproc_summary <- cbind(labels, preproc_summary)
  
  corrected_list <- lapply(preproc_res_list, \(x) x$corrected)
  if (tot_num > 1) {
    mean_dataset <- mean_mrs_list(corrected_list)
  } else {
    mean_dataset <- corrected_list
  }
  
  res <- list(res_list = preproc_res_list, summary = preproc_summary,
              mean_dataset = mean_dataset)
  
  rmarkdown::render("dataset_summary_fmrs_qa.Rmd",
                    params = list(data = res),
                    output_file = file.path(output_dir, "dataset_summary"))
  
  saveRDS(res, file.path(output_dir, "preproc_full.rds"))
  
  cut_res <- list(corrected = corrected_list, mean_dataset = mean_dataset,
                  labels = labels)
  
  saveRDS(cut_res, file.path(output_dir, "preproc_corrected.rds"))
  
  return(res)
}