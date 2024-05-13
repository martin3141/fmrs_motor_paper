library(spant)
library(ggplot2)
library(cowplot)

theme_set(theme_bw())

# only good data
preproc <- readRDS(file.path("good_fmrs_preproc", "preproc_corrected.rds"))

# SPECTROGRAM PLOTS

p1 <- function() preproc$corrected[[4]] |> lb(4) |> sub_mean_dyns() |> crop_spec() |> 
  bc_poly(1) |> image(legend = TRUE, xlim = c(4, 0.5),
                      legend.lab = "Spectral Intensity", legend.mar = 6.5,
                      legend.line = 3.5, plot_dim = "time_sec", cex.axis = 0.8,
                      legend.cex = 0.8, axis.args = list(cex.axis = 0.7),
                      hline = c(3 * 60, 11 * 60))

p2 <- function() preproc$mean_dataset |> lb(4) |> sub_mean_dyns() |> crop_spec() |>
  bc_poly(1) |> image(legend = TRUE, xlim = c(4, 0.5),
                      legend.lab = "Spectral Intensity", legend.mar = 6.5,
                      legend.line = 3.5, plot_dim = "time_sec", cex.axis = 0.8,
                      legend.cex = 0.8, axis.args = list(cex.axis = 0.7),
                      hline = c(3 * 60, 11 * 60))

# LINEWIDTH PLOTS

# confirm spectra are comparable
append_dyns(preproc$corrected[[4]] |> mean_dyns(), 
            preproc$mean_dataset |> mean_dyns()) |>
  stackplot(xlim = c(4, 0.5), y_offset = 20)

# read in all preproc results
preproc_full <- readRDS(file.path("good_fmrs_preproc", "preproc_full.rds"))

# get the linewidth estimates
lw_mat <- sapply(preproc_full$res_list, (\(x) x$diag_table$lw_ppm))
lw_mat_smo <- apply(lw_mat, 2, \(x) smooth.spline(x, spar = 0.8)$y)
lw_mat_smo_perc_change <- apply(lw_mat_smo, 2, (\(x) 100 * ((x / x[1]) - 1)))

stim_box <- gen_trap_rf(3 * 60, 8 * 60, "stim", preproc$mean_dataset)

lw_smo_df      <- lw_mat_smo |> data.frame() |> tidyr::gather()
lw_smo_df$time <- stim_box$time

# from the mean spectrum
mean_spec_peak_info <- preproc$mean_dataset |> 
                       peak_info(xlim = c(xlim = c(1.8, 2.2)))
mean_spec_lw_smo    <- smooth.spline(mean_spec_peak_info$fwhm_ppm, spar = 0.8)$y
mean_spec_lw_smo_df <- data.frame(t = stim_box$time, mean_spec_lw_smo)

p3 <- ggplot(lw_smo_df) +
  geom_line(aes(x = time, y = value, col = key), show.legend = FALSE) +
  annotate("rect", xmin = 3 * 60, xmax = 3 * 60 + 8 * 60, ymin = -Inf,
           ymax = Inf, alpha = 0.4, fill = "red") + 
  geom_line(aes(x = t, y = mean_spec_lw_smo), mean_spec_lw_smo_df,
            linewidth = 1.2, linetype = "dashed") + xlab("Time (s)") + 
  ylab("Linewidth (ppm)")

lw_smo_pc_df      <- lw_mat_smo_perc_change |> data.frame() |> tidyr::gather()
lw_smo_pc_df$time <- stim_box$time

# calculate percent change
calc_pc <- function(x) 100 * ((x / x[1]) - 1)

mean_spec_lw_smo_df$mean_spec_lw_smo_pc <- calc_pc(mean_spec_lw_smo_df$mean_spec_lw_smo)

p4 <- ggplot(lw_smo_pc_df) + 
  geom_line(aes(x = time, y = value, col = key), show.legend = FALSE) +
  annotate("rect", xmin = 3 * 60, xmax = 3 * 60 + 8 * 60, ymin = -Inf,
           ymax = Inf, alpha = 0.4, fill = "red") + 
  geom_line(aes(x = t, y = mean_spec_lw_smo_pc), mean_spec_lw_smo_df,
            linewidth = 1.2, linetype = "dashed") + xlab("Time (s)") +
  ylab("Linewidth change (%)")

tiff(file.path("FIGURES", "Fig4.tiff"), width = 1800, height = 1400, res = 200)
plot_grid(p1, p2, p3, p4, labels = c('A', 'B', 'C', 'D'), label_size = 12)
dev.off()
