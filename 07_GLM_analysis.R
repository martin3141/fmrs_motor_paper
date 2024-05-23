library(spant)
library(ggplot2)
library(cowplot)

theme_set(theme_bw())

source("extra_functions.R")

# only good data
preproc <- readRDS(file.path("good_fmrs_preproc", "preproc_corrected.rds"))

group_mrs_unproc <- preproc$corrected |> mean_mrs_list()

group_mrs <- group_mrs_unproc |> lb(2) |> crop_spec(xlim = c(4.3, 0.0)) |> 
             bc_als(lambda = 1000)

drift <- group_mrs |> spec_op(xlim = c(1.97, 2.04), operator = "sum") |>
         smooth.spline(spar = 0.7) |> (\(x) x$y)()

stim <- gen_trap_reg(3 * 60, 8 * 60, "stim", group_mrs)

glm_spec_res <- glm_spec_paper(group_mrs, stim)

p1 <- function() plot(glm_spec_res$p_value_log_mrs, y_scale = TRUE,
                      mar = c(2, 3, 1, 3.5), yaxis_lab = "-log10(p-value)",
                      xlim = c(4, 0.2))

stim_nuisance         <- cbind(stim, drift)
glm_spec_nuisance_res <- glm_spec_paper(group_mrs, stim_nuisance)

p2 <- function() plot(glm_spec_nuisance_res$p_value_log_mrs, y_scale = TRUE,
                      mar = c(2, 3, 1, 3.5), yaxis_lab = "-log10(p-value)",
                      xlim = c(4, 0.2))

metab_sim <- sim_basis(c("lac", "glu", "asp"), acq_paras = group_mrs_unproc,
                       pul_seq = seq_slaser_ideal) |> lb(2) |> basis2mrs_data()

ref_spec <- append_dyns(metab_sim, mean_dyns(group_mrs_unproc)) |> 
            scale_spec(operator = "max", xlim = c(4, 0.5), mean_dyns = FALSE) |>
            scale_mrs_amp(c(0.4, 0.4, 0.4, 1))

vlines <- c(0.87, 1.28, 1.35, 2.35, 2.54, 2.74, 2.78, 3.59, 3.23)

p3 <- function() ref_spec |> get_dyns(c(1, 2, 4)) |> 
  stackplot(labels = c("Lac", "Glu", "Mean"),
            y_offset = 35, mar = c(3, 3, 1, 3.5),
            xlim = c(4.0, 0.2))

tiff(file.path("FIGURES", "Fig5.tiff"), width = 1000, height = 1600, res = 200)
plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12, nrow = 3)
dev.off()

stim_lagged <- gen_trap_reg(5 * 60, 8 * 60, "stim", group_mrs) # lagged by 2 mins
stim_lagged <- cbind(stim_lagged, drift)

glm_spec_lagged_res <- glm_spec_paper(group_mrs, stim_lagged)

p4 <- function() plot(glm_spec_lagged_res$p_value_log_mrs, y_scale = TRUE,
                      mar = c(2, 3, 1, 3.5), yaxis_lab = "-log10(p-value)",
                      vline = vlines, vline_lty = 3, xlim = c(4, 0.2))

p5 <- function() ref_spec |> stackplot(labels = c("Lac", "Glu", "Asp", "Mean"),
                                       y_offset = 35, mar = c(3, 3, 1, 3.5),
                                       xlim = c(4, 0.2), vline = vlines,
                                       vline_lty = 3)

tiff(file.path("FIGURES", "Fig6.tiff"), width = 1000, height = 1100, res = 200)
plot_grid(p4, p5, labels = c('A', 'B'), label_size = 12, nrow = 2)
dev.off()

p6 <- function() plot(glm_spec_lagged_res$beta_weight_mrs, y_scale = TRUE,
                      mar = c(2, 3, 1, 3.5), yaxis_lab = "beta weights (au)",
                      vline = vlines, vline_lty = 3, xlim = c(4, 0.2))

tiff(file.path("FIGURES", "FigS4.tiff"), width = 1000, height = 1100, res = 200)
plot_grid(p6, p5, labels = c('A', 'B'), label_size = 12, nrow = 2)
dev.off()
