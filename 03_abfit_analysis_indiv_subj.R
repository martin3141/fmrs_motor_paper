library(spant)
library(ggplot2)
library(cowplot)
library(ggsignif)

theme_set(theme_bw())

# only good data
preproc <- readRDS(file.path("good_fmrs_preproc", "preproc_corrected.rds"))

block_size <- 50

mol_names <- get_1h_brain_basis_names(add = c("gly", "peth"))

basis   <- sim_basis(mol_names, pul_seq = seq_slaser_ideal,
                     acq_paras = preproc$corrected[[1]])

res_file <- "abfit_all_subj.rds"

if (!file.exists(res_file)) {
  fit_res_list <- vector("list", length(preproc$corrected))
  for (n in 1:length(fit_res_list)) {
    dataset <- preproc$corrected[[n]] |> mean_dyn_blocks(block_size)
    fit_res_list[[n]] <- fit_mrs(dataset, basis = basis) |>
      scale_amp_ratio("tCr", use_mean_value = TRUE)
  }
  saveRDS(fit_res_list, res_file)
}

fit_res_list <- readRDS(res_file)

stim_box <- gen_trap_reg(3 * 60, 8 * 60, "stim", fit_res_list[[1]]$data)

fit_res_tab_list <- lapply(fit_res_list, \(x) x$res_tab)

glu_mat  <- sapply(fit_res_tab_list, \(x) x$Glu)
glu_mean <- apply(glu_mat, 1, mean)
glu_sd   <- apply(glu_mat, 1, sd) / (19 ^ 0.5)

glu_mean_perc_change <- glu_mean / glu_mean[1] * 100 - 100
glu_perc_change_sd   <- glu_sd   / glu_mean[1] * 100

glu_df <- data.frame(time = stim_box$time, glu_mean_perc_change,
                     glu_perc_change_sd)

p1 <- ggplot(glu_df, aes(x = time, y = glu_mean_perc_change)) + 
  geom_point() + geom_line() + ylab("Glutamate change (%)") + xlab("Time (s)") +
  annotate("rect", xmin = 3 * 60, xmax = 3 * 60 + 8 * 60, ymin = -Inf,
           ymax = Inf, alpha = 0.4, fill = "red") + ylim(c(-3.5, 6)) +
  geom_errorbar(aes(ymin = glu_mean_perc_change - glu_perc_change_sd,
                    ymax = glu_mean_perc_change + glu_perc_change_sd))

lac_mat  <- sapply(fit_res_tab_list, \(x) x$Lac)
lac_mean <- apply(lac_mat, 1, mean)
lac_sd   <- apply(lac_mat, 1, sd) / (19 ^ 0.5)

lac_mean_perc_change <- lac_mean / lac_mean[1] * 100 - 100
lac_perc_change_sd   <- lac_sd   / lac_mean[1] * 100

lac_df <- data.frame(time = stim_box$time, lac_mean_perc_change,
                     lac_perc_change_sd)

p2 <- ggplot(lac_df, aes(x = time, y = lac_mean_perc_change)) + 
  geom_point() + geom_line() + ylab("Lactate change (%)") + xlab("Time (s)") +
  annotate("rect", xmin = 3 * 60, xmax = 3 * 60 + 8 * 60, ymin = -Inf,
           ymax = Inf, alpha = 0.4, fill = "red") + ylim(c(-40, 40)) +
  geom_errorbar(aes(ymin = lac_mean_perc_change - lac_perc_change_sd,
                    ymax = lac_mean_perc_change + lac_perc_change_sd))

tiff(file.path("FIGURES", "Fig3.tiff"), width = 1500, height = 800, res = 200)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
dev.off()

# Welsh
t.test(lac_mean[3:7], lac_mean[c(1, 2, 8:15)], var.equal = FALSE)
t.test(glu_mean[3:7], glu_mean[c(1, 2, 8:15)], var.equal = FALSE)

# Student's
# n.b. reviewer 1 prefers Student's t-tests, so here they are. I'm not p-hacking!
t.test(lac_mean[3:7], lac_mean[c(1, 2, 8:15)], var.equal = TRUE)
t.test(glu_mean[3:7], glu_mean[c(1, 2, 8:15)], var.equal = TRUE)

asp_mat  <- sapply(fit_res_tab_list, \(x) x$Asp)
asp_mean <- apply(asp_mat, 1, mean)
asp_sd   <- apply(asp_mat, 1, sd) / (19 ^ 0.5)

asp_mean_perc_change <- asp_mean / asp_mean[1] * 100 - 100
asp_perc_change_sd   <- asp_sd   / asp_mean[1] * 100

asp_df <- data.frame(time = stim_box$time, asp_mean_perc_change,
                     asp_perc_change_sd)

state      <- rep("REST", 15)
state[3:7] <- "TASK"
lac_tab <- data.frame(Lac = lac_mean_perc_change, state = state)

p3 <- ggplot(lac_tab, aes(x = state, y = Lac)) + geom_point() +
  geom_signif(comparisons = list(c("REST", "TASK")), test = "t.test",
              map_signif_level = function(p) sprintf("p = %.2g", p),
              textsize = 3) + xlab(NULL) + ylab("Lactate change (%)")

glu_tab <- data.frame(Glu = glu_mean_perc_change, state = state)
p4 <- ggplot(glu_tab, aes(x = state, y = Glu)) + geom_point() +
  geom_signif(comparisons = list(c("REST", "TASK")), test = "t.test",
              map_signif_level = function(p) sprintf("p = %.2g", p),
              textsize = 3) + xlab(NULL) + ylab("Glutamate change (%)")

tiff(file.path("FIGURES", "FigS3.tiff"), width = 1500, height = 700, res = 200)
plot_grid(p3, p4, labels = c('A', 'B'), label_size = 12)
dev.off()

# linear mixed effects model and plots suggested by reviewer 1
lac_df_long <- data.frame(lac = as.numeric(lac_mat), time = stim_box$time,
                          state = factor(state),
                          sub = factor(rep(preproc$labels, each = 15)))

p5 <- ggplot(lac_df_long, aes(x = time, y = lac, col = sub)) + geom_line() +
      theme(legend.position="none") + ylab("Lac / tCr") + xlab("Time (s)") +
      annotate("rect", xmin = 3 * 60, xmax = 3 * 60 + 8 * 60, ymin = -Inf,
               ymax = Inf, alpha = 0.4, fill = "red") 

glu_df_long <- data.frame(glu = as.numeric(glu_mat), time = stim_box$time,
                          state = factor(state),
                          sub = factor(rep(preproc$labels, each = 15)))

p6 <- ggplot(glu_df_long, aes(x = time, y = glu, col = sub)) + geom_line() + 
      theme(legend.position="none") + ylab("Glu / tCr") + xlab("Time (s)") +
      annotate("rect", xmin = 3 * 60, xmax = 3 * 60 + 8 * 60, ymin = -Inf,
               ymax = Inf, alpha = 0.4, fill = "red") 

tiff(file.path("FIGURES", "FigS4.tiff"), width = 1500, height = 700, res = 200)
plot_grid(p5, p6, labels = c('A', 'B'), label_size = 12)
dev.off()


library(nlme)
lme(lac ~ state, random =~1|sub, lac_df_long) |> anova()
lme(glu ~ state, random =~1|sub, glu_df_long) |> anova()


tiff(file.path("FIGURES", "mean_spectra_good.tiff"), width = 1000, height = 1200, res = 200)
preproc$corrected |> mean_dyns() |> stackplot(xlim = c(4, 0.5), y_offset = 12,
                                              label = preproc$labels)
dev.off()



break

lac_mat |> image()


mean_proc <- preproc$mean_dataset |> mean_dyn_blocks(block_size)


fit_res <- fit_mrs(mean_proc, basis = basis)

stim_box  <- gen_stim(3 * 60, 8 * 60, "stim_box", mean_proc, "box")
fit_res$res_tab$Time <- stim_box$time

fit_res$res_tab$Lac_perc_change <- fit_res$res_tab$Lac /
                                   fit_res$res_tab$Lac[1] * 100 - 100

fit_res$res_tab$Lac_sd <- fit_res$res_tab$Lac.sd /
                          fit_res$res_tab$Lac[1] * 100

fit_res$res_tab$Glu_perc_change <- fit_res$res_tab$Glu /
                                   fit_res$res_tab$Glu[1] * 100 - 100

fit_res$res_tab$Glu_sd <- fit_res$res_tab$Glu.sd /
                          fit_res$res_tab$Glu[1] * 100

fit_res$res_tab$tCho_perc_change <- fit_res$res_tab$tCho /
                                    fit_res$res_tab$tCho[1] * 100 - 100

fit_res$res_tab$tCho_sd <- fit_res$res_tab$tCho.sd /
                           fit_res$res_tab$tCho[1] * 100

p1 <- ggplot(fit_res$res_tab, aes(x = Time, y = Glu_perc_change)) + 
  geom_point() + geom_line() + ylab("Glutamate change (%)") + xlab("Time (s)") +
  annotate("rect", xmin = 3 * 60, xmax = 3 * 60 + 8 * 60, ymin = -Inf,
           ymax = Inf, alpha = 0.4, fill = "red") + ylim(c(-20, 35)) +
  geom_errorbar(aes(ymin = Glu_perc_change - Glu_sd,
                    ymax = Glu_perc_change + Glu_sd))

p2 <- ggplot(fit_res$res_tab, aes(x = Time, y = Lac_perc_change)) + 
  geom_point() + geom_line() + ylab("Lactate change (%)") + xlab("Time (s)") +
  annotate("rect", xmin = 3 * 60, xmax = 3 * 60 + 8 * 60, ymin = -Inf,
           ymax = Inf, alpha = 0.4, fill = "red") + ylim(c(-20, 35)) +
  geom_errorbar(aes(ymin = Lac_perc_change - Lac_sd,
                    ymax = Lac_perc_change + Lac_sd))


tiff(file.path("FIGURES", "Fig3.tiff"), width = 1500, height = 800, res = 200)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
dev.off()



