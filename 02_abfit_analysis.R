library(spant)
library(ggplot2)
library(cowplot)

theme_set(theme_bw())

# only good data
preproc <- readRDS(file.path("good_fmrs_preproc", "preproc_corrected.rds"))

block_size <- 50

mean_proc <- preproc$mean_dataset |> mean_dyn_blocks(block_size)

mol_names <- get_1h_brain_basis_names(add = c("gly", "peth"))

basis   <- sim_basis(mol_names, pul_seq = seq_slaser_ideal,
                     acq_paras = mean_proc)

fit_res <- fit_mrs(mean_proc, basis = basis)

stim_box <- gen_trap_reg(3 * 60, 8 * 60, "stim", mean_proc)
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
           ymax = Inf, alpha = 0.4, fill = "red") + ylim(c(-5, 6)) +
  geom_errorbar(aes(ymin = Glu_perc_change - Glu_sd,
                    ymax = Glu_perc_change + Glu_sd))

p2 <- ggplot(fit_res$res_tab, aes(x = Time, y = Lac_perc_change)) + 
  geom_point() + geom_line() + ylab("Lactate change (%)") + xlab("Time (s)") +
  annotate("rect", xmin = 3 * 60, xmax = 3 * 60 + 8 * 60, ymin = -Inf,
           ymax = Inf, alpha = 0.4, fill = "red") + ylim(c(-20, 35)) +
  geom_errorbar(aes(ymin = Lac_perc_change - Lac_sd,
                    ymax = Lac_perc_change + Lac_sd))

tiff(file.path("FIGURES", "Fig2.tiff"), width = 1500, height = 800, res = 200)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
dev.off()

t.test(fit_res$res_tab$Lac[3:7], fit_res$res_tab$Lac[c(1, 2, 8:15)])
t.test(fit_res$res_tab$Glu[3:7], fit_res$res_tab$Glu[c(1, 2, 8:15)])

t.test(fit_res$res_tab$Asp[3:7], fit_res$res_tab$Asp[c(1, 2, 8:15)])
t.test(fit_res$res_tab$Asp[4:8], fit_res$res_tab$Asp[c(1, 2, 3, 9:15)])


fit_res$res_tab$Asp_perc_change <- fit_res$res_tab$Asp /
                                   fit_res$res_tab$Asp[1] * 100 - 100

fit_res$res_tab$Asp_sd <- fit_res$res_tab$Asp.sd /
                          fit_res$res_tab$Asp[1] * 100

tiff(file.path("FIGURES", "S3.tiff"), width = 1000, height = 800, res = 200)
ggplot(fit_res$res_tab, aes(x = Time, y = Asp_perc_change)) + 
  geom_point() + geom_line() + ylab("Asp change (%)") + xlab("Time (s)") +
  annotate("rect", xmin = 3 * 60, xmax = 3 * 60 + 8 * 60, ymin = -Inf,
           ymax = Inf, alpha = 0.4, fill = "red") + ylim(c(-15, 15)) +
  geom_errorbar(aes(ymin = Asp_perc_change - Asp_sd,
                    ymax = Asp_perc_change + Asp_sd))
dev.off()