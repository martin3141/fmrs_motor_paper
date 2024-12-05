library(spant)
library(ggplot2)
library(cowplot)
library(ggsignif)

theme_set(theme_bw())

# only good data
preproc <- readRDS(file.path("good_fmrs_preproc", "preproc_corrected.rds"))

block_size <- 50

mean_proc <- preproc$mean_dataset |> mean_dyn_blocks(block_size)

mol_names <- get_1h_brain_basis_names(add = c("gly", "peth"))

basis   <- sim_basis(mol_names, pul_seq = seq_slaser_ideal,
                     acq_paras = mean_proc)

fit_res <- fit_mrs(mean_proc, basis = basis, method = "LCMODEL",
                   opts = c("NSIMUL=0"))

stim_box <- gen_trap_reg(3 * 60, 8 * 60, "stim_box", mean_proc)
fit_res$res_tab$Time <- stim_box$time

fit_res$res_tab$Lac_perc_change <- fit_res$res_tab$Lac /
                                   fit_res$res_tab$Lac[1] * 100 - 100

fit_res$res_tab$Lac_sd <- fit_res$res_tab$Lac.sd /
                          fit_res$res_tab$Lac[1] * 100

fit_res$res_tab$Glu_perc_change <- fit_res$res_tab$Glu /
                                   fit_res$res_tab$Glu[1] * 100 - 100

fit_res$res_tab$Glu_sd <- fit_res$res_tab$Glu.sd /
                          fit_res$res_tab$Glu[1] * 100

p1 <- ggplot(fit_res$res_tab, aes(x = Time, y = Glu_perc_change)) + 
  geom_point() + geom_line() + ylab("Glutamate change (%)") + xlab("Time (s)") +
  annotate("rect", xmin = 3 * 60, xmax = 3 * 60 + 8 * 60, ymin = -Inf,
           ymax = Inf, alpha = 0.4, fill = "red") + ylim(c(-6, 5)) +
  geom_errorbar(aes(ymin = Glu_perc_change - Glu_sd,
                    ymax = Glu_perc_change + Glu_sd))

p2 <- ggplot(fit_res$res_tab, aes(x = Time, y = Lac_perc_change)) + 
  geom_point() + geom_line() + ylab("Lactate change (%)") + xlab("Time (s)") +
  annotate("rect", xmin = 3 * 60, xmax = 3 * 60 + 8 * 60, ymin = -Inf,
           ymax = Inf, alpha = 0.4, fill = "red") + ylim(c(-20, 40)) +
  geom_errorbar(aes(ymin = Lac_perc_change - Lac_sd,
                    ymax = Lac_perc_change + Lac_sd))

tiff(file.path("FIGURES", "S1.tiff"), width = 1500, height = 800, res = 200)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
dev.off()

# Welsh
t.test(fit_res$res_tab$Lac[3:7], fit_res$res_tab$Lac[c(1, 2, 8:15)], var.equal = FALSE)
t.test(fit_res$res_tab$Glu[3:7], fit_res$res_tab$Glu[c(1, 2, 8:15)], var.equal = FALSE)

# Student's
# n.b. reviewer 1 prefers Student's t-tests, so here they are. I'm not p-hacking!
t.test(fit_res$res_tab$Lac[3:7], fit_res$res_tab$Lac[c(1, 2, 8:15)], var.equal = TRUE)
t.test(fit_res$res_tab$Glu[3:7], fit_res$res_tab$Glu[c(1, 2, 8:15)], var.equal = TRUE)

# state      <- rep("REST", 15)
# state[3:7] <- "TASK"
# lac_tab <- data.frame(Lac = fit_res$res_tab$Lac_perc_change, state = state)
# 
# p3 <- ggplot(lac_tab, aes(x = state, y = Lac)) + geom_point() +
#   geom_signif(comparisons = list(c("REST", "TASK")), test = "t.test",
#               map_signif_level = function(p) sprintf("p = %.2g", p),
#               textsize = 3) + xlab(NULL) + ylab("Lactate change (%)")
# 
# glu_tab <- data.frame(Glu = fit_res$res_tab$Glu_perc_change, state = state)
# p4 <- ggplot(glu_tab, aes(x = state, y = Glu)) + geom_point() +
#   geom_signif(comparisons = list(c("REST", "TASK")), test = "t.test",
#               map_signif_level = function(p) sprintf("p = %.2g", p),
#               textsize = 3) + xlab(NULL) + ylab("Glutamate change (%)")
# 
# tiff(file.path("FIGURES", "FigSY.tiff"), width = 1500, height = 700, res = 200)
# plot_grid(p3, p4, labels = c('A', 'B'), label_size = 12)
# dev.off()

