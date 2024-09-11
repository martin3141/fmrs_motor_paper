library(spant)
library(cowplot)

theme_set(theme_bw())

# only good data
preproc <- readRDS(file.path("good_fmrs_preproc", "preproc_corrected.rds"))

block_size <- 50

mean_proc <- preproc$mean_dataset |> mean_dyn_blocks(block_size)

mol_names <- get_1h_brain_basis_names(add = c("gly", "peth"))

basis     <- sim_basis(mol_names, pul_seq = seq_slaser_ideal,
                       acq_paras = mean_proc)

first_dyn <- get_dyns(mean_proc, 1)

fit_res_abfit <- fit_mrs(first_dyn, basis = basis)

fit_res_lcm   <- fit_mrs(first_dyn, basis = basis, method = "LCMODEL",
                         opts = c("NSIMUL=0"))

tiff(file.path("FIGURES", "FigS1.tiff"), width = 1500, height = 700, res = 200)
plot_grid(~plot(fit_res_abfit), ~plot(fit_res_lcm), labels = c('A', 'B'),
          label_size = 12)
dev.off()
