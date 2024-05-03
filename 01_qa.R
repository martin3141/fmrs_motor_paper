library(spant)

source("extra_functions.R")

paths  <- Sys.glob(file.path("DATA", "sub-??", "mrs", "sub-??_svs.nii.gz"))
labels <- strtrim(basename(paths), 6)

# all data
preproc_all <- preproc_fmrs_dataset_paper(paths, labels = labels,
                                          output_dir = "all_fmrs_preproc")

# remove sub-02 due to drift
paths <- grep("sub-02", paths, invert = TRUE, value = TRUE)
# remove sub-03 due to movement
paths <- grep("sub-03", paths, invert = TRUE, value = TRUE)
# remove sub-10 due to movement
paths <- grep("sub-10", paths, invert = TRUE, value = TRUE)
# remove sub-23 due to movement
paths <- grep("sub-23", paths, invert = TRUE, value = TRUE)

labels <- strtrim(basename(paths), 6)

# only good data
preproc_good <- preproc_fmrs_dataset_paper(paths, labels = labels,
                                           output_dir = "good_fmrs_preproc")

tiff("mean_spec.tiff", res = 200, width = 1000, height = 1000)
preproc_good$res_list[[17]]$mean_corr |> zf() |> plot(xlim = c(4, 0.5))
dev.off()
