library(spant)

paths  <- Sys.glob(file.path("DATA", "sub-??", "mrs", "sub-??_svs.nii.gz"))

# remove sub-02 due to drift
paths <- grep("sub-02", paths, invert = TRUE, value = TRUE)
# remove sub-03 due to movement
paths <- grep("sub-03", paths, invert = TRUE, value = TRUE)
# remove sub-10 due to movement
paths <- grep("sub-10", paths, invert = TRUE, value = TRUE)
# remove sub-23 due to movement
paths <- grep("sub-23", paths, invert = TRUE, value = TRUE)

# remove sub-14 due to no acquired fMRI
paths <- grep("sub-14", paths, invert = TRUE, value = TRUE)

labels <- strtrim(basename(paths), 6)

t1_paths <- file.path("DATA", labels, "func", paste0(labels, "_bold.feat"),
                      "hr", "background.nii.gz")

voxel_paths <- file.path("DATA", labels, "func", paste0(labels, "_bold.feat"),
                         "hr", "svs_voi.nii.gz")

for (n in 1:length(voxel_paths)) {
  svs_voi <- get_svs_voi(read_mrs(paths[n]), readNifti(t1_paths[n]))
  writeNifti(svs_voi, voxel_paths[n])
}
