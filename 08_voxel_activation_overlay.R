library(spant)
library(RNifti)
library(fslr)
library(cowplot)

options(fsl.path = "/Users/martinwilson/fsl")

# path to MNI152 template image
standard_path <- file.path("FEAT_HL_FIXED_EFFECTS.gfeat", "bg_image.nii.gz")
standard      <- readNifti(standard_path)

# paths to feat output
# (14 is missing due to participant having mid-scan nap)
feat_paths    <- Sys.glob(file.path("DATA", "*", "func", 
                                    "*bold_std_space.feat"))
highres_paths <- file.path(feat_paths, "reg", "highres.nii.gz")
mrs_paths     <- Sys.glob(file.path(dirname(dirname(feat_paths)), "mrs",
                                    "*svs.nii.gz"))
hr2std_paths  <- file.path(feat_paths, "reg", "highres2standard.mat")

n_voxels <- length(feat_paths)
for (n in 1:n_voxels) {
  mrs     <- read_mrs(mrs_paths[n])
  highres <- readNifti(highres_paths[n])
  voxel   <- get_svs_voi(mrs, highres)
  
  voxel_std <- flirt_apply(infile = voxel, reffile = standard_path,
                           initmat = hr2std_paths[n])
  
  if (n == 1) {
    mean_voxel_std <- voxel_std / n_voxels
  } else {
    mean_voxel_std <- mean_voxel_std + voxel_std / n_voxels
  }
}

orientation(standard)       <- "RAS"
orientation(mean_voxel_std) <- "RAS"

p1 <- \(x) ortho3(standard, mean_voxel_std, zlim_ol = c(0, 1), xyz = c(28, 52, 64))

stat <- readNifti(file.path("FEAT_HL_FIXED_EFFECTS.gfeat", "cope1.feat",
                            "thresh_zstat1.nii.gz"))

orientation(stat) <- "RAS"


p2 <- \(x) ortho3(standard, stat, zlim_ol = c(14, 32), xyz = c(28, 52, 64))


tiff(file.path("FIGURES", "SX.tiff"), width = 1500, height = 800, res = 200)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12, scale = 0.81,
          label_colour = "white")
dev.off()
