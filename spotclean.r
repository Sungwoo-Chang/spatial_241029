# install packages
# devtools::install_github("zijianni/SpotClean", build_manual = TRUE, build_vignettes = TRUE)

# Load packages
library(SpotClean)
library(S4Vectors)

library(SpatialExperiment)

# loac count matrix ============================================================
data(mbrain_raw)
str(mbrain_raw)

# read spatial metadata ========================================================
folderPath <- "/home/jsw/BI/Database/Ravi_Cancer_Cell_2022/10XVisium_2/#UKF260_T_ST/outs/"

count_matrix <- read10xRaw(file.path(folderPath, "raw_feature_bc_matrix"))
# count_matrix <- read10xRaw(file.path(folderPath, 'filtered_feature_bc_matrix'))

spe <- read10xSlide(
  tissue_csv_file = file.path(folderPath, "spatial/tissue_positions_list.csv"),
  tissue_img_file = file.path(folderPath, "spatial/tissue_lowres_image.png"),
  scale_factor_file = file.path(folderPath, "spatial/scalefactors_json.json")
)

str(spe)

# create the slide object ======================================================
slide_obj <- createSlide(
  count_mat = count_matrix,
  slide_info = spe
)
slide_obj

# vsualize the slide ===========================================================
x11(width = 5, height = 5)
visualizeSlide(slide_obj)

visualizeLabel(slide_obj, "tissue")

metadata(slide_obj)$slide$total_counts <- Matrix::colSums(count_matrix)
visualizeHeatmap(slide_obj, "total_counts")
visualizeHeatmap(slide_obj, rownames(count_matrix[1]))

visualizeHeatmap(slide_obj, "OLIG2")

decont_obj <- spotclean(slide_obj, maxit = 10, candidata_radius = 20)

names(metadata(decont_obj))
x11(width = 5, height = 5)
visualizeHeatmap(decont_obj, "OLIG2")
