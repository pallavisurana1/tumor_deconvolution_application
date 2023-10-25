##--------------------------------------------------------------------------------------------------------------
# Load Required Packages
##--------------------------------------------------------------------------------------------------------------

source("Preprocess/packages.R")

##--------------------------------------------------------------------------------------------------------------
## Load Data Functions
##--------------------------------------------------------------------------------------------------------------

load_bulk_data <- function(hlca_path) {
  readRDS(paste0(hlca_path, "/data/bulk_expr/", "gtex_tcga_filtered_counts_for_deconv.rds"))
}

load_single_cell_ref_data <- function(hlca_path) {
  hlca.sce = readRDS(paste0(hlca_path, "/data/sc_expr/lung_singlecell_healthy_lungcells_sce.rds"))
  list(data = counts(hlca.sce), meta = hlca.sce@colData)
}

##--------------------------------------------------------------------------------------------------------------
## InstaPrism 
##--------------------------------------------------------------------------------------------------------------

run_instaprism <- function(sc_data, bk_data_elem, cell_type_labels) {
  myPrism <- new.prism(
    reference = t(sc_data), 
    mixture = t(bk_data_elem),
    input.type = "count.matrix", 
    cell.type.labels = cell_type_labels, 
    cell.state.labels = NULL,
    key = NULL,
    outlier.cut = 0.01,
    outlier.fraction = 0.1
  )
  
  ip.res = InstaPrism(input_type = 'prism', prismObj = myPrism, update = FALSE, return.Z.cs = FALSE, return.Z.ct = TRUE)
  (ip.res@Post.ini.ct@theta) %>% as.data.frame %>% rownames_to_column("cell_types")
}

##--------------------------------------------------------------------------------------------------------------
## MUSIC 
##--------------------------------------------------------------------------------------------------------------

music_deconvolution <- function(bulk_data_elem, hlca_sce, fname_i) {
  est.prop = music_prop(
    bulk.mtx = bulk_data_elem, 
    sc.sce = hlca_sce,
    clusters = 'clusters', 
    samples = 'samples', 
    select.ct = NULL,
    verbose = TRUE
  )
  est.prop$Est.prop.weighted
}

##--------------------------------------------------------------------------------------------------------------
## Main Process
##--------------------------------------------------------------------------------------------------------------

hlca_path <- "YOUR_PATH_HERE"  # Replace with your actual path

# Load Data
bk.data <- load_bulk_data(hlca_path)
cle <- load_cell_line_data(hlca_path)
sc <- load_single_cell_ref_data(hlca_path)

# InstaPrism Deconvolution
for(i in 1:length(bk.data)) {
  props <- run_instaprism(sc$data, bk.data[[i]], as.vector(sc$meta$clusters))
  fwrite(props, paste0(hlca_path, "/results/props_instaprism/", names(bk.data)[i], ".csv"))
}


# MUSIC Deconvolution
for(i in 1:5) {
  est.prop <- music_deconvolution(bk.data[[i]], sc$hlca.sce, fname[i])
  fwrite(est.prop, paste0(hlca_path, "/results/props_music/", fname[i], ".csv"))
}

##--------------------------------------------------------------------------------------------------------------
