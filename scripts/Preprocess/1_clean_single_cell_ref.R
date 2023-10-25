##--------------------------------------------------------------------------------------------------------------
# Load Required Packages
##--------------------------------------------------------------------------------------------------------------

source("packages.R")

##-----------------------------------------------------------------------------------------------------------------
# Functions
##-----------------------------------------------------------------------------------------------------------------

load_data <- function(path){
  return(readRDS(path))
}

filter_healthy_cells <- function(data){
  return(subset(data, subset=condition %in% c("healthy", "nan")))
}

calculate_mito_genes <- function(data, genes, assay="RNA"){
  Percent_mito_genes <- PercentageFeatureSet(data, features = genes, assay=assay)
  data@meta.data$Percent_mito_genes <- Percent_mito_genes
  return(data)
}

apply_filters <- function(data){
  return(subset(data, 
                subset = Percent_mito_genes < 10 & 
                  nFeature_RNA > 200 & 
                  nFeature_RNA < 2500 & 
                  !anatomical_region_level_3 %in% c("2_trachea", "2_inferior_turbinate", 
                                                    "inferior turbinate right nostril", 
                                                    "inferior turbinate left nostril") & 
                  as.numeric(age) >= 25.0 & 
                  !cell_type %in% c("nasal mucosa goblet cell", 
                                    "tracheobronchial goblet cell", 
                                    "lung neuroendocrine cell", 
                                    "serous secreting cell", 
                                    "brush cell of trachebronchial tree", 
                                    "fibroblast", 
                                    "mucus secreting cell", 
                                    "stromal cell")))
}

create_SCE <- function(data){
  counts <- as.matrix(GetAssayData(object=data, slot="counts"))
  colData <- as.data.frame(matrix(nrow=nrow(counts), ncol=2))
  rownames(colData) <- colnames(data)
  colnames(colData) <- c("clusters", "samples")
  colData[,1] <- data@meta.data$cell_type
  colData[,2] <- data@meta.data$sample
  
  return(SingleCellExperiment(list(counts=counts), colData=colData))
}

##-----------------------------------------------------------------------------------------------------------------
# Main script
##-----------------------------------------------------------------------------------------------------------------

# Paths (Users can modify these)
input_path <- "Human_Lung_Cell_Atlas_cellxgene_local.rds"
output_path_dataset <- "lung_singlecell_healthy.rds"
output_path_sce <- "lung_singlecell_healthy_sce.rds"

# Mitochondrial genes
mito_genes <- c("ENSG00000198888", "ENSG00000198763", "ENSG00000198840", 
                "ENSG00000212907", "ENSG00000198886", "ENSG00000198786",
                "ENSG00000198695", "ENSG00000198727", "ENSG00000198804", 
                "ENSG00000198712", "ENSG00000198938", "ENSG00000198899", 
                "ENSG00000228253")

# Load dataset
lung_singlecell_data <- load_data(input_path)

# Process data
lung_singlecell_healthy <- lung_singlecell_data %>% 
  filter_healthy_cells() %>% 
  calculate_mito_genes(mito_genes) %>% 
  apply_filters()

# Save processed dataset
saveRDS(lung_singlecell_healthy, file=output_path_dataset)

# Create and save SingleCellExperiment object
lung_sce <- create_SCE(lung_singlecell_healthy)
saveRDS(lung_sce, file=output_path_sce)

##-----------------------------------------------------------------------------------------------------------------