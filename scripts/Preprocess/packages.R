
##--------------------------------------------------------------------------------------------------------------
# Installing and loading required libraries
##--------------------------------------------------------------------------------------------------------------

# Load the pacman library
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)


# CRAN Packages
cran_packages <- c(
  "dplyr", "tidyr", "tidyverse", "tibble", "data.table", "smplot", "ggplot2",
  "openxlsx", "reshape", "ggpubr", "cowplot", "ggstatsplot", "ggsignif",
  "rstatix", "stringr"
)

# Install CRAN packages if they are not already installed
p_load(char = cran_packages)



# Bioconductor Packages
bioc_packages <- c("Biobase", "BisqueRNA", "Seurat", "SingleCellExperiment")

# Install Bioconductor packages if they are not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(bioc_packages)

# Load Bioconductor packages
sapply(bioc_packages, requireNamespace, quietly = TRUE)



# GitHub Packages (For illustration, using assumed paths. Replace with actual paths if these are incorrect)
github_packages <- c("humengying0907/InstaPrism", "Danko-Lab/BayesPrism", "xuranw/MuSiC")

# Install GitHub packages if they are not already installed and load them
p_load_gh(github_packages)
