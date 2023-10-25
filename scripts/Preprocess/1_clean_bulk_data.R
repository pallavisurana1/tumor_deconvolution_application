##--------------------------------------------------------------------------------------------------------------
# Initialization and Paths
##--------------------------------------------------------------------------------------------------------------

ft <- "/base/path/"
sc <- paste0(ft, "data/")
hlca_path <- paste0(ft, "projects/cell_type_origin/hlca_data/")

##--------------------------------------------------------------------------------------------------------------
# Load Required Packages
##--------------------------------------------------------------------------------------------------------------

source("packages.R")

##--------------------------------------------------------------------------------------------------------------
# Functions
##--------------------------------------------------------------------------------------------------------------

read_bulk_counts <- function(file_path) {
  counts <- fread(file_path)
  counts$geneID <- sub("\\..*", "", counts$geneID)
  counts <- counts[!duplicated(counts$geneID), ]
  return(counts %>% as.data.frame %>% remove_rownames() %>% column_to_rownames("geneID"))
}

load_mutation_data <- function(path, patterns) {
  file_list <- list.files(path, pattern = patterns, full.names = TRUE)
  mut_list <- lapply(file_list, fread)
  names(mut_list) <- sapply(file_list, function(x) tools::file_path_sans_ext(basename(x)))
  return(mut_list)
}

subset_mutation_data <- function(mut_list) {
  return(lapply(mut_list, function(df) {
    df[df[["sample_type"]] %in% c("Primary Tumor", "FFPE Scrolls"), ]
  }))
}

convert_gene_ids <- function(data, ensembl) {
  gene_id <- rownames(data)
  conversion <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "external_gene_name",
                      values = gene_id,
                      mart = ensembl)
  data <- data %>% as.data.frame %>% rownames_to_column("geneID")
  data <- merge(conversion, data, by.x = "external_gene_name", by.y = "geneID")
  data$external_gene_name <- NULL
  return(data[!duplicated(data), ] %>% as.data.frame %>% remove_rownames() %>% column_to_rownames("ensembl_gene_id"))
}

##--------------------------------------------------------------------------------------------------------------
# Main Execution
##--------------------------------------------------------------------------------------------------------------

# Read Bulk Counts
counts_raw <- read_bulk_counts("LUAD_raw/LUAD_htseq_counts_ENSG.csv")

# Load Mutation Data
mut_list <- load_mutation_data("mutation/", "*.tsv")
mut_lis <- lapply(mut_list, function(df) df[, c("sample_submitter_id", "case_submitter_id", "sample_type")])

# Subset Mutation Data
mut_lis_tumor <- subset_mutation_data(mut_lis)

# Organize Count Data based on Mutation
counts_raw_egfr <- counts_raw[, colnames(counts_raw) %in% mut_lis_tumor[["egfr"]]$sample_submitter_id]
counts_raw_kras <- counts_raw[, colnames(counts_raw) %in% mut_lis_tumor[["kras"]]$sample_submitter_id]
counts_raw_egfr_wt <- counts_raw[, !(colnames(counts_raw) %in% mut_lis_tumor[["egfr"]]$sample_submitter_id)]
counts_raw_kras_wt <- counts_raw[, !(colnames(counts_raw) %in% mut_lis_tumor[["kras"]]$sample_submitter_id)]

# Convert Gene IDs
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gtex_bulk <- convert_gene_ids(list(counts_raw_egfr, counts_raw_kras, counts_raw_egfr_wt, counts_raw_kras_wt), ensembl)

# Save Data
save(counts_raw_egfr, counts_raw_kras, counts_raw_egfr_wt, counts_raw_kras_wt, gtex_bulk, file = paste0(hlca_path, "/data/bulk_expr/", "gtex_tcga_raw_counts.RData"))


##--------------------------------------------------------------------------------------------------------------
