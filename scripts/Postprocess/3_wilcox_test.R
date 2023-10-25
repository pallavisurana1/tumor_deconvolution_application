##--------------------------------------------------------------------------------------------------------------
# Load Required Packages
##--------------------------------------------------------------------------------------------------------------

source("Preprocess/packages.R")

##--------------------------------------------------------------------------------------------------------------
# Read deconvolution results
##--------------------------------------------------------------------------------------------------------------

read_results <- function(hlca_path) {
  names <- list.files(paste0(hlca_path, "/results/props_instaprism"), pattern = "*.csv", recursive = TRUE, full.names = TRUE)
  res <- lapply(names, fread)
  names <- basename(names)
  names <- sub(".csv", "", names)
  names(res) = names
  return(res)
}

##--------------------------------------------------------------------------------------------------------------
# Modify results
##--------------------------------------------------------------------------------------------------------------

modify_data <- function(res) {
  prop.est.music <- prop <- res
  for (i in 1:length(prop.est.music)) {
    prop.est.music[[i]] <- res[[i]] %>% as.data.frame() %>% column_to_rownames("cell_types") %>% t() %>% as.data.frame
    prop.est.music[[i]] <- prop.est.music[[i]] %>% rownames_to_column("sampleID")
    prop[[i]] <- prop.est.music[[i]] %>% melt() %>% as.data.frame()
  }
  
  nm <- names(prop)
  for(i in 1:length(prop)) {
    prop[[i]]$type <- nm[i]
    prop.est.music[[i]]$type <- nm[i]
  }
  return(list(prop_HLCA = prop.est.music))
}

##--------------------------------------------------------------------------------------------------------------
# Data viz before applying test
##--------------------------------------------------------------------------------------------------------------

plot_distribution <- function(df) {
  plt <- ggplot(df, aes(bcell)) + 
    geom_histogram(bins = 10) + 
    facet_wrap(~type, scales = 'free_x') +
    ggtitle('Distribution of proportions')
  print(plt)
}

##--------------------------------------------------------------------------------------------------------------
# Wilcoxon test applied to our data
##--------------------------------------------------------------------------------------------------------------

perform_wilcox_test <- function(list_plot) {
  iter <- length(list_plot)
  stat.test <- vector(mode = "list", length = iter)
  eff_size <- vector(mode = "list", length = iter)
  for (i in 1:iter) {
    stat.test[[i]] <- list_plot[[i]] %>% wilcox_test(value ~ type) %>% add_significance()
    eff_size[[i]] <- list_plot[[i]] %>% wilcox_effsize(value ~ type)
  }
  stat <- do.call(rbind.data.frame, stat.test)
  eff <- do.call(rbind.data.frame, eff_size)
  merge_results <- merge(stat, eff, by = "tissue")
  return(merge_results)
}


run_analysis <- function() {
  load_libraries()
  
  paths <- get_paths()
  results <- read_results(paths$hlca_path)
  modified_data <- modify_data(results)
  
  plot_distribution(modified_data$prop_HLCA)
  
  melted.df <- plyr::ldply(modified_data$prop_HLCA)
  list.plot <- split(melted.df, f = melted.df$variable)
  
  wilcox_results <- perform_wilcox_test(list.plot)
  return(wilcox_results)
}

##--------------------------------------------------------------------------------------------------------------
# Call the main function
result <- run_analysis()
print(result)

# save if needed ......
##--------------------------------------------------------------------------------------------------------------