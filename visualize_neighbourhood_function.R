library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)

# plot neighbour cells of an individual cell or a certain cell group (i.e., cluster)

#'@title visualize_neighbourhood.
#'@description The function to highlight the neighbour cells of certain cells.
#'@param meta_data_column Name of the column in seurat_object@meta.data slot to pull value from for highlighting.
#'@param meta_data_highlight Name of variable(s) within meta_data_name to highlight in the plot.
#'@param reduction The reduction map used to calculate the coordinate variance.
#'@param seu The seurat object.
#'@param density Whether to plot the density countour plot of neighbour cells. Default is FALSE.
#'@return A number, the average variance in coordinates of neighbour cells of a certain cell i based on the provided reduction map.
visualize_neighbourhood = function(seu, meta_data_column, meta_data_highlight, reduction, density = F, graph = "RNA_nn") {
  # all neighbour cells 
  n = colnames(seu@graphs[[graph]])[colSums(seu@graphs[[graph]][seu[[meta_data_column]] == meta_data_highlight,]) > 0] # please change "RNA_nn" to "SCT_nn" if you are using SCTranform
  
  # reduction emeddings of these neighbour cells
  Embeddings(seu, reduction = reduction) %>%
    as_tibble(rownames = "bc") %>%
    mutate(neighbour = bc %in% n) -> d
  
  # matching axis names to the reduction names, please adjust it according to the real case
  axis = seu@reductions[[reduction]]@key
  #axis <- gsub("umap","UMAP",gsub("dm", "DC", reduction))           
  axis_1 <- rlang::sym(paste0(axis, "1")) 
  axis_2 <- rlang::sym(paste0(axis, "2")) 
  
  if (!density){
    # plot by just highlighting these neighbour cells in the reduction map
    ggplot(d, aes(!!axis_1, !!axis_2)) +
      geom_point(colour = "grey") +
      geom_point(data = d[d$neighbour,], colour = "red", size = 0.6) +
      theme_classic()+ ggtitle(paste0(meta_data_column, meta_data_highlight))
    
  } else{
    # plot the density plot of these neighbour cells in the reduction map
    ggplot(d, aes(!!axis_1, !!axis_2)) +
      geom_point(colour = "grey") +
      geom_density_2d(data = d[d$neighbour,]) +
      theme_classic() + ggtitle(paste0(meta_data_column, meta_data_highlight))
  }
}
