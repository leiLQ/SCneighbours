# SCneighbours
Single-Cell Neighbourhood Exploration Tool 

## Technique 1: Visualization of Cell Neighbourhoods

For individual cells or cell groups, we visualize their cell neighbourhoods as either dots or density contours on a dimensionality reduction map, uncovering tightly interwoven or widely dispersed regions. This allows us to learn which clusters might need merging if they have largely overlapped neighbourhoods, and infer trajectory relationships of cells obscured by the reduction.
```
source("visualize_neighbourhood_function.R")

visualize_neighbourhood(seu, meta_data_column = "seurat_clusters", meta_data_highlight = "12", reduction = "umap", density = T)
```

## Technique 2: Quantification of Cell Neighbourhoods

For each cell, we quantify the distribution of its neighbouring cells by the variance in their coordinates, providing insights into the cell positioning in the reduction map. Large variance indicates diminished confidence in a cell’s positioning.

```
source("quantify_neighbourhood_function.R")

calculate_neighbour_distance_for_all_cells(seu, reduction = "umap", colname = "neighbourhood_variance")
```

## Technique 3: Complementary reductions

The cell neighbourhoods can be explored across various dimensionality reductions for insights. For example, diffusion map and force-directed layout can provide clues about cell development by preserving the continuous cell transitioning trajectories.

