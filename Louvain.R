library(Seurat)
library(anndata)
library(stringr)
library(FNN)

set.seed(2025)
options(future.globals.maxSize = 1000 * 1024^2)  # Increase limit to 1000 MiB (1 GB)

plot_colors=c("0" = "#6D1A9C", "1" = "#CC79A7","2"  = "#7495D3", "3" = "#59BE86", "4" = "#56B4E9", "5" = "#FEB915", 
              "6" = "#DB4C6C", "7" = "#C798EE", "8" = "#3A84E6", "9"= "#FF0099FF", "10" = "#CCFF00FF",
              "11" = "#268785", "12"= "#FF9900FF", "13"= "#33FF00FF", "14"= "#AF5F3C", "15"= "#DAB370", 
              "16" = "#554236", "17"= "#787878", "18"= "#877F6C")

path <- normalizePath("D:\\Pitagoras\\Spatia_seq\\Dataset\\ARTICULOS\\DLPFC", winslash = "\\", mustWork = FALSE)

file_paths <-list.files(path,
                        pattern = "\\.h5ad$",
                        full.names = TRUE)
# file_path = "D:\\Pitagoras\\Spatia_seq\\Dataset\\ARTICULOS\\DLPFC/151507.h5ad"
for (file_path in file_paths) {
  print(file_path)
  
  # Read anndata objetc
  data <- read_h5ad(file_path)
  seurat <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs) # start seurat object
  # Remove spots with Nan as annotation
  cells_to_keep <- rownames(seurat@meta.data)[!is.na(seurat@meta.data$Region)]
  seurat <- subset(seurat, cells = cells_to_keep)
  
  # Initialize normalizations and PCA
  seurat <- SCTransform(seurat, verbose = F)
  seurat <- RunPCA(seurat, assay = "SCT", verbose = F, npcs = 10)
  seurat <- FindNeighbors(seurat, reduction = "pca", verbose = TRUE)

  nCluster <- length(unique(seurat$Region)) # number of clusters
  # Iterate the resolution to get the expected number of clusters
  for(qq in seq(0.05,1.5,0.01)){
    
    seurat <- FindClusters(seurat, resolution = qq,  verbose = FALSE)
    if(length(table(Idents(seurat)))==nCluster)
    {
      found_cluster <- TRUE
      break
    }
  }
  if (!found_cluster) warning("No se encontró resolución exacta para ", basename(file_path))
  
  # Obtener coordenadas espaciales desde AnnData
  coords <- data$obs[, c("array_row", "array_col")]
  coords <- coords[match(colnames(seurat), rownames(coords)), ]
  
  # Suavizar clústeres usando coordenadas de AnnData
  clusters <- as.integer(Idents(seurat))
  k <- 6
  neighbors <- get.knn(coords, k = k)$nn.index
  smoothed_clusters <- clusters
  for (i in 1:nrow(neighbors)) {
    neighbor_clusters <- clusters[neighbors[i, ]]
    smoothed_clusters[i] <- as.integer(names(sort(table(neighbor_clusters), decreasing = TRUE)[1]))
  }
  
  # Add to the former anndata the classification
  data$obs$Louvain <- seurat@meta.data$seurat_clusters[match(rownames(data$obs), colnames(seurat))]
  data$obs$`Louvain smooth` <- smoothed_clusters[match(rownames(data$obs), colnames(seurat))]
  
  # Make new subfolder: Louvain
  louvain_path <- file.path(dirname(file_path), "Louvain")
  dir.create(louvain_path, recursive = T)
  save_path <- file.path(louvain_path, basename(file_path))
  
  # Save the new object in a subfolder: Louvain
  write_h5ad(data, save_path)
}
