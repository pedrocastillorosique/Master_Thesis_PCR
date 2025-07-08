library(Seurat)
library(anndata)
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(stringr)

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


for (file_path in file_paths) {
  print(file_path)
  
  data <- read_h5ad(file_path)
  
  # Convert to Seurat
  seurat <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs)
  cells_to_keep <- rownames(seurat@meta.data)[!is.na(seurat@meta.data$Region)]
  seurat <- subset(seurat, cells = cells_to_keep)

  # Normalize and apply PCA
  seurat <- SCTransform(seurat, verbose = F)
  seurat <- RunPCA(seurat, assay = "SCT", verbose = F, npcs = 50)
  
  #Convert to SCE object
  sce = as.SingleCellExperiment(seurat, assay = "SCT")
  coords = cbind(seurat@meta.data$array_col, seurat@meta.data$array_row) #convert seurat to SCE
  colnames(coords) = c("array_col", "array_row")
  colData(sce) = cbind(colData(sce), coords) #add spatial info to SCE
  
  #BayesSpace Workflow
  sce = spatialPreprocess(sce, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
  sce = spatialCluster(sce, nrep = 5000, burn.in = 500, 
                       q = length(unique(seurat$Region)), 
                       d = 25, gamma = 3, init.method = "kmeans")
  
  seurat@meta.data = cbind(seurat@meta.data, BayesSpace = sce$spatial.cluster) #add BayesSpace clusters to Seurat obj
  
  data$obs$BayesSpace <- sce$spatial.cluster[match(rownames(data$obs), colnames(sce))]
  
  bayes_path <- file.path(dirname(file_path), "BayesSpace")
  dir.create(bayes_path, recursive = TRUE)
  save_path <- file.path(bayes_path, basename(file_path))
  write_h5ad(data, save_path)
}
