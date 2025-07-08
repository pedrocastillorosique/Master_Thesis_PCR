library("Seurat")
library("anndata")
library('ggplot2')
library('Matrix')
library('ggrepel')
library('igraph')

plot_colors=c("0" = "#6D1A9C", "1" = "#CC79A7","2"  = "#7495D3", "3" = "#59BE86", "4" = "#56B4E9", "5" = "#FEB915", 
              "6" = "#DB4C6C", "7" = "#C798EE", "8" = "#3A84E6", "9"= "#FF0099FF", "10" = "#CCFF00FF",
              "11" = "#268785", "12"= "#FF9900FF", "13"= "#33FF00FF", "14"= "#AF5F3C", "15"= "#DAB370", 
              "16" = "#554236", "17"= "#787878", "18"= "#877F6C")

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(rsvd))
options(future.globals.maxSize = 1000 * 1024^2)  # Increase limit to 1000 MiB (1 GB)

Cell_modules <- function(basePath, filePath, nCluster = 7, save_path = NULL, pdf_file = NULL ){
  MP_rep         = as.matrix(read.table( paste0(basePath, "Semantic_representations.txt"), header = T, row.names = 1))
  SC_rep         = as.matrix(read.table( paste0(basePath, "Hierarchical_representations.txt"), header = T, row.names = 1))
  robust_rep = cbind(MP_rep, SC_rep)
  
  data = read_h5ad(filePath)
  
  seurat <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs)
  cells_to_keep <- rownames(seurat@meta.data)[!is.na(seurat@meta.data$Region)]
  seurat <- subset(seurat, cells = cells_to_keep)

  idc = SCTransform(seurat, verbose = TRUE)
  idc = RunPCA(idc, assay = "SCT", verbose = FALSE, npcs = 100)

  inter_c = intersect(row.names(robust_rep), colnames(idc))
  Cell_obj = subset(idc, cells =inter_c)      
  in_feas = robust_rep[match(colnames(Cell_obj), row.names(robust_rep)),]
  Cell_obj@reductions$pca@cell.embeddings[,1:dim(in_feas)[2]] = in_feas
  Cell_obj = FindNeighbors(Cell_obj, reduction = "pca", dims = 1:dim(in_feas)[2])

  nCluster <- length(unique(seurat$Region))

  for(qq in seq(0.05,1.5,0.01))
  {
    Cell_obj = FindClusters( Cell_obj, resolution = qq,  verbose = FALSE )
    if(length(table(Idents(Cell_obj)))==nCluster)
    {
      break
    }
  }
  data$obs$stKeep <- Cell_obj@meta.data$seurat_clusters[match(rownames(data$obs), colnames(Cell_obj))]
  save_path <- file.path(basePath, basename(filePath))
  
  write_h5ad(data, save_path)
}