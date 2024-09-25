library(Seurat)
library(paletteer) 
library(ggplot2)
library(ggsci)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
library(Scillus)
library(patchwork)
library(SeuratWrappers)

data_path <- "./Tumor_heterogeneity/raw_data/Matrix"

folders <- list.dirs(path = data_path, recursive = FALSE)
seurat_list <- list()
for(i in 1:length(folders)) {
      data <- Read10X(data.dir = folders[i])
      sample_name <- basename(folders[i])
      obj <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3, min.features = 200)
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
      obj <- subset(obj, subset = percent.mt < 25)
      seurat_list[[i]] <- obj
      names(seurat_list)[i] <- sample_name
}

result_obj <- merge(seurat_list[[1]], y=c(unlisted(seurat_list[-1])))
result_obj <- NormalizeData(result_obj) %>% FindVariableFeatures() %>% ScaleData(split.by = "orig.ident", do.center = FALSE)
result_obj <- RunOptimizeALS(result_obj, k = 20, lambda = 5, split.by = "orig.ident")
result_obj <- RunQuantileNorm(result_obj, split.by = "orig.ident")

split_obj <- SplitObject(result_obj, split.by = "sample")
split_obj <- lapply(split_obj, function(x) {
  x <- FindNeighbors(x, reduction = "iNMF", dims = 1:20) %>%
       FindClusters(resolution = 0.3)
  x <- RunUMAP(x, dims = 1:ncol(x[["iNMF"]]), reduction = "iNMF")
  return(x)
})

## marker anno
markers <- c("KRT7","HNF1B","CFTR", # Cholangiocytes
             "VWF","CD93","EMCN",  # Endothelial cells
             "ALB","SERPINA6","NR1I3", # Hepatocytes
             "CD163","TIMD4","VSIG4",  # Kuppffer cells
             "COL1A2","COL3A1","SPARC", # Stellate_cells
             "PIM2","MZB1","IGHG1",  # Plasma cells
             "CD3E","IL7R","LEF1") # T cells
DotPlot(split_obj[[1]],features = markers)
FeaturePlot(split_obj[[1]],features = markers)

new.cluster.ids <- c("Hepatocytes",
"T cells",
"T cells",
"Kupffer cells",
"Stellate cells",
"Hepatocytes",
"T cells",
"T cells",
"T cells",
"Endothelial cells",
"T cells",
"B/Plasma cells",
"Epithelial cells",
"Hepatocytes")

names(new.cluster.ids) <- levels(split_obj[[1]])
split_obj[[1]] <- RenameIdents(split_obj[[1]], new.cluster.ids)
split_obj[[1]]$celltype <- Idents(split_obj[[1]])

## split_umap(Refer Figure 1C)
pdf("./PaperResult/UMAP_split_sample_finial.pdf",13,8)
umap_plots <- lapply(seq_along(split_obj), function(i) {
  current_name <- names(split_obj)[i]
  CellDimPlot(srt = split_obj[[i]], 
              group.by = "celltype", 
              reduction = "UMAP", 
              theme_use = "theme_blank", 
              legend.position = "none",
              subtitle = current_name)
})
combined_plot <- wrap_plots(umap_plots, ncol = 6)
print(combined_plot)
dev.off()
