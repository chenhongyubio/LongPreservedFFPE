library(Seurat)
library(paletteer) 
library(ggplot2)
library(ggsci)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
library(Scillus)
library(patchwork)
library(rliger)

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

### liger
scRNA_liger <- merge(seurat_list[[1]], y=c(unlisted(seurat_list[-1])))
scRNA_liger <- NormalizeData(scRNA_liger) %>% FindVariableFeatures() %>% ScaleData(split.by = "orig.ident", do.center = FALSE)
scRNA_liger <- RunOptimizeALS(scRNA_liger, k = 20, lambda = 5, split.by = "orig.ident")
scRNA_liger <- RunQuantileNorm(scRNA_liger, split.by = "orig.ident")
scRNA_liger <- FindNeighbors(scRNA_liger, reduction = "iNMF", dims = 1:20) %>% FindClusters(resolution = 0.3)
scRNA_liger <- RunUMAP(scRNA_liger, dims = 1:ncol(scRNA_liger[["iNMF"]]), reduction = "iNMF")
saveRDS(scRNA_liger,paste0(data_path,"/obj_liger_unanno.rds"))

# vlnplot(Refer to Figure 1B)
scRNA_liger$sample <- scRNA_liger$orig.ident
scRNA_liger$sample <- gsub("A524","PP2_1",scRNA_liger$sample)
scRNA_liger$sample <- gsub("A525","PP2_2",scRNA_liger$sample)
scRNA_liger$sample <- gsub("A526","PP2_3",scRNA_liger$sample)
scRNA_liger$sample <- gsub("A527","PP2_4",scRNA_liger$sample)
scRNA_liger$sample <- gsub("A633","PP1_1",scRNA_liger$sample)
scRNA_liger$sample <- gsub("A634","PP1_2",scRNA_liger$sample)
scRNA_liger$sample <- gsub("A635","PP1_3",scRNA_liger$sample)
scRNA_liger$sample <- gsub("A636","PP1_4",scRNA_liger$sample)
scRNA_liger$sample <- gsub("A818","GP1_4",scRNA_liger$sample)
scRNA_liger$sample <- gsub("K_A630","GP1_1",scRNA_liger$sample)
scRNA_liger$sample <- gsub("K_A631","GP1_2",scRNA_liger$sample)
scRNA_liger$sample <- gsub("K_A632","GP1_3",scRNA_liger$sample)
scRNA_liger$sample <- gsub("K_A1047L","PP3_1",scRNA_liger$sample)
scRNA_liger$sample <- gsub("K_A1048L","PP3_2",scRNA_liger$sample)
scRNA_liger$sample <- gsub("K_A518","GP2_1",scRNA_liger$sample)
scRNA_liger$sample <- gsub("K_A519","GP2_2",scRNA_liger$sample)
scRNA_liger$sample <- gsub("K_A520","GP2_3",scRNA_liger$sample)
scRNA_liger$sample <- gsub("K_A521","GP2_4",scRNA_liger$sample)

scRNA_liger[["percent.mt"]] <- PercentageFeatureSet(scRNA_liger, pattern = "^MT-")

png("./BaseInfo1.png",600,300,res=120)
VlnPlot(scRNA_liger, features = c("nFeature_RNA"), group.by = "sample", pt.size = 0, y.max = 2000) +NoLegend()
dev.off()

png("./BaseInfo2.png",600,300,res=120)
VlnPlot(scRNA_liger, features = c("nCount_RNA"), group.by = "sample", pt.size = 0, y.max = 2000) +NoLegend()
dev.off()

png("./BaseInfo3.png",600,300,res=120)
VlnPlot(scRNA_liger, features = c("percent.mt"), group.by = "sample", pt.size = 0, y.max = 100) +NoLegend()
dev.off()

## marker anno
markers <- c("KRT7","HNF1B","CFTR", # Cholangiocytes
             "VWF","CD93","EMCN",  # Endothelial cells
             "ALB","SERPINA6","NR1I3", # Hepatocytes
             "CD163","TIMD4","VSIG4",  # Kuppffer cells
             "COL1A2","COL3A1","SPARC", # Stellate_cells
             "PIM2","MZB1","IGHG1",  # Plasma cells
             "CD3E","IL7R","LEF1") # T cells
DotPlot(scRNA_liger,features = markers)
FeaturePlot(scRNA_liger,features = markers)

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

names(new.cluster.ids) <- levels(scRNA_liger)
scRNA_liger <- RenameIdents(scRNA_liger, new.cluster.ids)
scRNA_liger$celltype <- Idents(scRNA_liger)

saveRDS(scRNA_liger,"./obj_all_liger_anno.rds")