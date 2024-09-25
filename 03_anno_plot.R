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
library(SCP)
library(Scillus)

setwd("/public/home/chenhy/Tumor_heterogeneity/Result")

result_obj <- readRDS("./obj_all_liger_unanno.rds")

# UMAP plot(Refer to Figure 1C)
pdf("./PaperResult/UMAP.pdf",7,6)
CellDimPlot(
  srt = result_obj, group.by = c("celltype"),
  reduction = "UMAP", theme_use = "theme_blank"
)
dev.off()

# plot(Refer to Figure 3A)
pdf("./PaperResult/BILI.pdf",10,7)
CellStatPlot(srt = result_obj, stat.by = "celltype", group.by = "sample")
dev.off()

# meta plot(Refer to Figure 3B)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
clinical.df=data.frame(
  Sample=c(rep("Patient1",4),rep("Patient2",4),rep("Patient3",4),rep("Patient4",4),rep("Patient5",2)),
  Sample2=c(paste0("GP1_",1:4),paste0("GP2_",1:4),paste0("PP1_",1:4),paste0("PP2_",1:4),paste0("PP3_",1:2)),
  Type=c("Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Normal","Tumor","Tumor","Tumor","Normal","Tumor","Tumor","Tumor","Normal","Tumor","Tumor"),
  Prognosis=c(rep("Good",8),rep("Poor",10))
)

clinical.df2=melt(clinical.df,id="Sample2")
clinical.df2$variable=factor(clinical.df2$variable,levels = c("Sample","Type","Prognosis"))

cols=c(
  "Patient1"="#14517C","Patient2"="#2F7FC1","Patient3"="#E7EFFA","Patient4"="#96C37D","Patient5"="#F3D266",
  "Tumor"="#66C2A5","Normal"="#FC8D62",
  "Good"="#377EB8","Poor"="#E41A1C"
)

pdf("./PaperResult/sample_info.pdf",10,3)
clinical.df2%>%ggplot(aes(x=Sample2,y=variable))+
  geom_tile(aes(fill=value),color="white",size=1)+ 
  scale_x_discrete("",expand = c(0,0))+ 
  scale_y_discrete("",expand = c(0,0))+
  scale_fill_manual(values = cols)+ 
  theme(
    axis.text.x.bottom = element_text(size=10,color = "black"),axis.text.y.left = element_text(size = 12,color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_blank(), 
    legend.title = element_blank() 
)
dev.off()

## CelltypePlot(Refer to Figure2)
biotype <- read.table("/public/home/chenhy/Tumor_heterogeneity/human_gene_biotype.txt", header = F, sep = "\t")
colnames(biotype) <- c("gene","biotype")

keep_biotype <- c("protein_coding","lncRNA","miRNA","snRNA","snoRNA","scaRNA","misc_RNA")

raw_data <- Read10X("./PP3_2")
obj <- CreateSeuratObject(counts = raw_data,  min.cell = 1)
result <- matrix(ncol=2, nrow = 8)
colnames(result) <- c("biotype","PP3_2")
gene_list <- list()
n <- 1
for(i in keep_biotype){
    gene <- biotype[biotype$biotype == i,1]
    count <- sum(rownames(obj) %in% gene)
    result[n,1] <- i
    result[n,2] <- count
    keep_gene <- rownames(obj)[rownames(obj) %in% gene]
    gene_list[[i]] <- keep_gene
    n <- n+1
}
result[n,1] <- "others"
result[n,2] <- dim(obj)[1]-sum(as.numeric(result[1:7,2]))
PP3_2 <- result
PP3_2_list <- gene_list

all <- merge(GP1_1, GP1_2, by="biotype")
all <- merge(all, GP1_3, by="biotype")
all <- merge(all, GP1_4, by="biotype")
all <- merge(all, GP2_1, by="biotype")
all <- merge(all, GP2_2, by="biotype")
all <- merge(all, GP2_3, by="biotype")
all <- merge(all, GP2_4, by="biotype")
all <- merge(all, PP1_1, by="biotype")
all <- merge(all, PP1_2, by="biotype")
all <- merge(all, PP1_3, by="biotype")
all <- merge(all, PP1_4, by="biotype")
all <- merge(all, PP2_1, by="biotype")
all <- merge(all, PP2_2, by="biotype")
all <- merge(all, PP2_3, by="biotype")
all <- merge(all, PP2_4, by="biotype")
all <- merge(all, PP3_1, by="biotype")
all <- merge(all, PP3_2, by="biotype")

all <- t(all)
colnames(all) <- all[1,]
all <- all[-1,]
data <- melt(all,id=rownames(all))
colnames(data) <- c("sample","biotype","count")
data$biotype <- factor(data$biotype, levels = c("protein_coding","lncRNA","miRNA","snRNA","snoRNA","scaRNA","misc_RNA","others"))
data$count <- as.numeric(data$count)

my20colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', 
               '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', 
               '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', 
               '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', 
               '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', 
               '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175',
               "#e5192c","#3a77b7","#3cac4c")

p1 <- ggplot(data, aes(x=biotype, y=count, fill=sample)) + 
        geom_bar(stat="identity",position = "stack") +
        scale_y_continuous(labels = scales::number_format(big.mark = ",")) +
        scale_fill_manual(values = my20colors) +
        theme(panel.background = element_rect(colour = "black", size = 1, fill = 'white', linetype = "solid"),
            text = element_text(family = "Arial"),
            axis.title = element_blank(),
            axis.text.x = element_text(family = 'Arial',size = 18,angle = 30, hjust = 1, vjust = 1,color = "black"),
            axis.text.y = element_text(family = 'Arial',size = 18,color = "black"),
            plot.title = element_text(family = 'Arial',size = 18, hjust = 0.5,color = "black"),
            plot.margin = margin(t=15,r=15,b=10,l=10),
            legend.position= c(0.93,0.73),
            legend.key.size = unit(20, "pt"),
            legend.text=element_text(size=20),
            legend.title=element_blank()) +
            theme( legend.key.size = unit(1, 'cm'), 
            legend.key.height = unit(0.5, 'cm'), 
            legend.key.width = unit(0.5, 'cm'), 
            legend.title = element_text(size=14), 
            legend.text = element_text(size=10)) +
            guides(col = guide_legend(nrow =6, byrow = TRUE)) + 
        annotate("text",x = 1,y=9000, label="16,888",size = 5, color="white") +
        annotate("text",x = 1,y=26000, label="18,457",size = 5, color="white") +
        annotate("text",x = 1,y=43000, label="16,730",size = 5, color="white") +
        annotate("text",x = 1,y=59500, label="17,509",size = 5, color="white") +
        annotate("text",x = 1,y=76500, label="17,137",size = 5, color="white") +
        annotate("text",x = 1,y=93500, label="17,445",size = 5, color="white") +
        annotate("text",x = 1,y=112000, label="18,204",size = 5, color="white") +
        annotate("text",x = 1,y=130500, label="17,367",size = 5, color="white") +
        annotate("text",x = 1,y=147500, label="17,534",size = 5, color="white") +
        annotate("text",x = 1,y=166000, label="16,090",size = 5, color="white") +
        annotate("text",x = 1,y=182500, label="17,811",size = 5, color="white") +
        annotate("text",x = 1,y=199500, label="17,578",size = 5, color="white") +
        annotate("text",x = 1,y=217000, label="18,005",size = 5, color="white") +
        annotate("text",x = 1,y=235000, label="18,366",size = 5, color="white") +
        annotate("text",x = 1,y=252500, label="17,421",size = 5, color="white") +
        annotate("text",x = 1,y=269500, label="17,897",size = 5, color="white") +
        annotate("text",x = 1,y=289000, label="18,178",size = 5, color="white") +
        annotate("text",x = 1,y=308000, label="18,054",size = 5, color="white") +
        annotate("text",x = 2,y=5000, label="10,647",size = 5, color="white") +
        annotate("text",x = 2,y=17800, label="15,671",size = 5, color="white") +
        annotate("text",x = 2,y=30600, label="9,940",size = 5, color="white") +
        annotate("text",x = 2,y=41900, label="12,606",size = 5, color="white") +
        annotate("text",x = 2,y=53700, label="11,405",size = 5, color="white") +
        annotate("text",x = 2,y=65200, label="12,303",size = 5, color="white") +
        annotate("text",x = 2,y=78700, label="14,762",size = 5, color="white") +
        annotate("text",x = 2,y=92200, label="12,053",size = 5, color="white") +
        annotate("text",x = 2,y=104200, label="12,462",size = 5, color="white") +
        annotate("text",x = 2,y=115200, label="9,506",size = 5, color="white") +
        annotate("text",x = 2,y=127700, label="13,243",size = 5, color="white") +
        annotate("text",x = 2,y=141500, label="12,602",size = 5, color="white") +
        annotate("text",x = 2,y=154800, label="14,065",size = 5, color="white") +
        annotate("text",x = 2,y=169300, label="15,186",size = 5, color="white") +
        annotate("text",x = 2,y=182800, label="12,151",size = 5, color="white") +
        annotate("text",x = 2,y=195500, label="13,588",size = 5, color="white") +
        annotate("text",x = 2,y=209500, label="14,652",size = 5, color="white") +
        annotate("text",x = 2,y=223800, label="14,269",size = 5, color="white")

        
data2 <- data[!(data$biotype %in% c("protein_coding","lncRNA","others")),]
data2 <- data2[data2$count != 0,]
p2 <- ggplot(data2, aes(x=biotype, y=count, fill=sample)) + 
        geom_bar(stat="identity",position = "stack") +
        coord_trans(y="sqrt") +
        scale_y_continuous(labels = scales::number_format(big.mark = ",")) +
        NoLegend() +
        scale_fill_manual(values = my20colors) +
        theme(panel.background = element_rect(colour = "black", size = 1, fill = 'white', linetype = "solid"),
            text = element_text(family = "Arial"),
            axis.title=element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1,family = 'Arial',size = 18,color = "black"),
            axis.text.y = element_text(family = 'Arial',size = 18,color = "black"),
            plot.title = element_text(family = 'Arial',size = 18, hjust = 0.5,color = "black"),
            plot.margin = margin(t=15,r=15,b=10,l=10))


p2_grob <- ggplotGrob(p2)
png("/public/home/chenhy/Tumor_heterogeneity/Result/All_gene_type.png", 1200, 1200, res=120)
p1 + annotation_custom(grob = p2_grob,xmin = 2.5,xmax = 7.5,ymin = 6000,ymax = 190000)
dev.off()

### gene and lncRNA marker(Refer to Figure 4A)
lncRNA <- read.table("./Human_gencode_v43_lncRNA.gtf", header = F, sep = '\t')
mRNA <- read.table("./Human_gencode_v43_mRNA.gtf", header = F, sep = '\t')

subset.matrix <- result_obj@raw.data[lncRNA$V4, ]
object2 <- CreateSeuratObject(subset.matrix)
orig.ident <- object1@meta.data
object2 <- AddMetaData(object = object2, metadata = orig.ident) 
object2 <- SetAllIdent(object = object2, id = ident.use)

all_marker <- FindAllMarkers(result_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

all_marker_mRNA <- all_marker[intersect(rownames(all_marker), unique(mRNA$V4)),]
all_marker_lncRNA <- all_marker[intersect(rownames(all_marker), unique(lncRNA$V4)),]
write.table(all_marker_lncRNA,"./PaperResult/lncRNA_marker.txt", sep = "\t")
write.table(all_marker_mRNA,"./PaperResult/mRNA_marker.txt", sep = "\t")

library(ClusterGVis)
library(RColorBrewer)
top_gene <- all_marker_mRNA %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC)

st.data <- prepareDataFromscRNA(object = result_obj,
                                diffData = top_gene,
                                showAverage = T)

markGenes = unique(top_gene$gene)[sample(1:length(unique(top_gene$gene)),40,replace = F)]
pdf("./mRNA_heatmap.pdf",width = 9,height = 8,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           ctAnno.col = brewer.pal(12, "Paired"),
           cluster.order = c(1:12),
           genes.gp = c('italic',12,NA),
           annnoblock.text = FALSE,
           add.sampleanno = FALSE,
           show_row_dend = F)
dev.off()

### Plot(Refer to Figure 5A)
data <- result_obj@meta.data %>%
   group_by(type,celltype) %>%
   summarise(n=n()) %>%
   mutate(freq = n/sum(n))
colnames(data) <- c("type","Cell_Type","Value","Proportion")
data$type <- factor(data$type, levels = c("Normal", "Tumor(Good Prognosis)", "Tumor(Poor Prognosis)"))

png("./PaperResult/celltype_proportion.png",700,500,res=120)
ggplot(data, aes(x = type,y=Proportion,fill = Cell_Type, 
                stratum = Cell_Type, alluvium = Cell_Type)) +
  geom_col(width = 0.5,color=NA)+
  geom_flow(width = 0.5,alpha = 0.4,knot.pos = 0)
  scale_fill_brewer(palette="Paired")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(expand = c(0.2,0)) +
  theme_bw()+
  theme(axis.text.x=element_text(size=9,family = "Arial",color = "black")) +
  theme(axis.text.x=element_text(angle=10,hjust=1),
    text = element_text(size=12,family = "Arial",color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(), 
    strip.text = element_text(size=12,family = "Arial",color = "black"),
    strip.background = element_blank(),
    axis.line = element_line(color = 'black'),
    axis.title.x = element_blank()
)
dev.off()


### Bhattacharyya distance plot(Refer to Figure 5B)
library(distdimscr)
result_obj_sub <- subset(result_obj, subset = type %in% c("Tumor(Poor Prognosis)","Tumor(Good Prognosis)"))
knitr::kable(table(result_obj$type,result_obj$celltype))

type <- c("Hepatocytes",
"T cells",
"Kupffer cells",
"Stellate cells",
"Endothelial cells",
"B/Plasma cells",
"Cholangiocytes")

pca_emb <- result_obj@reductions$iNMF@cell.embeddings

for(m in type){
      ln1 <- colnames(result_obj)[result_obj$celltype==m & result_obj$type=="Tumor(Poor Prognosis)"]
      ln0 <- colnames(result_obj)[result_obj$celltype==m & result_obj$type=="Tumor(Good Prognosis)"]

      ln1_pca <- pca_emb[ln1,]
      ln0_pca <- pca_emb[ln0,]

      bhatt.dist <- bhatt.dist.rand <- vector("logical",length=500)
      set.seed("1234")

      for (i in 1:500) {
            bhatt.dist[[i]] <- dim_dist(embed_mat_x=ln1_pca,embed_mat_y=ln0_pca,dims_use=1:20,num_cells_sample=500,distance_metric="bhatt_dist",random_sample=FALSE)
            bhatt.dist.rand[[i]] <- dim_dist(embed_mat_x=ln1_pca,embed_mat_y=ln0_pca,dims_use=1:20,num_cells_sample=500,distance_metric="bhatt_dist",random_sample=TRUE)
            }
      
      bhatt.dist <- data.frame(cells.distance=bhatt.dist,comparison="real")
      bhatt.dist.rand <- data.frame(cells.distance=bhatt.dist.rand,comparison="random")

      bhatt.res <- rbind(bhatt.dist,bhatt.dist.rand)
      
      path <- paste0("/public/home/chenhy/Tumor_heterogeneity/Result/PaperResult/",gsub(' ',"_",m),"_distance.png")
      png(path, 800, 700)
      p <- ggplot(bhatt.res,aes(x=comparison,y=cells.distance)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(size=0.5) +
      theme_bw() +
      xlab("Comparison type") +
      ylab("Bhattacharrya distance")
      print(p)
      dev.off()
      write.table(bhatt.res,paste0("/public/home/chenhy/Tumor_heterogeneity/Result/PaperResult/",gsub(' ',"_",m),"_distance.txt"),sep = "\t")
}

mean_diff <- dis %>%
  group_by(celltype, comparison) %>%
  summarise(mean = mean(cells.distance)) %>%
  spread(comparison, mean) %>%
  mutate(mean_diff = real / random)

library(gghalves)
ggplot(dis) +
  geom_half_violin(aes(celltype, cells.distance, split = comparison, fill = comparison),
                   position = "identity") +
  scale_fill_manual(values = c("#3182bd", "#e6550d")) +
  geom_half_boxplot(data = filter(dis, comparison == "real"),
                    aes(celltype, cells.distance),
                    width = 0.15,
                    side = "r") +
  geom_half_boxplot(data = filter(dis, comparison == "random"),
                    aes(celltype, cells.distance),
                    width = 0.15,
                    side = "l") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("/public/home/chenhy/Tumor_heterogeneity/Result/PaperResult/Distance.png",900,400,res=120)
ggplot(dis) +
  geom_half_violin(aes(celltype, cells.distance, split = comparison, fill = comparison),
                   position = "identity",width = 1) +
  scale_fill_manual(values = c("#3182bd", "#e6550d")) +
  theme(text = element_text(family = "Arial",color = "black")) +
  geom_half_boxplot(data = filter(dis, comparison == "real"),
                    aes(celltype, cells.distance),
                    width = 0.15,
                    side = "r") +
  geom_half_boxplot(data = filter(dis, comparison == "random"),
                    aes(celltype, cells.distance),
                    width = 0.15,
                    side = "l") +
  geom_text(data = mean_diff, aes(x = celltype, y = 0.4, label = sprintf("%.2f\nfold", mean_diff)),
            size = 4) +
  theme_minimal() +
  theme(text = element_text(family = "Arial",color = "black")) +
    theme(axis.line.x = element_line(color = "black", size = 0.5),  
        axis.line.y = element_line(color = "black", size = 0.5),  
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 10, color = "black"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_blank())
dev.off()

### DEG plot(Refer to Figure 5C)
result_obj <- RunDEtest(srt = result_obj, group_by = "celltype", fc.threshold = 1, only.pos = FALSE)
pdf("/public/home/chenhy/Tumor_heterogeneity/Result/PaperResult/DEG.pdf",12,10)
VolcanoPlot(srt = result_obj, group_by = "celltype", label.size = 2)
dev.off()

### Enrichment plot(Refer to Figure 5D)
result_obj <- RunEnrichment(
  srt = result_obj, group_by = "celltype", db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05"
)
EnrichmentPlot(
  srt = result_obj, group_by = "celltype", group_use = c("Hepatocytes", "T cells", "Kupffer cells","Stellate cells","Endothelial cells","B/Plasma cells","Cholangiocytes"),
  plot_type = "dot"
)