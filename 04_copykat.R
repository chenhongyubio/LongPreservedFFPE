# copykat
library(copykat)
library(Seurat)
library(tidyverse)
library(phylogram)
rm(list = ls())
source("/workcenters/workcenter2/chenhy/tumor_heterogeneity/copykat2.r")
setwd("/workcenters/workcenter2/chenhy/tumor_heterogeneity/")
obj <-readRDS("obj_all_liger_anno.rds")
counts <- as.matrix(obj@assays$RNA@counts)
normal_name <- colnames(subset(obj, subset = sample %in% c("GP2_4","PP1_4","PP2_4")))
setwd("/workcenters/workcenter2/chenhy/tumor_heterogeneity/Result/02_copykat/all_CNV/")
copykat <- copykat2(rawmat=counts, ngene.chr=5, sam.name="all_CNV", n.cores=8, norm.cell.names=normal_name)
saveRDS(copykat, "all_CNV_copykat.rds")
mallignant <- read.delim("all_CNV_copykat_prediction.txt", row.names = 1)
obj <- AddMetaData(obj, metadata = mallignant)
png("./p1.png",500,500)
DimPlot(obj, group.by = "copykat.pred")
dev.off()
ggsave("all_CNV_pred_mallignant.pdf", width = 12, height = 5)

pred.test <- data.frame(copykat$prediction)
CNA.test <- data.frame(copykat$CNAmat)

# CNV热图
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)
rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]
cells <- rbind(pred,pred)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
png("all_CNV_heatmap2.png", 2000,1200)
heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), 
          hclustfun = function(x) hclust(x, method="ward.D2"),
          ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",
          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")
dev.off()