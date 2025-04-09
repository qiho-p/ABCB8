
suppressPackageStartupMessages({
  require(Rsamtools)
  require(GenomicFeatures)
  require(GenomicAlignments)
  require(BiocParallel)
  require(pheatmap)
  require(RColorBrewer)
  require(PoiClaClu)
  require(org.Mm.eg.db)
  require(AnnotationDbi)
  require(DOSE)
  require(clusterProfiler)
  require(topGO)
  require(pathview)
  require(org.Hs.eg.db)
  require(AnnotationDbi)
  require(DOSE)
  require(clusterProfiler)
  require(topGO)
  require(ggplot2)
  require(DESeq2)
})
library(BuenColors)
library(future)
library(future.apply)
library(stringr)
library(reshape2)
library(ggpubr)
library(trqwe)
library(Seurat)

******************************************************************************************************
********************************** Fig3 GSEA ***************************************************
********************************** Fig3 GSEA ***************************************************
********************************** Fig3 GSEA ***************************************************
******************************************************************************************************

gsea_df <- read.csv(row.names=1,"./Abcb8_zjn/FAS_RNASeq_20240813/DEG/GSEA.shAbcb8.FAS_vs_Veh.10000.pdf")

filtered_gsea_df <- gsea_df[grep("DIFFERENTIATION",gsea_df$ID), ]
filtered_gsea_results <- gsea_ABCB8
filtered_gsea_results@result <- filtered_gsea_df

p5 <- gseaplot2(filtered_gsea_results, 
                geneSetID = c(4,6), #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                color = c("#F4A460", "#CD5C5C"),
                pvalue_table = TRUE,
                ES_geom = "line",
                rel_heights = c(1.5, 0.5, 1), #子图高度
                subplots = 1:2
)

pdf(width=6,height=4,"./Abcb8_zjn/FAS_RNASeq_20240813/DEG/gseaplot.GSEA.shAbcb8.FAS_vs_Veh.pdf")
p5
dev.off()


******************************************************************************************************
********************************** Fig3 ***************************************************
********************************** Fig3 ***************************************************
********************************** Fig3 ***************************************************
******************************************************************************************************


shA_vs_shR_all <- read.csv(row.names=1,file='./Abcb8_zjn/FAS_RNASeq_20240813/DEG/shAbcb8_vs_shRen.all.csv')
shA_vs_shR_all$group <- as.factor(ifelse(abs(shA_vs_shR_all$log2FoldChange) > 0 & shA_vs_shR_all$pvalue < 0.05, ifelse(shA_vs_shR_all$log2FoldChange > 0 ,'Up','Down'),'NS'))
table(shA_vs_shR_all$group)
 # Down    NS    Up
 # 2381 12672  2565
shA_vs_shR_down <- shA_vs_shR_all[which(shA_vs_shR_all$group == "Down"),]

shA.FAS_vs_shA.Veh_all <- read.csv(row.names=1,file="./Abcb8_zjn/FAS_RNASeq_20240813/DEG/shAbcb8.FAS_vs_Veh.all.csv")
shA.FAS_vs_shA.Veh_all$group <- as.factor(ifelse(abs(shA.FAS_vs_shA.Veh_all$log2FoldChange) > 0 & shA.FAS_vs_shA.Veh_all$pvalue < 0.05, ifelse(shA.FAS_vs_shA.Veh_all$log2FoldChange > 0 ,'Up','Down'),'NS'))
table(shA.FAS_vs_shA.Veh_all$group)
 # Down    NS    Up
 #  506 16049  1063
shA.FAS_vs_shA.Veh_up <- shA.FAS_vs_shA.Veh_all[which(shA.FAS_vs_shA.Veh_all$group == "Up"),]

library(Vennerable)
data<-Venn(list(shA_vs_shR_down=rownames(shA_vs_shR_down),shA.FAS_vs_shA.Veh_up=rownames(shA.FAS_vs_shA.Veh_up)))
plot(data,doWeight=T,show = list(Faces = FALSE, DarkMatter = FALSE))

pdf("./Abcb8_zjn/FAS_RNASeq_20240813/DEG/Venn.shA_vs_shA.shA_FAS_vs_Veh.p0.05.fc0.pdf")
plot(data,doWeight=T,show = list(Faces = FALSE, DarkMatter = FALSE))
dev.off()
enrich_pvalue(24421,2073,585,478)  #-overlap
enrich_pvalue(24421,1364,261,236)  #-overlap
enrich_pvalue(24421,2809,969,1596)  #-overlap


plot.gene <- intersect(rownames(shA_vs_shR_down),rownames(shA.FAS_vs_shA.Veh_up))

sample.counts <- read.csv(row.names=1,file='./Abcb8_zjn/FAS_RNASeq_20240813/DEG/norm.counts.3.group.csv')

final_marker <- c("Cebpb","Il18bp","Il18r1","Il2rg","Lcn2","Mmp25","Tlr2","Batf","Batf2","Lef1","Spib","Zfp36","Trib1","Lif","Inhba")

library(ComplexHeatmap)
library(BuenColors)
zscore <- log2(sample.counts[plot.gene,]+1)
zscore <- na.omit(zscore)
zscore <- t(apply(zscore, 1, function(x) (x-mean(x))/sd(x)))
zscore <- na.omit(zscore)
zscore[zscore < -2] <- -2
zscore[zscore > 2] <- 2
aa <- jdb_palette("brewer_yes")
color <- colorRampPalette(aa)(30)
library(circlize)
heatmap_colors <- colorRamp2(c(-2,-1,-0.5,0,0.1,0.25,0.5,1,2), c("#053061", "#2971B1", "#6AACD0",
  "#F4BAA2","#EFA68C","#E27B63","#CA4A46","#BF3238","#AC212F"))
show_col(aa)
show_col(color)

gene <- final_marker
genemark <- which(rownames(zscore) %in% gene)
labs <- rownames(zscore)[genemark]
ha1 <-  rowAnnotation(
  foo = anno_mark(at = genemark,
  labels = labs,
  labels_gp = gpar(fontsize = 12)
  ))

ha_t <- HeatmapAnnotation(
  sample = c(rep("shRen_Veh",3),rep("shAbcb8_Veh",3),rep("shAbcb8_FAS",3)),which="column",
  col = list( 
    sample = c("shRen_Veh" = "red", 
    "shAbcb8_Veh" = "blue",
    "shAbcb8_FAS" = "yellow") #分类附颜色
  ))

p <- Heatmap(zscore,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    right_annotation = ha1,
    left_annotation = NULL,
    show_row_names = FALSE,
    show_column_names = FALSE,
    col = heatmap_colors,
    row_names_gp = gpar(fontsize = 12),
    top_annotation = ha_t,
    bottom_annotation = NULL ,
    column_names_rot  = 45,
    row_names_side = "right",
    column_title = NULL
    )

pdf(width=5,height=5,"./Abcb8_zjn/FAS_RNASeq_20240813/DEG/FAS.heatmap.3.group.diff.genes.p0.05fc0.v2.pdf")
p
dev.off()


******************************************************************************************************
********************************** Fig3D ***************************************************
********************************** Fig3D ***************************************************
********************************** Fig3D ***************************************************
******************************************************************************************************


shA_vs_shR_all <- read.csv(row.names=1,file='./Abcb8_zjn/FAS_RNASeq_20240813/DEG/shAbcb8_vs_shRen.all.csv')
shA_vs_shR_all$group <- as.factor(ifelse(abs(shA_vs_shR_all$log2FoldChange) > 0 & shA_vs_shR_all$pvalue < 0.05, ifelse(shA_vs_shR_all$log2FoldChange > 0 ,'Up','Down'),'NS'))
table(shA_vs_shR_all$group)
shA_vs_shR_down <- shA_vs_shR_all[which(shA_vs_shR_all$group == "Down"),]
shA_vs_shR_up <- shA_vs_shR_all[which(shA_vs_shR_all$group == "Up"),]

shA.FAS_vs_shA.Veh_all <- read.csv(row.names=1,file="./Abcb8_zjn/FAS_RNASeq_20240813/DEG/shAbcb8.FAS_vs_Veh.all.csv")
shA.FAS_vs_shA.Veh_all$group <- as.factor(ifelse(abs(shA.FAS_vs_shA.Veh_all$log2FoldChange) > 0 & shA.FAS_vs_shA.Veh_all$pvalue < 0.05, ifelse(shA.FAS_vs_shA.Veh_all$log2FoldChange > 0 ,'Up','Down'),'NS'))
table(shA.FAS_vs_shA.Veh_all$group)
shA.FAS_vs_shA.Veh_up <- shA.FAS_vs_shA.Veh_all[which(shA.FAS_vs_shA.Veh_all$group == "Up"),]
shA.FAS_vs_shA.Veh_down <- shA.FAS_vs_shA.Veh_all[which(shA.FAS_vs_shA.Veh_all$group == "Down"),]

library(Vennerable)
data<-Venn(list(shA_vs_shR_down=rownames(shA_vs_shR_down),shA.FAS_vs_shA.Veh_up=rownames(shA.FAS_vs_shA.Veh_up)))
plot(data,doWeight=T,show = list(Faces = FALSE, DarkMatter = FALSE))


shA.FAS_vs_shA.Veh_up_gene <- as.character(rownames(shA.FAS_vs_shA.Veh_up))
infoo <- shA_vs_shR_all[intersect(shA.FAS_vs_shA.Veh_up_gene,rownames(shA_vs_shR_all)),]

df <- data.frame(table(infoo$group))
df$pct <- df$Freq / sum(df$Freq)
df$pct <- round(df$pct,4) * 100

p2 <- ggplot(df, aes(x = "", y = pct, fill = Var1)) + 
  geom_col(color = "grey") + 
  geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 7) + 
  scale_fill_brewer(palette = c("Set1")) + 
  coord_polar("y", start = 0)  + 
  theme_void()

pdf("./Abcb8_zjn/FAS_RNASeq_20240813/DEG/FAS_shA_pie.pdf")
p2
dev.off()
