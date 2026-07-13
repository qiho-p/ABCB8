# ABCB8_Fig5 & FigS5

In this pipeline, we showed the detail codes in Fig5 & FigS5.

```
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(future.apply)
library(dplyr)
library(reticulate)
library(ReductionWrappers)
library(s2a)
library(cowplot)
library(ggplot2)
library(trqwe)
library(patchwork)
library(BuenColors)
library(paletteer)
library(ggsci)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(KEGG.db)
library(scCustomize)
library(Nebulosa)
library(SCP)
library(RColorBrewer)
library(ggpubr)
library(reshape2)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 200000 * 1024^2) 
```


#### Fig5A
```r

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)

UTX_CHIP <- read.csv(row.names=1,"./DEG/DESeq2_sizefactor_true/GSM2304478_UTX-IP.anno.csv")
length(unique(UTX_CHIP$SYMBOL))

me3_shA_vs_shR_up_promoter <- read.csv(row.names=1,file="./DEG/DESeq2_sizefactor_true/fc2.up.promoter.DI.peaks.shA_322.vs.shR.rmIgg.csv")

H3k27me3_up_gene <- as.character(unique(me3_shA_vs_shR_up_promoter$SYMBOL))
UTX_CHIP_gene <- as.character(unique(UTX_CHIP$SYMBOL))
ol_gene <- intersect(H3k27me3_up_gene,UTX_CHIP_gene)
nobinding <- setdiff(H3k27me3_up_gene,ol_gene)
length(H3k27me3_up_gene)
length(UTX_CHIP_gene)
length(ol_gene)
length(nobinding)

df <- data.frame(Var1=c("UTX_binding","UTX_non_binding"),
  Freq=c(length(ol_gene),length(nobinding)))
df$pct <- df$Freq / sum(df$Freq)
df$pct <- round(df$pct,4) * 100

p1 <- ggplot(df, aes(x = "", y = pct, fill = Var1)) + 
  geom_col(color = "grey") + 
  geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 7) + 
  scale_fill_brewer(palette = c("Set1")) + 
  coord_polar("y", start = 0)  + 
  theme_void()


shA_vs_shR <- read.csv(row.names=1,"./DEG/DESeq2_sizefactor_true/res_ckit_sh_ren_all.csv")
shA_vs_shR$group <- ifelse(abs(shA_vs_shR$log2FoldChange) > 0.1, ifelse(shA_vs_shR$log2FoldChange > 0.1 ,'Up','Down'),'NS')
table(shA_vs_shR$group)
shA_vs_shR_epi_gene <- intersect(rownames(shA_vs_shR),ol_gene)
length(shA_vs_shR_epi_gene)
shA_vs_shR_epi_gene_info <- shA_vs_shR[shA_vs_shR_epi_gene,]
dim(shA_vs_shR_epi_gene_info)
table(shA_vs_shR_epi_gene_info$group)

df <- data.frame(Var1=c("Down","NS","Up"),
  Freq=c(117,39,70))
df$pct <- df$Freq / sum(df$Freq)
# df$pct <- round(df$pct,3) * 100
df$pct <- round(df$pct,4) * 100

p2 <- ggplot(df, aes(x = "", y = pct, fill = Var1)) + 
  geom_col(color = "grey") + 
  geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 7) + 
  scale_fill_brewer(palette = c("Set1")) + 
  coord_polar("y", start = 0)  + 
  theme_void()


```

#### Fig5B
```r


GO_ol_gene <- read.csv(row.names=1,"./DESeq2_sizefactor_true/GO.Utx.binding.K27.up.csv")
GO_nobinding <- read.csv(row.names=1,"./DESeq2_sizefactor_true/GO.Utx.nonbinding.K27.up.csv")

GO_ol_gene <- as.data.frame(GO_ol_gene)
GO_nobinding <- as.data.frame(GO_nobinding)

GO_ol_gene$Description <- as.character(GO_ol_gene$Description)
GO_nobinding$Description <- as.character(GO_nobinding$Description)
ptw <- c("regulation of leukocyte differentiation","positive regulation of leukocyte differentiation","regulation of myeloid cell differentiation","myeloid leukocyte differentiation","regulation of myeloid leukocyte differentiation")

GO_ol_gene <- GO_ol_gene[which(GO_ol_gene$Description %in% ptw),]
GO_nobinding <- GO_nobinding[which(GO_nobinding$Description %in% ptw),]
GO_ol_gene$Cluster <- "Utx_binding"
GO_nobinding$Cluster <- "Utx_non_binding"
myeloid_info <- rbind(GO_ol_gene,GO_nobinding)
myeloid_info$Description <- factor(myeloid_info$Description,levels=rev(c("regulation of leukocyte differentiation","myeloid leukocyte differentiation","positive regulation of leukocyte differentiation","regulation of myeloid cell differentiation","regulation of myeloid leukocyte differentiation")))

p <- ggplot(myeloid_info, aes(Cluster, Description)) +
  geom_point(aes(fill = pmin(pvalue, 0.2), size = Count), shape = 21) +
  scale_fill_gradientn(
    colours = c("#E41A1C", "#FB8072", "#FDB462", "#FEE6CE", "#F0F0F0"),  # 添加更浅的颜色
    values = scales::rescale(c(0, 0.01, 0.05, 0.1, 0.2)),  # 更细分的阈值
    breaks = c(0.01, 0.05, 0.1, 0.15, 0.2),
    labels = c("0.01", "0.05", "0.10", "0.15", "≥0.20"),
    guide = guide_colorbar(reverse = TRUE, title = "pvalue")
  ) +
  scale_size_continuous(range = c(2, 10), name = "Count") +   # 调整点的大小范围
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text = element_text(color = 'black', size = 10)
  )


```



#### Fig5D
```r


AML_Utx_KO_v_WT_result <- read.csv(file="./Utx_05_KO_v_WT_count_DEseqnormalize_symbol_and_anno.csv",row.names=1)
AML_Utx_KO_v_WT_result_p <- AML_Utx_KO_v_WT_result[which(AML_Utx_KO_v_WT_result$KO_v_WT_pvalue < 0.05),]
AML_Utx_KO_v_WT_DN_all <- AML_Utx_KO_v_WT_result_p[which(AML_Utx_KO_v_WT_result_p$KO_v_WT_log2FoldChange < 0),]
AML_Utx_KO_v_WT_DN_gene <- rownames(AML_Utx_KO_v_WT_DN_all)

true.3.group.norm.counts <- read.csv(row.names=1,"./FAS_RNASeq_20240813/DEG/norm.counts.FAS.3.group.csv")
true.3.group.norm.counts$shA_vs_shR_FC <- rowMeans(true.3.group.norm.counts[,4:6]) / rowMeans(true.3.group.norm.counts[,1:3])
true.3.group.norm.counts$FAS_vs_Veh_FC <- rowMeans(true.3.group.norm.counts[,7:9]) / rowMeans(true.3.group.norm.counts[,4:6])
true.3.group.norm.counts$FC_score <- true.3.group.norm.counts$FAS_vs_Veh_FC - true.3.group.norm.counts$shA_vs_shR_FC

tmp <- true.3.group.norm.counts[intersect(AML_Utx_KO_v_WT_DN_gene,rownames(true.3.group.norm.counts)),]
tmp <- tmp[order(tmp$FC_score,decreasing=TRUE),]

library(circlize)
# 自定义颜色函数
col_fun <- colorRamp2(
  breaks = c(-1.5, -0.4, -0.1, 0.4, 1.5),   # break 点
  colors = c("#2166ac", "#f7fbff", "#fff5f0", "#d6604d", "#b2182b")  # 对应颜色
)

library(ComplexHeatmap)
library(BuenColors)
zscore <- log2(tmp[,1:9]+1)
zscore <- na.omit(zscore)
zscore <- t(apply(zscore, 1, function(x) (x-mean(x))/sd(x)))
zscore <- na.omit(zscore)
zscore[zscore < -1.5] <- -1.5
zscore[zscore > 1.5] <- 1.5
gene <- c("Ets1","Batf","Meis2","Spib","Lef1","Il17rb","Zfp36","Csf1r","Spi1","Gfi1b")
genemark <- which(rownames(zscore) %in% gene)
labs <- rownames(zscore)[genemark]
ha1 <-  rowAnnotation(
  foo = anno_mark(at = genemark,
  labels = labs,
  labels_gp = gpar(fontsize = 10)
  ))
Heatmap(zscore, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    right_annotation = ha1,
    left_annotation = NULL,
    show_row_names = FALSE,
    show_column_names = TRUE,
    col = col_fun,
    row_names_gp = gpar(fontsize = 12),
    top_annotation = NULL,
    bottom_annotation = NULL ,
    column_names_rot  = 90,
    row_names_side = "right",
    column_title = "AML_Utx_KO_v_WT_DN_fc0_all_gene"
    )


```


#### Fig5O
```r

GO_up <- read.csv(row.names=1,"./EZH2i_RNASeq_20250429/DEG/GO_mmu_AML2_EZH2i_vs_DMSO_up.csv")

GO_shA_show <- GO_up[which(GO_up$Description %in% c("stem cell differentiation","myeloid leukocyte differentiation",
  "regulation of hemopoiesis","regulation of leukocyte differentiation","myeloid cell differentiation","leukocyte cell-cell adhesion","regulation of leukocyte cell-cell adhesion")),]
GO_shA_show$log10p.adjust <- -log10(GO_shA_show$p.adjust)
unique(GO_shA_show$Description)
p1 <- ggbarplot(GO_shA_show,x="Description",y="log10p.adjust",orientation = "horiz",color = "red",
  fill="red",sort.val = "asc",title="GO_mmu_AML2_EZH2i_vs_DMSO_up",alpha=0.7)


```

<img src="./Fig5_images/GO_EZH2i_AML2.png" alt="GO_EZH2i_AML2.png" width="400" />


#### Fig5P
```r

geneList <- fread("./EZH2i_RNASeq_20250429/DEG/GSEA/mmu.AML2.EZH2i_vs_DMSO.results.rnk")
geneList <- as.data.frame(geneList)
geneList_ <- geneList$V2
names(geneList_) <- geneList$V1

shA_gmt <- read.gmt("./EZH2i_RNASeq_20250429/DEG/GSEA/shA_gene_sets.gmt")
FAS_gmt <- read.gmt("./EZH2i_RNASeq_20250429/DEG/GSEA/FAS_gene_sets.gmt")
UTX_gmt <- read.gmt("./EZH2i_RNASeq_20250429/DEG/GSEA/UTX_gene_sets.gmt")

all_gmt <- rbind(shA_gmt,FAS_gmt,UTX_gmt)

gsea_ABCB8 <- GSEA(geneList_,  #待富集的基因列表
    TERM2GENE = all_gmt,  #基因集
    pvalueCutoff = 0.05,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000,
    minGSSize = 15,
     by = "fgsea")  #指定 p 值校正方法

library(enrichplot)

p5 <- gseaplot2(gsea_ABCB8, 
                geneSetID = c("top_200AML_Utx_KO_v_WT_down_no_p_filt.grp","top_200shA.FAS_vs_shA.Veh_up.grp","top_200shA_vs_shR_down.grp"), #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                color = c("#39a26d","#36638f"),
                pvalue_table = FALSE,
                ES_geom = "line",
                rel_heights = c(1.5, 0.5, 1), #子图高度
                subplots = 1:2
)

```

<img src="./Fig5_images/gseaplot.GSEA.EZH2i_AML2.png" alt="gseaplot.GSEA.EZH2i_AML2.png" width="350" />



#### FigS5A
```r

GOupres_1_all <- mcreadRDS("./GO.utx_dn.rds")
GO_ckit_sha_dn <- mcreadRDS("./GO.shA_dn.rds")

utx_go_dn <- GOupres_1_all@result
utx_go_dn <- utx_go_dn[which(utx_go_dn$p.adjust < 0.05),]
sha_go_dn <- GO_ckit_sha_dn@result
sha_go_dn <- sha_go_dn[which(utx_go_dn$p.adjust < 0.05),]

ol_pathway <- intersect(sha_go_dn$ID,utx_go_dn$ID)

utx_ol_go_dn <- utx_go_dn[which(utx_go_dn$ID %in% ol_pathway),]
utx_ol_go_dn$log10p.adjust <- -log10(utx_ol_go_dn$p.adjust)

p <- ggbarplot(utx_ol_go_dn[1:20,],x="Description",y="log10p.adjust",orientation = "horiz",color = "red",
  fill="red",sort.val = "asc")


```

<img src="./Fig5_images/shUtx_shA_common_dn.jpg" alt="shUtx_shA_common_dn.jpg" width="350" />


#### FigS5F
```r


mmu.EZH2i_vs_DMSO.results <- read.csv(row.names=1,file="./EZH2i_RNASeq_20250429/DEG/mmu.AML2.EZH2i_vs_DMSO.results.csv")
mmu.EZH2i_vs_DMSO.results$group <- ifelse(mmu.EZH2i_vs_DMSO.results$pvalue < 0.05 & abs(mmu.EZH2i_vs_DMSO.results$log2FoldChange) > 1, 
  ifelse(mmu.EZH2i_vs_DMSO.results$log2FoldChange > 1, "Up", "Down"), "NS")
table(mmu.EZH2i_vs_DMSO.results$group)

show_data <- mmu.EZH2i_vs_DMSO.results[!is.na(mmu.EZH2i_vs_DMSO.results$group),]
show_data$log2baseMean <- log2(show_data$baseMean)
show_data$log10pvalue <- -log10(show_data$pvalue)
show_data[which(show_data$log10pvalue > 250),]$log10pvalue <- 250
show_data[which(show_data$log2FoldChange > 12),]$log2FoldChange <- 12

final_marker <- c("Il6ra","Gata2","Gli3","Tal1","Cebpa","Cebpb")

show_data$label <- ""
show_data[which(rownames(show_data) %in% final_marker & show_data$group != "NS"),]$label <- rownames(show_data[which(rownames(show_data) %in% final_marker & show_data$group != "NS"),])

p <- ggplot(show_data, aes(log2FoldChange,log10pvalue,colour=group))+ labs(x="log2FoldChange",y="log10pvalue")
p <- p + geom_point(alpha=1, size=2)  +
   theme_pubr()+ scale_color_manual(values=c("#2971B1","grey","#BB2933")) + 
  geom_text_repel(data = show_data, aes(x = log2FoldChange,
  y = log10pvalue,
  label = label),
  size = 4,box.padding = unit(3, "lines"),color="black",
  point.padding = unit(0.5, "lines"),
  show.legend = FALSE, max.overlaps = Inf) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) + ggtitle("mmu.AML2: EZH2i vs DMSO") + xlim(-12,12)


```

<img src="./Fig5_images/maplot_AML2.png" alt="maplot_AML2.png" width="350" />



















