# ABCB8_Fig4 & S4

In this pipeline, we showed the detail codes in Fig4 & S4.

```
#load packages
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
options(future.globals.maxSize = 200000 * 1024^2) # for 200 Gb RAM
```


#### Fig4E
```

computeMatrix reference-point --referencePoint center -b 5000 -a 5000 \
-R./sorted_ckit_H3K27me3/DEG/peaks.shA_322.vs.shR.rmIgg.bed \
-S H3K27me3_IgG_sorted_ckit_2.Nonenorm.50.bw H3K27me3_shRen_2.merge.scaleFactor.Nonenorm.50.real.bw H3K27me3_shAbcb8_322.merge.scaleFactor.Nonenorm.50.real.bw \
--numberOfProcessors 30 --skipZeros \
-o ./plots/center.peaks.shA_322.vs.shR.rmIgg.real.mat.gz 

plotHeatmap -m ./plots/center.peaks.shA_322.vs.shR.rmIgg.real.mat.gz \
 -out ./plots/center.peaks.shA_322.vs.shR.rmIgg.real.pdf \
 --colorList 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red'

```


#### Fig4F
```

computeMatrix reference-point --referencePoint center -b 5000 -a 5000 \
-R./sorted_ckit_H3K27me3/DEG/DESeq2_sizefactor_true/peaks.shA_322.vs.shR.rmIgg.up.p0.05fc1.5.bed \
-S H3K27me3_shRen_2.merge.scaleFactor.Nonenorm.50.real.bw H3K27me3_shAbcb8_322.merge.scaleFactor.Nonenorm.50.real.bw \
--numberOfProcessors 30 --skipZeros \
-o ./plots/center.correct.shA_322.vs.shR.rmIgg.up.p0.05fc1.5.real.mat.gz 

plotHeatmap -m ./plots/center.correct.shA_322.vs.shR.rmIgg.up.p0.05fc1.5.real.mat.gz \
 -out ./plots/center.correct.shA_322.vs.shR.rmIgg.up.p0.05fc1.5.real.pdf \
 --colorList 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red'

```


#### FigS4B
```

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)

shA_vs_shR_res_anno <- read.csv(row.names=1,file="./sorted_ckit_H3K27me3/DEG/DESeq2_sizefactor_true/res.peaks.shA_322.vs.shR.rmIgg.csv")

me3_shA_vs_shR_up <- shA_vs_shR_res_anno[which(shA_vs_shR_res_anno$pvalue < 0.05 & shA_vs_shR_res_anno$log2FoldChange > 1.5),]

me3_shA_vs_shR_up_anno <- as(me3_shA_vs_shR_up,"GRanges")
me3_shA_vs_shR_up_anno <- annotatePeak(me3_shA_vs_shR_up_anno, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db",verbose=FALSE)
plotAnnoPie(me3_shA_vs_shR_up_anno)


```


#### Fig4G
```


computeMatrix reference-point --referencePoint TSS -b 5000 -a 5000 \
-R ./sorted_ckit_H3K27me3/DEG/Fig4I.match.shA_vs_shR_down.sig.fc1.5.bed \
-S H3K27me3_shRen_2.merge.scaleFactor.Nonenorm.50.real.bw H3K27me3_shAbcb8_322.merge.scaleFactor.Nonenorm.50.real.bw \
--numberOfProcessors 30 --skipZeros \
-o ./sorted_ckit_H3K27me3/bw/plots/TSS.me3.shA_322.vs.shR.Fig4I.match.shA_vs_shR_down.sig.fc1.5.real.mat.gz 

plotHeatmap -m ./sorted_ckit_H3K27me3/bw/plots/TSS.me3.shA_322.vs.shR.Fig4I.match.shA_vs_shR_down.sig.fc1.5.real.mat.gz \
 -out ./sorted_ckit_H3K27me3/bw/plots/TSS.me3.shA_322.vs.shR.Fig4I.match.shA_vs_shR_down.sig.fc1.5.real.pdf \
 --colorList 'white,purple' 'white,purple' 'white,purple' 'white,purple' 'white,purple' 'white,purple' 'white,purple' 


```



#### Fig4H
```


library(ComplexHeatmap)
library(circlize)

res_ckit_all <- read.csv(file="./shABCB8_vs_shRen_all.csv",row.names=1)
res_ckit_all$group <- ifelse(res_ckit_all$pvalue < 0.05 & abs(res_ckit_all$log2FoldChange) > 1.5, 
  ifelse(res_ckit_all$pvalue < 0.05 & res_ckit_all$log2FoldChange > 1.5, "UP", "DOWN"), "NS")
RNA.down <- res_ckit_all[which(res_ckit_all$group == "DOWN"),]
nrow(RNA.down)

shA_vs_shR_res_anno <- read.csv(row.names=1,file="./DESeq2_sizefactor_true/res.peaks.shA_322.vs.shR.rmIgg.csv")
shA_vs_shR_res_promoter <- shA_vs_shR_res_anno[grepl("Promoter",shA_vs_shR_res_anno$annotation),]
shA_vs_shR_res_promoter_gene <- shA_vs_shR_res_promoter %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = abs(baseMean), n = 1)
shA_vs_shR_res_promoter_gene <- as.data.frame(shA_vs_shR_res_promoter_gene)
nrow(shA_vs_shR_res_promoter_gene)
length(unique(shA_vs_shR_res_promoter_gene$SYMBOL))

df_merge <- merge(
  RNA.down,
  shA_vs_shR_res_promoter_gene,
  by.x = "row.names",
  by.y = "SYMBOL"
)

df_merge <- df_merge[,c("Row.names","log2FoldChange.x","log2FoldChange.y","pvalue.x","pvalue.y")]
colnames(df_merge) <- c("Gene","log2FoldChange.RNA","log2FoldChange.H3K27me3","pvalue.RNA","pvalue.H3K27me3")

mat <- df_merge
mat <- mat[,c("log2FoldChange.H3K27me3","log2FoldChange.RNA")]
mat <- mat[order(mat[, "log2FoldChange.H3K27me3"], decreasing = TRUE), ]

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166ac", "white", "#b2182b")
)

Heatmap(
  mat,
  name = "log2FC",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = "RNA: shA vs shR DOWN genes",
  row_names_gp = gpar(fontsize = 6),
  show_row_names = FALSE
)


```


#### Fig4I
```

GO_up <- read.csv("./DEG/DESeq2_sizefactor_true/GO_H3k27me3_promoter_shA.322_vs_shR_up.fc1.5p0.05.csv",row.names=1)
GO_up$log10p.adjust <- -log10(GO_up$p.adjust)
stem_diff <- c("stem cell differentiation","myeloid leukocyte migration","regulation of leukocyte differentiation","positive regulation of leukocyte differentiation","myeloid leukocyte differentiation","regulation of myeloid cell differentiation","positive regulation of myeloid leukocyte differentiation")

p <- ggbarplot(GO_up[which(GO_up$Description %in% stem_diff),],x="Description",y="log10p.adjust",orientation = "horiz",color = "red",
  fill="red",sort.val = "asc",title="differentiation related pathways")


```


#### Fig4J
```

shA_vs_shR_res_anno <- read.csv(row.names=1,file="./sorted_ckit_H3K27me3/DEG/DESeq2_sizefactor_true/res.peaks.shA_322.vs.shR.rmIgg.csv")
shA_vs_shR_res_promoter <- shA_vs_shR_res_anno[grepl("Promoter",shA_vs_shR_res_anno$annotation),]
shA_vs_shR_res_promoter_gene <- shA_vs_shR_res_promoter %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = abs(baseMean), n = 1)
shA_vs_shR_res_promoter_gene <- as.data.frame(shA_vs_shR_res_promoter_gene)
nrow(shA_vs_shR_res_promoter_gene)
length(unique(shA_vs_shR_res_promoter_gene$SYMBOL))

shA_vs_shR_res_promoter_gene$GSEA.score <- (shA_vs_shR_res_promoter_gene$log2FoldChange)^3 * shA_vs_shR_res_promoter_gene$baseMean
rnk_shA_vs_shR_res_promoter_gene <- shA_vs_shR_res_promoter_gene[order(shA_vs_shR_res_promoter_gene$GSEA.score,decreasing=TRUE),]

geneList_ <- rnk_shA_vs_shR_res_promoter_gene$GSEA.score
names(geneList_) <- rnk_shA_vs_shR_res_promoter_gene$SYMBOL
geneList_ <- geneList_[!is.na(names(geneList_))]

all_GSEA_GMT <- read.gmt("./msigdb.v2022.1.Mm.symbols.gmt")

gsea_ABCB8 <- GSEA(geneList_,  #待富集的基因列表
    TERM2GENE = all_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

p5 <- gseaplot2(gsea_ABCB8, 
                geneSetID = c("GOBP_STEM_CELL_DIFFERENTIATION","GOBP_REGULATION_OF_MYELOID_CELL_DIFFERENTIATION"), #或直接输入基因集ID向量名，如c("hsa04610","hsa00260")
                pvalue_table = TRUE,
                ES_geom = "line",
                rel_heights = c(1.5, 0.5, 1), #子图高度
                subplots = 1:2
)


```


#### Fig4N-O
```


computeMatrix reference-point --referencePoint center -b 5000 -a 5000 \
-R ./FAS_H3K27me3/DEG/peaks.6.samples.bed \
-S ./bw/H3K27me3_IgG_sorted_ckit_2.Nonenorm.50.bw shRen_DMSO.2sample.merge.scaleFactor.Nonenorm.50.bw shRen_DFP.2sample.merge.scaleFactor.Nonenorm.50.bw shRen_DFP_FAS.2sample.merge.scaleFactor.Nonenorm.50.bw \
--numberOfProcessors 30 --skipZeros \
-o ./FAS_H3K27me3/bw/plots/center.peaks.6.samples.mat.gz 

plotHeatmap -m ./FAS_H3K27me3/bw/plots/center.peaks.6.samples.mat.gz \
 -out ./FAS_H3K27me3/bw/plots/center.peaks.6.samples.pdf \
 --colorList 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' 'white,red' \
 --samplesLabel IgG shRen_DMSO shRen_DFP shRen_DFP_FAS \
 --zMax 25


```


#### Fig4P
```

all.samples.norm.counts <- mcreadRDS("./FAS_H3K27me3/DEG/only.6.samples.norm.counts.rds")

DFP_vs_DMSO_res_anno <- read.csv(row.names=1,file="./FAS_H3K27me3/DEG/res.all.peaks.DFP_vs_DMSO_2samples.csv")
DFP_FAS_vs_DFP_res_anno <- read.csv(row.names=1,file="./FAS_H3K27me3/DEG/res.all.peaks.DFP_FAS_vs_DFP_2samples.csv")
DFP_vs_DMSO_res_anno$ID <- paste0("ID_",DFP_vs_DMSO_res_anno$seqnames,"_",DFP_vs_DMSO_res_anno$start,"_",DFP_vs_DMSO_res_anno$end)
DFP_FAS_vs_DFP_res_anno$ID <- paste0("ID_",DFP_FAS_vs_DFP_res_anno$seqnames,"_",DFP_FAS_vs_DFP_res_anno$start,"_",DFP_FAS_vs_DFP_res_anno$end)

DFP_vs_DMSO_res_anno$group <- ifelse(abs(DFP_vs_DMSO_res_anno$log2FoldChange) > 1 & DFP_vs_DMSO_res_anno$pvalue < 0.05, ifelse(DFP_vs_DMSO_res_anno$log2FoldChange > 1 ,'Up','Down'),'NS')
DFP_FAS_vs_DFP_res_anno$group <- ifelse(abs(DFP_FAS_vs_DFP_res_anno$log2FoldChange) > 1 & DFP_FAS_vs_DFP_res_anno$pvalue < 0.05, ifelse(DFP_FAS_vs_DFP_res_anno$log2FoldChange > 1 ,'Up','Down'),'NS')

DFP_vs_DMSO_DEG <- DFP_vs_DMSO_res_anno[which(DFP_vs_DMSO_res_anno$group != "NS"),]
DFP_FAS_vs_DFP_DEG <- DFP_FAS_vs_DFP_res_anno[which(DFP_FAS_vs_DFP_res_anno$group != "NS"),]
table(DFP_vs_DMSO_DEG$group)
table(DFP_FAS_vs_DFP_DEG$group)

candidate_IDs <- union(DFP_vs_DMSO_DEG$ID,DFP_FAS_vs_DFP_DEG$ID)
length(candidate_IDs)

candidate_IDs_exp <- all.samples.norm.counts[candidate_IDs,]
candidate_IDs_exp <- na.omit(candidate_IDs_exp)
dim(candidate_IDs_exp)

aa <- merge(DFP_vs_DMSO_res_anno[,c("pvalue","log2FoldChange","ID")],DFP_FAS_vs_DFP_res_anno[,c("pvalue","log2FoldChange","ID")],by="ID")
aa$score <- aa$log2FoldChange.x * aa$log2FoldChange.y
aa <- aa[order(aa$score),]
aa <- aa[which(aa$ID %in% candidate_IDs),]

candidate_IDs_exp_order <- candidate_IDs_exp[aa$ID,]

candidate_IDs_exp_order$peak_ID <- gsub('ID_','',rownames(candidate_IDs_exp_order))
library(tidyr)
tmp <- do.call(rbind, strsplit(candidate_IDs_exp_order$peak_ID, "_"))
candidate_IDs_exp_order$chr   <- tmp[, 1]
candidate_IDs_exp_order$start <- as.integer(tmp[, 2])
candidate_IDs_exp_order$end   <- as.integer(tmp[, 3])

candidate_IDs_exp_order_gr <- GRanges(candidate_IDs_exp_order)
candidate_IDs_exp_order_gr <- annotatePeak(candidate_IDs_exp_order_gr, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db",verbose=FALSE)
candidate_IDs_exp_order_ <- as.data.frame(candidate_IDs_exp_order_gr)

library(dplyr)

candidate_IDs_exp_order_$baseMean <- rowMeans(candidate_IDs_exp_order_[,6:11])

df_top <- candidate_IDs_exp_order_ %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = baseMean, n = 1, with_ties = FALSE) %>%
  ungroup()
df_top <- as.data.frame(df_top)

library(ComplexHeatmap)
library(circlize)

df <- df_top

expr_cols <- grep("^shRen_", colnames(df), value = TRUE)
mat <- as.matrix(df[, expr_cols])
rownames(mat) <- df$peak_ID

mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0
mat_z[mat_z < -2] <- -2
mat_z[mat_z > 2] <- 2
mat_z <- mat_z[,c(5:6,1:2,3:4)]

genes_interest <- c("Ets1","Batf","Meis2","Spib","Lef1","Il17rb","Zfp36","Csf1r","Spi1","Gfi1b","Mafb")

idx <- which(df$SYMBOL %in% genes_interest)

# 同一基因多个peak时，标签加序号
dup_id <- ave(df$SYMBOL[idx], df$SYMBOL[idx], FUN = seq_along)
labs <- paste0(df$SYMBOL[idx], "_peak", dup_id)

ha_mark <- rowAnnotation(
  mark = anno_mark(
    at = idx,
    labels = labs,
    side = "right"
  )
)

ht <- Heatmap(
  mat_z,
  name = "z-score",
  col = colorRampPalette(jdb_palette("brewer_yes")[1:length(jdb_palette("brewer_spectra"))])(100),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  right_annotation = ha_mark
)

draw(ht)


```


#### Fig4Q
```


GO_up <- read.csv("./GO_ol_peaks.6samples.fc1p0.05.csv",row.names=1)
GO_up$log10p.adjust <- -log10(GO_up$p.adjust)
stem_diff <- c("leukocyte migration","myeloid leukocyte differentiation","myeloid leukocyte migration","hematopoietic stem cell migration","regulation of leukocyte differentiation","stem cell differentiation")

p <- ggbarplot(GO_up[which(GO_up$Description %in% stem_diff),],x="Description",y="log10p.adjust",orientation = "horiz",color = "red",
  fill="red",sort.val = "asc",title="differentiation related pathways")


```





