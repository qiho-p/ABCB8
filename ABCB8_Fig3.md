# ABCB8_Fig3

In this pipeline, we showed the detail codes in Fig3.

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
#### Fig3G
```
shA_all <- read.csv(row.names=1,file='/mnt/data/user_data/yiman/workshop/RNAseq/Abcb8_zjn/FAS_RNASeq_20240813/DEG/shAbcb8.FAS_vs_Veh.all.csv')
mmu.aml.rnk <- shA_all[order(shA_all$log2FoldChange,decreasing=T),]
mmu.aml.rnk$SYMBOL <- rownames(mmu.aml.rnk)
mmu.aml.rnk <- mmu.aml.rnk[,c("SYMBOL","log2FoldChange")]
write.table(mmu.aml.rnk,"/mnt/data/user_data/yiman/workshop/RNAseq/Abcb8_zjn/FAS_RNASeq_20240813/DEG/shA_FAS_vs_shA_Veh.rnk.txt",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)

shA_all <- read.csv(row.names=1,file='/mnt/data/user_data/yiman/workshop/RNAseq/Abcb8_zjn/new_RNASeq_20231228/shABCB8_vs_shRen_all.csv')
shA_vs_shR_p <- shA_all[which(shA_all$pvalue < 0.05),]
shA_vs_shR_DN_400 <- shA_vs_shR_p[order(shA_vs_shR_p$log2FoldChange,decreasing=FALSE),][1:400,]

##Then, Run GSEA prerank.

```
<img src="./Fig3_images/enplot_top_400shA_vs_shR_down.png" alt="enplot_top_400shA_vs_shR_down.png" width="300" />

#### Fig3H
```

shA_all <- read.csv(row.names=1,file='/mnt/data/user_data/yiman/workshop/RNAseq/Abcb8_zjn/FAS_RNASeq_20240813/DEG/shAbcb8.FAS_vs_Veh.all.csv')

shA_all$entrez <- mapIds(x = org.Mm.eg.db,
                        keys = rownames(shA_all),
            keytype ="SYMBOL",
            column ="ENTREZID",
            multiVals="first")

shA_all$group <- as.factor(ifelse(abs(shA_all$log2FoldChange) >= 0.5 & shA_all$pvalue < 0.05, ifelse(shA_all$log2FoldChange >= 0.5 ,'Up','Down'),'NS'))
shA_up <- shA_all[which(shA_all$group == "Up"),]
nrow(shA_up)

GO_up <- enrichGO(gene = shA_up$entrez, 
             OrgDb = org.Mm.eg.db,
        ont = "all", 
                 pvalueCutoff = 0.05, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 1,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GO_up,"/mnt/data/user_data/yiman/workshop/RNAseq/Abcb8_zjn/FAS_RNASeq_20240813/DEG/GO_shA.FAS_vs_shA.Veh_up.csv")

selected_ptw <- c("regulation of hemopoiesis","positive regulation of hemopoiesis","myeloid leukocyte activation","positive regulation of leukocyte differentiation","myeloid cell differentiation","regulation of leukocyte differentiation","T cell differentiation","lymphocyte differentiation","regulation of myeloid cell differentiation","myeloid leukocyte differentiation","positive regulation of myeloid cell differentiation","myeloid cell homeostasis","positive regulation of myeloid leukocyte differentiation","regulation of myeloid leukocyte differentiation")
GO_up <- read.csv(row.names=1,"/mnt/data/user_data/yiman/workshop/RNAseq/Abcb8_zjn/FAS_RNASeq_20240813/DEG/GO_shA.FAS_vs_shA.Veh_up.csv")
GO_up$log10_padj <- -log10(GO_up$p.adjust)

tmp <- GO_up[which(GO_up$Description %in% selected_ptw),]

p1 <- ggbarplot(tmp, x = "Description", y = "log10_padj",
          fill = "ONTOLOGY",           # change fill color by mpg_level
          sort.val = "asc",          # Sort the value in descending order
          ylab = "-log10(padj)",
          legend.title = "Enrichment Group",
          rotate = TRUE,
          ggtheme = theme_pubr(),
          palette="nejm",
          title = "GO BP: shA.FAS_vs_shA.Veh UP"
          )

ggsave(height=7,width=8,"/mnt/data/user_data/yiman/workshop/RNAseq/Abcb8_zjn/github_code/GO.shA.FAS_vs_shA.Veh_up.png")


```
<img src="./Fig3_images/GO.shA.FAS_vs_shA.Veh_up.png" alt="GO.shA.FAS_vs_shA.Veh_up.png" width="400" />
