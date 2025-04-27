# ABCB8_Fig1

In this pipeline, we showed the detail codes in Fig1.


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
#### Fig1F
```
bao_myeloid <- mcreadRDS("/mnt/data/user_data/yiman/project/zjn_abcb8_screening/Bao.patient.normal.BM.myeloid.final.rds", mc.cores = 20)
p <- DimPlot(object = bao_myeloid, reduction = "umap",repel=TRUE,label=FALSE,group.by="CellType",raster = FALSE) + ggtitle("Normal_BM")
ggsave(height=3.5,width=5,"/mnt/data/user_data/yiman/project/zjn_abcb8_screening/github_code/Dimplot_Bao.patient.normal.BM.myeloid.final.png")

```
![Fig1F_Dimplot](.\Fig1_images\Dimplot_Bao.patient.normal.BM.myeloid.final.png)

#### Fig1G
```

data_tmp <- as.matrix(GetAssayData(object = bao_myeloid, slot = "counts",assay="RNA"))
results <- CytoTRACE(data_tmp, ncores = 20, subsamplesize = 1000)
bao_myeloid$CytoTRACE <- results$CytoTRACE[rownames(bao_myeloid[[]])]
bao_myeloid$CytoTRACErank <- results$CytoTRACErank[rownames(bao_myeloid[[]])]

AML.exp <- FetchData(object = bao_myeloid, vars = c("CytoTRACE","CellType","ABCB8"),slot="data")
AML.exp <- na.omit(AML.exp)
AML.exp$ABCB8 <- as.numeric(AML.exp$ABCB8)
AML.exp$CytoTRACE <- as.numeric(AML.exp$CytoTRACE)

AML.exp <- AML.exp[order(AML.exp$CytoTRACE),]
AML.exp$CytoTRACE_order <- c(1:nrow(AML.exp))

set.seed(6)
AML.exp__ <- AML.exp
CellType <- c("Erythroid","Mono&Macro","HSPC","GMP","Ery.Pro","Megak.")
All_gsva_seura_ <- future_lapply(CellType,function(x) {
    sel_tmp <- AML.exp__[which(AML.exp__$CellType == x),]
    cell_name <- rownames(sel_tmp)
    cell_name <- sample(cell_name,length(cell_name))
    cell.sets1 <- split(cell_name, sort(1:length(cell_name)%%70))
    profile.set1 = matrix(, nrow = length(cell.sets1), ncol = 3)
      for (i in 1:length(cell.sets1)) {
        this.set <- cell.sets1[[i]]
        sub.matrix <- sel_tmp[this.set, ]
        ABCB8 <- mean(sub.matrix$ABCB8)
        CytoTRACE <- mean(sub.matrix$CytoTRACE)
        profile.set1[i,1 ] <- ABCB8
        profile.set1[i,2 ] <- CytoTRACE
        profile.set1[i,3 ] <- x
      }
    rownames(profile.set1) <- paste(x, 1:length(cell.sets1),sep="_")
    colnames(profile.set1) <- c("ABCB8","CytoTRACE","CellType")
    return(profile.set1)
    }, future.seed = TRUE)
All_gsva_seura <- do.call(rbind,All_gsva_seura_)
All_gsva_seura <- as.data.frame(All_gsva_seura)
All_gsva_seura$ABCB8 <- as.numeric(All_gsva_seura$ABCB8)
All_gsva_seura$CytoTRACE <- as.numeric(All_gsva_seura$CytoTRACE)

myeloid.lin <- All_gsva_seura[which(All_gsva_seura$CellType %in% c("Mono&Macro","HSPC","GMP")),]
myeloid.lin <- myeloid.lin[order(myeloid.lin$CytoTRACE,decreasing=TRUE),]
myeloid.lin$CytoTRACE_order <- c(1:nrow(myeloid.lin))
myeloid.lin[which(myeloid.lin$ABCB8 > 0.1),]$ABCB8 <- 0.1
myeloid.lin$ABCB8.zscore <- (myeloid.lin$ABCB8 - mean(myeloid.lin$ABCB8)) / sd(myeloid.lin$ABCB8)
mcsaveRDS(myeloid.lin,"/mnt/data/user_data/yiman/project/zjn_abcb8_screening/github_code/Bao.patient.normal.BM.Cyto.ABCB8.70.rds")

p <- ggscatter(myeloid.lin, 
  x = "CytoTRACE_order", y = "ABCB8.zscore",color="CellType",
  add = "loess", conf.int = TRUE,palette="npg",
  alpha=0.5, fullrange = TRUE, rug = TRUE,size=2,
  title="ABCB8 exp. and CytoTRACE",add.params=list(color = "red", fill = "lightgray")) + 
stat_cor(color="black", method = c("pearson"))

ggsave(height=5,width=5,"/mnt/data/user_data/yiman/project/zjn_abcb8_screening/github_code/Cyto.ABCB8.Bao.patient.normal.BM.myeloid.70.seed6.png")

```
![Fig1G](.\Fig1_images\Cyto.ABCB8.Bao.patient.normal.BM.myeloid.70.seed6.png)


#### Fig1H

```
table(bao_myeloid$CellType)
SCLC.only.GSVA <- bao_myeloid
Idents(SCLC.only.GSVA) <- SCLC.only.GSVA$CellType
SCLC.only.GSVA$CellType <- factor(SCLC.only.GSVA$CellType,levels=
  c("HSPC","GMP","Ery.Pro","Erythroid","Megak.","Mono&Macro"))
Idents(SCLC.only.GSVA) <- SCLC.only.GSVA$CellType
SCLC.only.GSVA$CellType <- factor(SCLC.only.GSVA$CellType,levels=
  c("HSPC","GMP","Ery.Pro","Erythroid","Megak.","Mono&Macro"))
All_gsva_seura_ <- future_lapply(1:length(levels(SCLC.only.GSVA$CellType)),function(i) {
    sel_tmp <- subset(SCLC.only.GSVA,idents=levels(SCLC.only.GSVA$CellType)[i])
    sel_tmp <- pseudo_bulk_seurat_mean_random(seurat_obj=sel_tmp,num_split=20,seed.use=1,slot="data",prefix=levels(SCLC.only.GSVA$CellType)[i],assay="RNA")
    metadata <- data.frame(cell_type=c(rep(levels(SCLC.only.GSVA$CellType)[i],20)),
    row.names=colnames(sel_tmp))
    sel_gsva_seurat <- CreateSeuratObject(counts = sel_tmp,assay = 'RNA',project = 'RNA',min.cells = 0,meta.data = metadata)
    message(levels(SCLC.only.GSVA$CellType)[i], " is done")
    return(sel_gsva_seurat)
})
All_gsva_seura <- merge(x = All_gsva_seura_[[1]], y = All_gsva_seura_[c(2:length(All_gsva_seura_))])
All_gsva_seura$cell_type <- factor(All_gsva_seura$cell_type,levels=c("HSPC","GMP","Ery.Pro","Erythroid","Megak.","Mono&Macro"))
Idents(All_gsva_seura) <- All_gsva_seura$cell_type
mcsaveRDS(All_gsva_seura,"/mnt/data/user_data/yiman/project/zjn_abcb8_screening/Bao.patient.normal.BM.myeloid.RNA.20.rds")

All_gsva_seura <- mcreadRDS("/mnt/data/user_data/yiman/project/zjn_abcb8_screening/Bao.patient.normal.BM.myeloid.RNA.20.rds")
All_gsva_seura$cell_type <- factor(All_gsva_seura$cell_type,levels=c("HSPC","GMP","Ery.Pro","Erythroid","Megak.","Mono&Macro"))

all_data <- FetchData(object = All_gsva_seura, vars = c("ABCB8","cell_type"),slot="data")
all_data <- all_data[which(all_data$cell_type %in% c("HSPC","GMP","Megak.","Mono&Macro")),]

p <- ggboxplot(all_data,x="cell_type",y="ABCB8",fill="cell_type",group="cell_type",alpha=0.8,
  xlab="Cell type",ylab="ABCB8 normalized counts")+RotatedAxis()+NoLegend()+ylim(0,0.15)

ggsave(height=4,width=4,"/mnt/data/user_data/yiman/project/zjn_abcb8_screening/github_code/boxplot.Bao.patient.normal.BM.myeloid.ABCB8.20.png")

```
![Fig1H](.\Fig1_images\boxplot.Bao.patient.normal.BM.myeloid.ABCB8.20.png)








