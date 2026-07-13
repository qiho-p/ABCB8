# ABCB8_Fig7

In this pipeline, we showed the detail codes in Fig7.

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


#### Fig7A

```r



computeMatrix reference-point --referencePoint center -b 3000 -a 3000 \
-R peaks.H3k27me3.final.3vs3.sample.bed \
-S ABCB8.high.3.samples.sorted.bw ABCB8.low.3.samples.sorted.bw \
--numberOfProcessors 30 --skipZeros \
-o peaks.H3k27me3.merge.mat.gz 

plotHeatmap -m peaks.H3k27me3.merge.mat.gz \
 -out peaks.H3k27me3.merge.pdf \
 --colorList 'white, green' 'white, green' 'white, green' 'white, green' 'white, green' 'white, green' 'white, green' 'white, green' 'white, green' 'white, green' 'white, green' \
 --plotFileFormat "pdf" \
 --missingDataColor "white"


```


#### Fig7C & 7D

```r


all.counts <- read.csv(row.names=1,"BeatAML.7C.data.csv")

p <- ggplot(data=all.counts,aes(x=ABCB8.zscore,y=UTX.sig.final.zscore)) 
p <- p + stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon')+
  scale_fill_viridis_c(option = "F")+
  geom_smooth(method=lm,colour="#cb181d",se=TRUE)+
  guides(alpha="none")+
  geom_point(color="black",size=1) + theme_pubr()+
  stat_cor(color="black", method = c("spearman")) + ylim(-2.7,2.7)

final_table <- read.csv(row.names=1,"GSE67040.7D.data.csv")

p <- ggplot(data=final_table,aes(x=ABCB8.zscore,y=UTX_binding.zscore)) 
p <- p + stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon')+
  scale_fill_viridis_c(option = "F")+
  geom_smooth(method=lm,colour="#cb181d",se=TRUE)+
  guides(alpha="none")+
  geom_point(color="black",size=1) + theme_pubr()+
  stat_cor(color="black", method = c("spearman")) + ggtitle("GSE67040")



```


#### FigS7A

```r

dds_AML_normalize <- read.csv(row.names=1,file='./Nature_2022/outs/DESeq2.norm.AML.reorder.S7A.csv')

p <- ggboxplot(dds_AML_normalize,x="group",y="ABCB8",fill="group",palette="lancet",add="jitter")+RotatedAxis()+
  stat_compare_means(comparisons=list(c("High","Low")),method="t.test")


```

