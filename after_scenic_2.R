suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(SeuratObject)
  library(SeuratData)
  library(loomR)
  library(hdf5r)
  library(LoomExperiment)
  library(SCopeLoomR)
  library(patchwork)
  library(dplyr)
  library(purrr)
  library(graphics)
  library(grDevices)
  library(utils)
  library(datasets)
  library(methods)
  library(base)
  library(stats)
  library(tidyverse)
  library(ggplot2)
})
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
library(AUCell)
library(reshape2)
#### 1.提取 out_SCENIC.loom 信息
path = "D:/Windows/Share/FCAdata/"
ds <-
  Connect(filename = paste0(path, "s_fca_biohub_body_10x.loom"), mode = "r")
loom <- open_loom('out_SCENIC.loom') 

regulons_incidMat <- get_regulons(ds, column.attr.name="MotifRegulonGeneOccurrences")
regulons_incidMat[1:40,1:40] 

regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(ds,column.attr.name='MotifRegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(ds)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(ds)  
close_loom(ds)

rownames(regulonAUC)
names(regulons)



