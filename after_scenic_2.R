# 加载必要的库
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
  library(ggplot2)
  library(AUCell)
  library(SCENIC)
  library(KernSmooth)
  library(RColorBrewer)
  library(plotly)
  library(BiocParallel)
  library(grid)
  library(ComplexHeatmap)
  library(data.table)
  library(scRNAseq)
  library(stringr)
  library(circlize)
  library(reshape2)
})

# 1. 提取 out_SCENIC.loom 信息
path <- "D:/Windows/Share/FCAdata/"
loom_file <- paste0(path, "s_fca_biohub_body_10x.loom")

# 连接loom文件并提取信息
ds <- Connect(filename = loom_file, mode = "r")
regulons_incidMat <- get_regulons(ds, column.attr.name = "MotifRegulonGeneOccurrences")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(ds, column.attr.name = 'MotifRegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(ds)
embeddings <- get_embeddings(ds)
close_loom(ds)

# 查看部分结果
print(regulons_incidMat[1:40, 1:40])
print(tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))]))
print(rownames(regulonAUC))
print(names(regulons))