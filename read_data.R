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
rm(list = ls())
Convert(
  'adata_body_S_v1.0.h5ad',
  dest = "h5seurat",
  overwrite = TRUE
)
# ds5 <-
#   Connect(file = "adata_headBody_S_v1.0.h5seurat", mode = 'r')

# ### 调用python读入（不全）
# # 导入所需的包和库
# library(anndata)
# library(reticulate)
# # use_condaenv(condaenv = "C:/python/envs/py39-bio", required = TRUE)
# use_python("C:/Users/admin/anaconda3/python.exe")
# ad <- import("anndata")
# 
# # 读取数据集
# final_ad <- ad$read_h5ad("adata_headBody_S_v1.0.h5ad")
# 
# # 设置变量名
# colnames(final_ad$X) <- final_ad$var$features
# 
# # 创建 Seurat 对象
# adata_headBody_seurat <- Seurat::CreateSeuratObject(counts = final_ad$X, assay = "RNA",
#                                          meta.data = final_ad$obs)

### Seurat导入
afca_headBody_seurat <- LoadH5Seurat("adata_headBody_S_v1.0.h5seurat",meta.data = T)
afca_head_seurat <- LoadH5Seurat("adata_head_S_v1.0.h5seurat",meta.data = T)
afca_body_seurat <- LoadH5Seurat("adata_body_S_v1.0.h5seurat",meta.data = T)

saveRDS(afca_headBody_seurat, "afca_headBody_seurat.RDS")
saveRDS(afca_head_seurat, "afca_head_seurat.RDS")
saveRDS(afca_body_seurat, "afca_body_seurat.RDS")























