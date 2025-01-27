# 加载所需的R包
suppressPackageStartupMessages({
  library(Seurat)            # 单细胞数据分析
  library(SeuratDisk)        # Seurat数据格式转换
  library(SeuratObject)      # Seurat对象操作
  library(SeuratData)        # Seurat示例数据集
  library(loomR)             # loom文件操作
  library(hdf5r)             # HDF5文件格式支持
  library(LoomExperiment)    # loom实验数据操作
  library(SCopeLoomR)        # SCope loom文件支持
  library(patchwork)         # 图形拼接
  library(dplyr)             # 数据操作
  library(purrr)             # 函数式编程工具
  library(graphics)          # 图形系统
  library(grDevices)         # 图形设备
  library(utils)             # 实用工具
  library(datasets)          # 示例数据集
  library(methods)           # 方法定义
  library(base)              # 基础函数
  library(stats)             # 统计函数
  library(tidyverse)         # 数据科学工具集
  library(ggplot2)           # 数据可视化
})

# 清空当前环境
rm(list = ls())

# 将H5AD文件转换为H5Seurat格式
Convert(
  'adata_body_S_v1.0.h5ad',
  dest = "h5seurat",
  overwrite = TRUE
)

# 使用Seurat加载H5Seurat文件
afca_headBody_seurat <- LoadH5Seurat("adata_headBody_S_v1.0.h5seurat", meta.data = TRUE)
afca_head_seurat <- LoadH5Seurat("adata_head_S_v1.0.h5seurat", meta.data = TRUE)
afca_body_seurat <- LoadH5Seurat("adata_body_S_v1.0.h5seurat", meta.data = TRUE)

# 保存Seurat对象为RDS文件
saveRDS(afca_headBody_seurat, "afca_headBody_seurat.RDS")
saveRDS(afca_head_seurat, "afca_head_seurat.RDS")
saveRDS(afca_body_seurat, "afca_body_seurat.RDS")