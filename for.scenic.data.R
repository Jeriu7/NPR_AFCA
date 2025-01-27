# 加载必要的R包
library(Seurat)  # 单细胞数据分析
library(dplyr)   # 数据操作

# 读取转录因子 (TFs) 列表
TFs <- read.csv("cisTarget_databases/allTFs_dmel.txt", header = FALSE)$V1

# 定义神经肽受体 (NPR) 和神经肽 (NP) 列表
NPR_HR <- c("AkhR", "InR", "AstA-R1", "AstA-R2", "AstC-R1", "AstC-R2", "CapaR", "CCAP-R", 
            "CCHa1-R", "CCHa2-R", "CCKLR-17D1", "CCKLR-17D3", "CG4313", "CG12290", 
            "CG13229", "CG13575", "CG13995", "CG30340", "CG32547", "CG33639", "CNMaR", 
            "CrzR", "ETHR", "FMRFaR", "Lgr1", "Lgr3", "Lgr4", "Lkr", "moody", "MsR1", 
            "MsR2", "NPFR", "PK1-R", "PK2-R1", "PK2-R2", "Proc-R", "rk", "RYa-R", 
            "SIFaR", "sNPF-R", "SPR", "TkR86C", "TkR99D", "Tre1", "TrissinR")

NP <- c("Akh", "amn", "AstA", "AstC", "AstCC", "Burs", "Capa", "CCAP", "CCHa1", "CCHa2", 
        "CNMa", "Crz", "Dh31", "Dh44", "Dsk", "Eh", "ETH", "FMRFa", "Gpa2", "Gpb5", 
        "Hug", "Ilp1", "Ilp2", "Ilp3", "Ilp4", "Ilp5", "Ilp6", "Ilp7", "Ilp8", "ITP", 
        "Lk", "Mip", "Ms", "NPF", "Nplp1", "Nplp2", "Nplp3", "Nplp4", "Orcokinin", 
        "Pburs", "Pdf", "Proc", "Ptth", "RYa", "SIFa", "sNPF", "SP", "spab", "Tk")

# 读取并整合NPR列表
NPR <- read.csv("id_validation_table_NPR-activity.txt", sep = "\t")$current_symbol
NPR <- union(NPR, NPR_HR)

# 整合TFs、NP和NPR基因列表
genes_TF_NP_NPR <- union(TFs, NP)
genes_TF_NP_NPR <- union(genes_TF_NP_NPR, NPR)

# 根据基因列表对Seurat对象进行子集化
afca_head_seurat_n <- subset(afca_head_seurat, features = genes_TF_NP_NPR)
afca_body_seurat_n <- subset(afca_body_seurat, features = genes_TF_NP_NPR)

# 根据年龄筛选细胞
Cells.head.5 <- WhichCells(afca_head_seurat_n, expression = age == '5')
Cells.head.30 <- WhichCells(afca_head_seurat_n, expression = age == '30')
Cells.head.50 <- WhichCells(afca_head_seurat_n, expression = age == '50')
Cells.head.70 <- WhichCells(afca_head_seurat_n, expression = age == '70')

Cells.body.5 <- WhichCells(afca_body_seurat_n, expression = age == '5')
Cells.body.30 <- WhichCells(afca_body_seurat_n, expression = age == '30')
Cells.body.50 <- WhichCells(afca_body_seurat_n, expression = age == '50')
Cells.body.70 <- WhichCells(afca_body_seurat_n, expression = age == '70')

# 导出头部数据
write.csv(t(as.matrix(afca_head_seurat_n[, Cells.head.5]@assays$RNA@counts)), 
          file = "afca_head_seurat_n_5.for.scenic.data.csv")
write.csv(t(as.matrix(afca_head_seurat_n[, Cells.head.30]@assays$RNA@counts)), 
          file = "afca_head_seurat_n_30.for.scenic.data.csv")
write.csv(t(as.matrix(afca_head_seurat_n[, Cells.head.50]@assays$RNA@counts)), 
          file = "afca_head_seurat_n_50.for.scenic.data.csv")
write.csv(t(as.matrix(afca_head_seurat_n[, Cells.head.70]@assays$RNA@counts)), 
          file = "afca_head_seurat_n_70.for.scenic.data.csv")

# 导出身体数据
write.csv(t(as.matrix(afca_body_seurat_n[, Cells.body.5]@assays$RNA@counts)), 
          file = "afca_body_seurat_n_5.for.scenic.data.csv")
write.csv(t(as.matrix(afca_body_seurat_n[, Cells.body.30]@assays$RNA@counts)), 
          file = "afca_body_seurat_n_30.for.scenic.data.csv")
write.csv(t(as.matrix(afca_body_seurat_n[, Cells.body.50]@assays$RNA@counts)), 
          file = "afca_body_seurat_n_50.for.scenic.data.csv")
write.csv(t(as.matrix(afca_body_seurat_n[, Cells.body.70]@assays$RNA@counts)), 
          file = "afca_body_seurat_n_70.for.scenic.data.csv")