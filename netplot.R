# 加载所需的包
library(igraph)      # 网络分析
library(ggraph)      # 网络可视化
library(tidygraph)   # 图形数据处理
library(ggplot2)     # 数据可视化
library(dplyr)       # 数据操作

#### 1. 读取 TF-gene 数据并计算调控重要性 ####
# 定义文件路径
files <- c("scenic_results/adj.afca_body_seurat_n_5.tsv",
           "scenic_results/adj.afca_head_seurat_n_5.tsv",
           "scenic_results/adj.afca_body_seurat_n_30.tsv",
           "scenic_results/adj.afca_head_seurat_n_30.tsv",
           "scenic_results/adj.afca_body_seurat_n_50.tsv",
           "scenic_results/adj.afca_head_seurat_n_50.tsv",
           "scenic_results/adj.afca_body_seurat_n_70.tsv",
           "scenic_results/adj.afca_head_seurat_n_70.tsv")

# 定义 NP 和 NPR 列表
NPR_HR <- c("AkhR", "InR", "AstA-R1", "AstA-R2", "AstC-R1", "AstC-R2", "CapaR", "CCAP-R", 
            "CCHa1-R", "CCHa2-R", "CCKLR-17D1", "CCKLR-17D3", "CG4313", "CG12290", "CG13229", 
            "CG13575", "CG13995", "CG30340", "CG32547", "CG33639", "CNMaR", "CrzR", "ETHR", 
            "FMRFaR", "Lgr1", "Lgr3", "Lgr4", "Lkr", "moody", "MsR1", "MsR2", "NPFR", "PK1-R", 
            "PK2-R1", "PK2-R2", "Proc-R", "rk", "RYa-R", "SIFaR", "sNPF-R", "SPR", "TkR86C", 
            "TkR99D", "Tre1", "TrissinR")
NP <- c("Akh", "amn", "AstA", "AstC", "AstCC", "Burs", "Capa", "CCAP", "CCHa1", "CCHa2", 
        "CNMa", "Crz", "Dh31", "Dh44", "Dsk", "Eh", "ETH", "FMRFa", "Gpa2", "Gpb5", "Hug", 
        "Ilp1", "Ilp2", "Ilp3", "Ilp4", "Ilp5", "Ilp6", "Ilp7", "Ilp8", "ITP", "Lk", "Mip", 
        "Ms", "NPF", "Nplp1", "Nplp2", "Nplp3", "Nplp4", "Orcokinin", "Pburs", "Pdf", "Proc", 
        "Ptth", "RYa", "SIFa", "sNPF", "SP", "spab", "Tk")

# 从文件中读取 NPR 列表
NPR <- read.csv("id_validation_table_NPR-activity.txt", sep = "\t")
NPR <- NPR$current_symbol
NPR <- union(NPR, NPR_HR)

# 创建保存结果的目录
dir.create("ave_importance", showWarnings = FALSE)

# 循环处理每个文件
for (file in files) {
  # 读取数据
  tf_gene_data <- read.table(file, header = TRUE)
  
  # 筛选 NP 和 NPR 的目标基因
  tf_gene_data_NP <- tf_gene_data[tf_gene_data$target %in% NP, ]
  tf_gene_data_NPR <- tf_gene_data[tf_gene_data$target %in% NPR, ]
  
  # 分别选择前 50 个 NP 和前 50 个 NPR
  tf_gene_data_NP <- tf_gene_data_NP[1:50, ]
  tf_gene_data_NPR <- tf_gene_data_NPR[1:50, ]
  
  # 合并选取的 NP 和 NPR 数据
  tf_gene_data_selected <- rbind(tf_gene_data_NP, tf_gene_data_NPR)
  
  # 添加 NP 和 NPR 的标记
  tf_gene_data_selected$gene_type <- ifelse(tf_gene_data_selected$target %in% NP, "NP", 
                                            ifelse(tf_gene_data_selected$target %in% NPR, "NPR", NA))
  
  # 计算每个 TF 对 NP 和 NPR 的平均调控重要性
  importance_summary <- aggregate(importance ~ TF + gene_type, data = tf_gene_data_selected, FUN = mean)
  
  # 打印结果
  print(importance_summary)
  
  # 提取文件名以用于保存
  file_name <- gsub("scenic_results/", "", file)
  file_name <- gsub(".tsv", ".importance_summary.csv", file_name)
  
  # 保存结果
  write.csv(importance_summary, paste0("ave_importance/", file_name), row.names = FALSE)
}

#### 2. 设置阈值对比调控重要性 ####
for (file in files) {
  # 读取数据
  tf_gene_data <- read.table(file, header = TRUE)
  
  # 筛选 NP 和 NPR 的目标基因
  tf_gene_data_NP <- tf_gene_data[tf_gene_data$target %in% NP, ]
  tf_gene_data_NPR <- tf_gene_data[tf_gene_data$target %in% NPR, ]
  
  # 合并选取的 NP 和 NPR 数据
  tf_gene_data_selected <- rbind(tf_gene_data_NP, tf_gene_data_NPR)
  tf_gene_data_selected <- tf_gene_data_selected[order(-as.numeric(tf_gene_data_selected$importance)), ]
  
  # 选择前 100 个调控关系
  tf_gene_data_selected <- tf_gene_data_selected[1:100, ]
  
  # 添加 NP 和 NPR 的标记
  tf_gene_data_selected$gene_type <- ifelse(tf_gene_data_selected$target %in% NP, "NP", 
                                            ifelse(tf_gene_data_selected$target %in% NPR, "NPR", NA))
  
  # 计算每个 TF 对 NP 和 NPR 的平均调控重要性
  importance_summary <- aggregate(importance ~ TF + gene_type, data = tf_gene_data_selected, FUN = mean)
  
  # 打印结果
  print(importance_summary)
  
  # 提取文件名以用于保存
  file_name <- gsub("scenic_results/", "", file)
  file_name <- gsub(".tsv", ".importance_summary.csv", file_name)
  
  # 保存结果
  write.csv(importance_summary, paste0("ave_importance/Top100_", file_name), row.names = FALSE)
}

#### 3. 可视化 TF 对 NP 和 NPR 的调控重要性 ####
# 绘制柱状图
ggplot(importance_summary, aes(x = TF, y = importance, fill = gene_type)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "TF对NP和NPR的调控重要性对比", x = "TF", y = "平均调控重要性") +
  theme_minimal()

# 创建有向图对象
tf_g <- tf_gene_data[as.numeric(tf_gene_data$importance) > 50, ]
g <- graph_from_data_frame(tf_g, directed = TRUE)

# 使用 ggraph 绘制网络
ggraph(g, layout = 'fr') + 
  geom_edge_link(aes(width = tf_g$importance), alpha = 0.8) + 
  geom_node_point(size = 5, color = "lightblue") + 
  geom_node_text(aes(label = name), vjust = 1, hjust = 1) + 
  theme_void()

#### 4. 读取 regulons 数据并分析 ####
# 读取 regulons 数据
regulons_incidMat <- read.csv("scenic_results/out.afca_body_seurat_n_30.csv")
regulons_incidMat[1:4, 1:4]

# 将 regulons 转换为基因列表
regulons_body_30 <- regulonsToGeneLists(regulons_incidMat)

# 获取 regulons AUC 和阈值
regulonAUC <- get_regulons_AUC(loom_body_30, column.attr.name = 'MotifRegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom_body_30)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

# 获取嵌入信息并关闭 loom 文件
embeddings <- get_embeddings(loom_body_30)
close_loom(loom_body_30)

# 打印 regulons AUC 的行名和 regulons 名称
rownames(regulonAUC)
names(regulons_body_30)

#### 5. 分析 Body 数据 ####
# 连接 loom 文件
loom_body <- Connect(filename = paste0(path, "s_fca_biohub_body_10x.loom"), mode = "r")

# 获取 regulons 数据
regulons_incidMat <- get_regulons(loom_body, column.attr.name = "MotifRegulonGeneOccurrences")
regulons_incidMat[1:4, 1:4]

# 将 regulons 转换为基因列表
regulons_body <- regulonsToGeneLists(regulons_incidMat)

# 获取 regulons AUC 和阈值
regulonAUC <- get_regulons_AUC(loom_body, column.attr.name = 'MotifRegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom_body)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

# 获取嵌入信息并关闭 loom 文件
embeddings <- get_embeddings(loom_body)
close_loom(loom_body)

# 打印 regulons AUC 的行名和 regulons 名称
rownames(regulonAUC)
names(regulons_body)