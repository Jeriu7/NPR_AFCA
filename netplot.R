# 加载所需的包
library(igraph)
library(ggraph)
library(tidygraph)

#### importance 作图 ####
# 读取 TF-gene 数据
# tf_gene_data <- read.csv("path_to_file.csv")
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

### 设置阈值对比importance
for (file in files) {
  
  # 读取数据
  tf_gene_data <- read.table(file, header = TRUE)
  
  # 筛选 NP 和 NPR 的目标基因
  tf_gene_data_NP <- tf_gene_data[tf_gene_data$target %in% NP, ]
  tf_gene_data_NPR <- tf_gene_data[tf_gene_data$target %in% NPR, ]

  # 合并选取的 NP 和 NPR 数据
  tf_gene_data_selected <- rbind(tf_gene_data_NP, tf_gene_data_NPR)
  tf_gene_data_selected <- tf_gene_data_selected[order(-as.numeric(tf_gene_data_selected$importance)),]
  
  # 合并选取的 NP 和 NPR 数据
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






######### ####
library(ggplot2)

# 绘制 TF 对 NP 和 NPR 的调控重要性的柱状图
ggplot(importance_summary, aes(x = TF, y = importance, fill = gene_type)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "TF对NP和NPR的调控重要性对比", x = "TF", y = "平均调控重要性") +
  theme_minimal()

# 创建有向图对象
tf_g <- tf_gene_data[(as.numeric(tf_gene_data$importance)>50),]
g <- graph_from_data_frame(tf_g, directed = TRUE)

# 使用 ggraph 绘制网络
ggraph(g, layout = 'fr') + 
  geom_edge_link(aes(width = tf_g$importance), alpha = 0.8) + 
  geom_node_point(size = 5, color = "lightblue") + 
  geom_node_text(aes(label = name), vjust = 1, hjust = 1) + 
  theme_void()






regulons_incidMat <- read_csv("scenic_results/out.afca_body_seurat_n_30.csv")
regulons_incidMat[1:4,1:4] 
regulons_incidMat_head[1:4,1:4] 
regulons_body_30 <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom_body_30,column.attr.name='MotifRegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom_body_30)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom_body_30)  
close_loom(loom_body_30)

rownames(regulonAUC)
names(regulons_body_30)
## Body
loom_body <-
  Connect(filename = paste0(path,"s_fca_biohub_body_10x.loom"), mode = "r")
regulons_incidMat <- get_regulons(loom_body, column.attr.name="MotifRegulonGeneOccurrences")
regulons_incidMat[1:4,1:4] 
regulons_body <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom_body,column.attr.name='MotifRegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom_body)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom_body)  
close_loom(loom_body)

rownames(regulonAUC)
names(regulons_body)






























# grn <- read.table(ingrn, sep='\t', header=T, stringsAsFactors=F)
# inregulons1 = gsub('[(+)]', '', inregulons)
# c1 <- which(grn$TF %in% inregulons1)
# grn <- grn[c1,]
# 
# pdf(paste0(ntop2, '_regulon_netplot.pdf'), 10, 10)
# for (tf in unique(grn$TF)) {
#   tmp <- subset(grn, TF == tf)
#   if (dim(tmp)[1] > ntop2) {
#     tmp <- tmp[order(tmp$importance, decreasing=T),]
#     tmp <- tmp[1:ntop2,]
#   }
#   
#   node2 <- data.frame(tmp$target)
#   node2$node.size = 1.5
#   node2$node.colour = 'black'
#   colnames(node2) <- c('node', 'node.size', 'node.colour')
#   df1 <- data.frame(node = tf, node.size = 2, node.colour = '#FFDA00')
#   node2 <- rbind(df1, node2)
#   
#   edge2 <- tmp
#   colnames(edge2) <- c('from', 'to', 'edge.width')
#   edge2$edge.colour <- "#1B9E77"
#   torange = c(0.1, 1)
#   edge2$edge.width <- scales::rescale(edge2$edge.width, to=torange)
#   
#   graph_data <- tidygraph::tbl_graph(nodes = node2, edges = edge2, directed = T)
#   p1 <- ggraph(graph = graph_data, layout = "stress", circular = TRUE) + 
#     geom_edge_arc(aes(edge_colour = edge.colour, edge_width = edge.width)) +
#     scale_edge_width_continuous(range = c(1, 0.2)) +
#     geom_node_point(aes(colour = node.colour, size = node.size)) + 
#     theme_void() +
#     geom_node_label(aes(label = node, colour = node.colour), size = 3.5, repel = TRUE)
#   
#   p1 <- p1 + scale_color_manual(values = c('#FFDA00', 'black')) +
#     scale_edge_color_manual(values = c("#1B9E77"))
#   
#   print(p1)
# }
# dev.off()


grn <- read.table("scenic_results/adj.afca_body_seurat_n_30.tsv",sep='\t',header=T,stringsAsFactors=F)
inregulons1=gsub('[(+)]','',inregulons)
c1 <- which(grn$TF %in% inregulons1)
grn <- grn[c1,]
#edge1 <- data.frame()
#node1 <- data.frame()
pdf(paste0(ntop2,'_regulon_netplot.pdf'),10,10)
for (tf in unique(grn$TF)) {
  tmp <- subset(grn,TF==tf)
  if (dim(tmp)[1] > ntop2) {
    tmp <- tmp[order(tmp$importance,decreasing=T),]
    tmp <- tmp[1:ntop2,]
  }
}
  