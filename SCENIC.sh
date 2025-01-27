# 使用 PySCENIC 进行基因调控网络 (GRN) 推断
# 输入文件: sample1.loom (单细胞数据), all_TFs.txt (转录因子列表)
# 输出文件: adj.sample.tsv (基因调控网络)
# 使用 GRNBoost2 方法，并行工作线程数为 20
pyscenic grn \
  --num_workers 20 \
  --output adj.sample.tsv \
  --method grnboost2 \
  sample1.loom \
  all_TFs.txt

# 使用 PySCENIC 进行调控子上下文分析
# 输入文件: adj.sample.tsv (基因调控网络), dm6_v10_clust.genes_vs_motifs.rankings.feather (基因组数据库)
# 输入文件: motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl (motif 注释文件), sample1.loom (单细胞数据)
# 输出文件: reg.csv (调控子结果)
# 使用 dask_multiprocessing 模式，并行工作线程数为 3，并屏蔽 dropout 事件
pyscenic ctx \
  adj.sample.tsv \
  dm6_v10_clust.genes_vs_motifs.rankings.feather \
  --annotations_fname motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl \
  --expression_mtx_fname sample1.loom \
  --mode "dask_multiprocessing" \
  --output reg.csv \
  --num_workers 3 \
  --mask_dropouts

# 使用 PySCENIC 进行 AUCell 分析
# 输入文件: sample1.loom (单细胞数据), reg.csv (调控子结果)
# 输出文件: sample1_SCENIC.loom (包含 AUCell 分析结果的 loom 文件)
# 并行工作线程数为 3
pyscenic aucell \
  sample1.loom \
  reg.csv \
  --output sample1_SCENIC.loom \
  --num_workers 3