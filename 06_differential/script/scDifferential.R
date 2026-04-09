#!/usr/bin/env Rscript

# ==========================================================
# 脚本名称：scDifferential.R
# 脚本功能：单细胞差异分析增强模块 (支持 Global DE 和 Subtype Pairwise DE)
# ==========================================================

suppressMessages(library(argparse))
suppressMessages(library(Seurat))
suppressMessages(library(yaml))
suppressMessages(library(this.path))

parser <- ArgumentParser(description = 'scMax 差异分析增强模块: 支持全标签组差异(Global DE)和亚群内组间比较(Group DE)')

parser$add_argument('-i', '--rds', dest='inputrds', required=TRUE, help='输入的 Seurat RDS 文件')
parser$add_argument('-o', '--outdir', dest='outdir', required=TRUE, help='结果输出总目录')
parser$add_argument('-m', '--method', dest='method', default='all_markers', choices=c('all_markers', 'group_comparison'), help='差异模式: all_markers (方案A: 找指定列的所有组Marker), group_comparison (方案B: 在亚群内部进行组间两两比较)')

parser$add_argument('-g', '--groupby', dest='groupby', default='CellType', help='差异比较的分类列 (如: CellType 或 Group)')
parser$add_argument('-s', '--split_by', dest='split_by', default='CellType', help='[模式B专用] 亚群拆分依据列 (通常是 CellType 或 seurat_clusters)')
parser$add_argument('-p', '--cmp_file', dest='cmp_file', default='', help='[模式B专用] Tab分隔的比较对列表文件 (每行两列: GroupA GroupB)')

parser$add_argument('--assay', dest='assay', default='RNA', help='分析所用的测量矩阵 (默认 RNA)')
parser$add_argument('--orgdb', dest='orgdb', default='org.Hs.eg.db', help='功能富集分析的 OrgDb 库')
parser$add_argument('--organism_kegg', dest='organism_kegg', default='hsa', help='功能富集分析的 KEGG 物种简称')
parser$add_argument('--do_enrich', dest='do_enrich', action='store_true', help='是否在差异分析后执行 GO/KEGG 全自动化富集')

opt <- parser$parse_args()

# --- 1. 环境准备与函数库加载 ---
script_dir <- this.path::this.dir()

# 加载 scMax 核心公共库 (复用富集与绘图底层逻辑)
# 寻找路径: 向上跳一级到 05_celltype/script/
func_common <- file.path(dirname(script_dir), "05_celltype", "script", "func_scRNA_celltype_anno.R")
if (file.exists(func_common)) {
  source(func_common)
} else {
  cat("! 警告: 未能找到公共函数库 func_scRNA_celltype_anno.R，部分富集功能可能受限。\n")
}

# 加载本模块私有核心库
func_private <- file.path(script_dir, "func_Differential.R")
if (file.exists(func_private)) {
  source(func_private)
} else {
  stop("错误: 缺失本模块私有核心库 func_Differential.R，请检查安装路径!")
}

# --- 2. 载入数据对象 ---
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
}

cat("\n========================================\n")
cat("=> [1] 正在加载 Seurat 数据对象...\n")
if (!file.exists(opt$inputrds)) stop(paste0("错误: 输入文件不存在 -> ", opt$inputrds))
data <- readRDS(opt$inputrds)

# --- 3. 执行差异分析业务逻辑 ---
if (opt$method == "all_markers") {
  # [模式 A]: 寻找所有组别的 Marker 基因
  # 往往用于描述不同 CellType 之间的个性特征
  cluster_de_all(data, assay=opt$assay, groupby=opt$groupby, outdir=opt$outdir)
  
  if (isTRUE(opt$do_enrich)) {
    cluster_enrich_auto(derds = file.path(opt$outdir, "DEGs.rds"), 
                        orgdb = opt$orgdb, 
                        organism_kegg = opt$organism_kegg, 
                        outdir = opt$outdir)
  }
  
} else if (opt$method == "group_comparison") {
  # [模式 B]: 在各个亚群(CellType)内部进行两两组间比较
  # 往往用于比较同一类细胞在 Case 和 Control 组的不同响应
  if (opt$cmp_file == "" || !file.exists(opt$cmp_file)) {
    stop("错误: 在 'group_comparison' 模式下，必须通过 --cmp_file 提供有效的比较对映射表!")
  }
  
  group_de_by_subtype(data, 
                      split_by = opt$split_by, 
                      condition = opt$groupby, 
                      cmp_file = opt$cmp_file, 
                      outdir = opt$outdir, 
                      assay = opt$assay)
  
  # 注意：模式 B 产生的差异文件较为分散(每个比较对/每个亚群一个目录)，
  # 如果需要执行富集，建议在此处增加对输出目录的遍历调用或由用户后续手动触发。
}

cat("\n========================================\n")
cat("=> 06_differential (差异挖掘模块) 运行结束！\n")
cat("=> 结果报告已存入:", opt$outdir, "\n")
