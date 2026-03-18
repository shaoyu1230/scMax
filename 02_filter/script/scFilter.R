#!/usr/bin/env Rscript

# ==========================================================
# 脚本名称：scFilter.R
# 脚本功能：基于 scQC 质量评估结果，对单细胞数据进行阈值过滤
# ==========================================================

suppressMessages(library(argparse))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(SCP))

parser <- argparse::ArgumentParser(description = '单细胞转录组数据过滤脚本')

parser$add_argument('--rds', dest='inputrds', required=TRUE, help='输入已完成质控打分的 Seurat RDS 格式文件', type='character')
parser$add_argument('--meta', dest='metafile', default="", help='可选的新 metadata 补充文件 (CSV)', type='character')
parser$add_argument('--outdir', dest='outdir', required=TRUE, help='输出结果存储目录')
parser$add_argument('--col_sample', dest='col_sample', default="Sample", help='样本所在原属性列')

parser$add_argument('--min_feat', dest='min_feat', default=200, type='integer', help='基因表达数目(nFeature_RNA)最小值')
parser$add_argument('--max_feat', dest='max_feat', default=1000000000, type='integer', help='基因表达数目(nFeature_RNA)最大值')
parser$add_argument('--min_count', dest='min_count', default=200, type='integer', help='转录本表达量(nCount_RNA)最小值')
parser$add_argument('--max_count', dest='max_count', default=1000000000, type='integer', help='转录本表达量(nCount_RNA)最大值')
parser$add_argument('--max_mt', dest='max_mt', default=20.0, type='double', help='线粒体基因比例(percent.mt)最大值')
parser$add_argument('--max_hb', dest='max_hb', default=5.0, type='double', help='血红蛋白基因比例(percent.hb)最大值')
parser$add_argument('--rm_doublet', dest='rm_doublet', action='store_true', default=FALSE, help='是否根据 DF.classifications 剔除双细胞')
parser$add_argument('--max_decontX', dest='max_decontX', default=1.0, type='double', help='DecontX 游离污染(Contamination)最大容忍度')

opt <- parser$parse_args()

outdir <- opt$outdir
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

cat("=> [1] 读取输入数据...\n")
if (!file.exists(opt$inputrds)) {
  stop(paste0("输入文件 ", opt$inputrds, " 不存在!"))
}

object <- readRDS(opt$inputrds)

if (opt$metafile != "") {
  new_meta <- read.csv(opt$metafile, row.names = 1, stringsAsFactors = FALSE)
  if (!identical(sort(rownames(new_meta)), sort(rownames(object@meta.data)))){
    stop("提供的 Metadata 数据与对象细胞无法对齐，停止分析程序。")
  }
  new_meta <- new_meta[rownames(object@meta.data), , drop = FALSE]
  object@meta.data <- new_meta
}

# 检查必要的列
req_cols <- c("nCount_RNA", "nFeature_RNA")
for(c in req_cols) {
  if (!(c %in% names(object@meta.data))) {
    stop(paste0("必要的数据列 ", c, " 不在您的 RDS 中！请确认是否经过上一步 scQC 质控！"))
  }
}

cat("=> [2] 开始执行阈值过滤...\n")
A <- opt$min_feat
B <- opt$max_feat
CC <- opt$min_count
DD <- opt$max_count
mt_percent <- opt$max_mt
HB_percent <- opt$max_hb
rm_doub <- opt$rm_doublet
decontX <- opt$max_decontX
col_sample <- opt$col_sample

object1 <- subset(object, subset = nFeature_RNA > A & nFeature_RNA < B & nCount_RNA > CC & nCount_RNA < DD)

if ("percent.mt" %in% names(object@meta.data)) {
  object1 <- subset(object1, subset = percent.mt < mt_percent)
} else {
  warning("percent.mt 列不在您的 metadata 中，将跳过线粒体过滤！")
}

if ("percent.hb" %in% names(object@meta.data)) {
  object1 <- subset(object1, subset = percent.hb < HB_percent)
} else {
  warning("percent.hb 列不在您的 metadata 中，将跳过红细胞过滤！")
}

if (rm_doub && "DF.classifications" %in% names(object@meta.data)) {
  object1 <- subset(object1, subset = DF.classifications == "Singlet")
}

if ("Contamination" %in% names(object@meta.data)) {
  object1 <- subset(object1, subset = Contamination < decontX)
}

cat("=> [3] 统计各个样本的过滤细胞数...\n")
samples <- unique(object@meta.data[[col_sample]])
del_cells_list <- list()

for (sample in samples) {
  sample_cells <- rownames(object@meta.data)[object@meta.data[[col_sample]] == sample]
  sample_meta <- object@meta.data[sample_cells, , drop=FALSE]
  
  sample_cells_filtered <- intersect(sample_cells, rownames(object1@meta.data))
  
  sample_data <- data.frame(
      sample = sample,
      total_cell = length(sample_cells),
      remaining_cell = length(sample_cells_filtered),
      low_nFeature = sum(sample_meta$nFeature_RNA <= A, na.rm = TRUE),
      high_nFeature = sum(sample_meta$nFeature_RNA >= B, na.rm = TRUE),
      low_nCount = sum(sample_meta$nCount_RNA <= CC, na.rm = TRUE),
      high_nCount = sum(sample_meta$nCount_RNA >= DD, na.rm = TRUE)
  )
  if ("percent.mt" %in% names(sample_meta)) {
    sample_data$high_percent_mt <- sum(sample_meta$percent.mt >= mt_percent, na.rm=TRUE)
  } else {
    sample_data$high_percent_mt <- 0
  }
  
  if ("percent.hb" %in% names(sample_meta)) {
    sample_data$high_HB <- sum(sample_meta$percent.hb >= HB_percent, na.rm=TRUE)
  } else {
    sample_data$high_HB <- 0
  }
  
  if (rm_doub && "DF.classifications" %in% names(sample_meta)) {
    sample_data$Doublet <- sum(sample_meta$DF.classifications == "Doublet", na.rm=TRUE)
  } else {
    sample_data$Doublet <- 0
  }
  
  if ("Contamination" %in% names(sample_meta)){
    sample_data$high_Contamination <- sum(sample_meta$Contamination >= decontX, na.rm=TRUE)
  } else {
    sample_data$high_Contamination <- 0
  }
  
  sample_data$all_filtered_cell <- length(sample_cells) - length(sample_cells_filtered)
  del_cells_list[[as.character(sample)]] <- sample_data
}

del_cells <- do.call(rbind, del_cells_list)
rownames(del_cells) <- NULL

write.csv(del_cells, file.path(outdir, "filter_stat_cell.csv"), quote=FALSE, row.names=FALSE)
# 直接在这里转置并覆盖原 python 脚本实现功能，消灭跨环境脚本依赖
stat_t <- t(del_cells)
write.table(stat_t, file.path(outdir, "filter_stat_cell.xls"), sep="\t", col.names=FALSE, quote=FALSE)

cat("=> [4] 绘制质控前后对比的小提琴图...\n")
plot_cols <- c("nFeature_RNA", "nCount_RNA")
if("percent.mt" %in% names(object@meta.data)) plot_cols <- c(plot_cols, "percent.mt")
if("percent.hb" %in% names(object@meta.data)) plot_cols <- c(plot_cols, "percent.hb")
if("Contamination" %in% names(object@meta.data)) plot_cols <- c(plot_cols, "Contamination")

# Plot Before
qc_b <- FeatureStatPlot(object, stat.by = plot_cols, group.by = col_sample, ncol=1, stack = TRUE, add_box=TRUE) & NoLegend() & theme(axis.title.y = element_blank())
wd <- max(7, ceiling(length(samples) * 0.2) + 7)
ggsave(file.path(outdir, "Vln_beforeQC.pdf"), plot=qc_b, width = wd, height = 7)
ggsave(file.path(outdir, "Vln_beforeQC.png"), plot=qc_b, width = wd, height = 7)

# Plot After
qc_a <- FeatureStatPlot(object1, stat.by = plot_cols, group.by = col_sample, ncol=1, stack = TRUE, add_box=TRUE) & NoLegend() & theme(axis.title.y = element_blank())
ggsave(file.path(outdir, "Vln_afterQC.pdf"), plot=qc_a, width = wd, height = 7)
ggsave(file.path(outdir, "Vln_afterQC.png"), plot=qc_a, width = wd, height = 7)

cat("=> [5] 保存最终过滤结果...\n")
saveRDS(object1, file.path(outdir, 'final_obj.rds'))
write.csv(object1@meta.data, file = file.path(outdir, "final_metadata.csv"), row.names = TRUE, quote = FALSE)

cat("=> 质控过滤分析完成！\n")
