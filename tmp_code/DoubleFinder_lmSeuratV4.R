#!/usr/bin/env Rscript
#/annogene/data2/share/software/SCV/SC_tools/Miniforge3/envs/lmSeuratV4/bin/Rscript

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(DoubletFinder)
  library(ggplot2)
  library(dplyr)
})

##############################
#     参数输入模块
##############################

option_list <- list(
  make_option(c("--outdir"), type="character", help="输出目录"),
  make_option(c("--input"), type="character", help="输入 rds 或 10X 目录"),
  make_option(c("--sample_col"), type="character", help="meta.data 中样本列名"),
  make_option(c("--samples"), type="character", default=NULL,
              help="可选：指定一个文件（每行一个样本名），只分析这些样本")
)

opt <- parse_args(OptionParser(option_list = option_list))

outdir      <- opt$outdir
input       <- opt$input
sample_col  <- opt$sample_col
sample_file <- opt$samples

if (is.null(outdir) | is.null(input) | is.null(sample_col)) {
  stop("❌ 参数不足：必须包含 --outdir --input --sample_col")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)


##############################
#      DoubletFinder 参数
##############################

pcs <- 1:20
pN <- 0.25
expectedDoubletRate <- 0.075
pK_search_pcs <- pcs
cluster_resolution <- 0.5


##############################
#        加载数据
##############################

sample_indir <- file.path(input, "outs/filtered_gene_bc_matrices/ref/")

if (file.exists(sample_indir)) {
  message("Detected 10X folder, reading with Read10X...")
  mat <- Read10X(sample_indir)
  seu <- CreateSeuratObject(mat)
} else {
  message("Reading Seurat rds object...")
  seu <- readRDS(input)
}

if (!(sample_col %in% colnames(seu@meta.data))) {
  stop(paste0("❌ meta.data 中找不到列：", sample_col))
}

# 创建 stim 列用于 subset
seu$stim <- seu@meta.data[[sample_col]]

message("全部检测到的样本：")
print(unique(seu$stim))


##############################
#     获取最终要分析的样本
##############################

if (!is.null(sample_file)) {
  if (!file.exists(sample_file)) stop("❌ 指定的 --samples 文件不存在")

  samples <- readLines(sample_file)
  samples <- samples[samples != ""]  # 去掉空行

  message("使用 --samples 指定的样本列表：")
  print(samples)

  # 过滤不存在的样本
  missing <- samples[!samples %in% seu$stim]
  if (length(missing) > 0) {
    message("⚠ 以下样本不在 Seurat 中被找到，将被忽略：")
    print(missing)
  }

  samples <- samples[samples %in% seu$stim]
  if (length(samples) == 0) stop("❌ 没有有效样本可供分析")
} else {
  samples <- unique(seu$stim)
  message("未提供 --samples，将分析全部样本：")
  print(samples)
}

seu$DF_result <- NA
seu$PANN_score <- NA


##############################
#        单个样本 DF 函数
##############################

run_DF_for_sample <- function(seu, sample_name) {

  message("\n==============================")
  message("开始处理样本：", sample_name)
  message("==============================")

  sub <- subset(seu, stim == sample_name)
  message("细胞数：", ncol(sub))

  sub <- NormalizeData(sub)
  sub <- FindVariableFeatures(sub)
  sub <- ScaleData(sub)
  sub <- RunPCA(sub)
  sub <- FindNeighbors(sub, dims = pcs)
  sub <- RunUMAP(sub, dims = 1:10)
  sub <- FindClusters(sub, resolution = cluster_resolution)

  #### pK sweep
  message("  pK sweep ...")
  sweep.res <- paramSweep(sub, PCs = pK_search_pcs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  pK_table <- find.pK(sweep.stats)
  best.pK <- as.numeric(as.character(
    pK_table$pK[which.max(pK_table$BCmetric)]
  ))
  message("最佳 pK = ", best.pK)

  #### homotypic correction
  nCells <- ncol(sub)
  nExp <- round(expectedDoubletRate * nCells)
  homotypic.prop <- modelHomotypic(sub$seurat_clusters)
  nExp.adj <- round(nExp * (1 - homotypic.prop))

  #### DoubletFinder
  message("  运行 DoubletFinder ...")
  sub <- doubletFinder(
    sub,
    PCs = pcs,
    pN = pN,
    pK = best.pK,
    nExp = nExp.adj,
    reuse.pANN = FALSE,
    sct = FALSE
  )

  df.col <- grep("^DF.classifications", colnames(sub@meta.data), value = TRUE)
  pANN.col <- grep("^pANN", colnames(sub@meta.data), value = TRUE)

  cells <- colnames(sub)
  seu@meta.data[cells, "DF_result"] <- sub[[df.col]]
  seu@meta.data[cells, "PANN_score"] <- sub[[pANN.col]]

  return(seu)
}

##############################
#        循环样本分析
##############################

for (s in samples) {
  message(paste0("开始分析：", s))
  seu <- run_DF_for_sample(seu, s)
}


##############################
#        输出结果
##############################

write.table(seu@meta.data, "DoubletFinder_meta.data.xls",
            sep = "\t", quote = FALSE, row.names = FALSE)

# UMAP 结果
 P <- DimPlot(seu, group.by = "DF_result") + ggtitle("DoubletFinder Predictions")
 ggsave("DoubletFinder_Predictions_umap.pdf", P, width=6, height=5)

 pdf("DF_Score_FeaturePlot.pdf", width = 6, height = 5)
 p <- FeaturePlot(seu, features = "PANN_score", min.cutoff = "q9",
                # cols = c("lightgrey", "red"), order = TRUE)
 print(p)
 dev.off()

# summary statistics
df_summary <- seu@meta.data %>%
  group_by(!!as.symbol(sample_col), DF_result) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(!!as.symbol(sample_col)) %>%
  mutate(prop = count / sum(count))

write.csv(df_summary, "DoubletFinder_summary_by_sample.csv",
          row.names = FALSE, quote = FALSE)

message("🎉 DoubletFinder 完成！")
