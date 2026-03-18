#!/usr/bin/env Rscript

# ==========================================================
# 脚本名称：scQC.R
# 脚本功能：单细胞/单核RNA质量控制、批次效应初步展示、双细胞和污染预估
# ==========================================================

suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(R.utils))
suppressMessages(library(Seurat))
suppressMessages(library(DoubletFinder))
suppressMessages(library(celda))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(ggpubr))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(SCP))
suppressMessages(library(ggridges))
suppressMessages(library(yaml))
suppressMessages(library(jsonlite))

bindir <- this.path::this.dir()
source(file.path(bindir,'func_scQC.R'))

# -----------------
# 参数定义 (暴露给外部调用)
# -----------------
parser <- argparse::ArgumentParser(description = '单细胞转录组 QC 和初步批次处理脚本')

parser$add_argument('--rds', dest='inputrds', default="", help='提供合并好或预先存好的 Seurat RDS 格式输入文件', type='character')
parser$add_argument('--samples', dest='samples_file', default="", help='输入样本的 CSV 文件，包含 SampleName 和 FilePath 列 (与 rds 互斥)', type='character')
parser$add_argument('--method', dest='method', default="merge", help='多样本整合方法(merge/cca/harmony)', type='character')
parser$add_argument('--meta', dest='metafile', default="", help='可选的新 metadata 补充文件 (CSV)', type='character')

parser$add_argument('--col_sample', dest='col_sample', default="Sample", help='样本所在原属性列')
parser$add_argument('--col_group', dest='col_group', default="Group", help='分组所在原属性列')
parser$add_argument('--species', dest='species', default="mouse", help='物种类别（human/mouse/other）')
parser$add_argument('--doub_rate', dest='doub_rate', default=0.076, type="double", help='预计双细胞比率')
parser$add_argument('--mt_genes', dest='mt_genes', default="", type="character", help='指定的线粒体基因逗号分隔表，留空则自动')
parser$add_argument('--hb_genes', dest='hb_genes', default="", type="character", help='指定的血红蛋白基因逗号分隔表，留空则自动')
parser$add_argument('--run_doublet', dest='run_doublet', action='store_true', default=FALSE, help='是否运行 DoubletFinder 双细胞预测')
parser$add_argument('--run_decontx', dest='run_decontx', action='store_true', default=FALSE, help='是否运行 DecontX 游离 RNA 污染估算')

parser$add_argument('-c', '--config', dest='config', default="", help='(兼容旧版) 传统 JSON/YAML 文件读取。如果提供则覆盖上面的参数', type='character')
parser$add_argument('-o', '--outdir', dest='outdir', required=TRUE, help='输出结果存储目录')

opt <- parser$parse_args()

# 如果有配置文件传入，则兼容原本的流程设置，覆盖通过命令行指定的变体（保持向下兼容）
if (opt$config != "") {
  suppressMessages(library(configr))
  config <- read.config(opt$config) 
  inputrds <- config$data$f01_rdsfile
  metafile <- config$data$f02_metadatafile
  col_sample <- config$param$param_ana$col_sample
  col_group <- config$param$param_ana$col_group
  species <- config$param$param_ana$species
  mt_genes <- config$param$param_ana$mt_genes
  hb_genes <- config$param$param_ana$hb_genes
  doub_rate <- config$param$param_ana$doub_rate
  run_doublet <- if(!is.null(config$param$param_ana$run_doublet)) config$param$param_ana$run_doublet else TRUE
  run_decontx <- if(!is.null(config$param$param_ana$run_decontx)) config$param$param_ana$run_decontx else TRUE
} else {
  inputrds <- opt$inputrds
  samples_file <- opt$samples_file
  method <- opt$method
  metafile <- opt$metafile
  col_sample <- opt$col_sample
  col_group <- opt$col_group
  species <- opt$species
  mt_genes <- opt$mt_genes
  hb_genes <- opt$hb_genes
  doub_rate <- opt$doub_rate
  run_doublet <- opt$run_doublet
  run_decontx <- opt$run_decontx
}

outdir <- opt$outdir

##### 创建相关的目录架构 ####
cat("=> [1] 初始化输出文件目录...\n")
dir.create(file.path(outdir, "tmp"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "1_Batch"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "2_QC"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "3_Doublet"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "4_Contamination"), showWarnings = FALSE, recursive = TRUE)

##### 读取数据或整合合并数据并核准 Metadata #####
cat("=> [2] 获取输入数据...\n")
if (samples_file != "") {
  if (!file.exists(samples_file)) {
    stop(paste0("错误：找不到样本信息文件 ", samples_file))
  }
  samples_info <- read.csv(samples_file, stringsAsFactors = FALSE)
  cat(sprintf("=> 共找到 %d 个样本信息，开始读取数据并整合...\n", nrow(samples_info)))
  
  seurat_list <- list()
  for (i in 1:nrow(samples_info)) {
    sample_id <- samples_info$SampleName[i]
    file_path <- samples_info$FilePath[i]
    message(paste0("读取样本: ", sample_id))
    
    if (grepl("\\.rds$", file_path, ignore.case = TRUE)) {
      obj <- readRDS(file_path)
    } else if (dir.exists(file_path)) {
      counts <- Read10X(data.dir = file_path)
      obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 200)
    } else {
      warning("未知的文件格式，跳过该样本。")
      next
    }
    
    obj$orig.ident <- sample_id
    obj$Sample <- sample_id
    
    for (col in colnames(samples_info)) {
      if (!(col %in% c("SampleName", "FilePath"))) {
        obj[[col]] <- samples_info[i, col]
      }
    }
    seurat_list[[sample_id]] <- obj
  }
  
  cat(sprintf("=> 开始使用 %s 策略进行多样本数据整合...\n", method))
  if (method == "merge") {
    data <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "Data")
  } else {
    stop("目前仅支持基于 merge 的直接合并，后续将扩展。")
  }
} else if (inputrds != "") {
  if (!file.exists(inputrds)) {
    stop(paste0("输入RDS文件 ", inputrds, " 不存在!"))
  }
  data <- readRDS(inputrds)
} else {
  stop("错误：必须提供 --rds 或者 --samples 作为数据输入源！")
}
if(is.null(metafile) || metafile == ""){
  warning("提示：本次分析未提供额外的 metadata 覆盖，将使用 RDS 内部原始元数据信息\n")
} else {
  new_meta <- read.csv(metafile, row.names = 1, stringsAsFactors = FALSE)
  if (!identical(sort(rownames(new_meta)), sort(rownames(data@meta.data)))){
    missing_in_new <- setdiff(rownames(data@meta.data), rownames(new_meta))
    missing_in_old <- setdiff(rownames(new_meta), rownames(data@meta.data))
    cat("错误：细胞 ID 不匹配！\n")
    cat("新 metadata 中缺失的细胞前几位:", paste(head(missing_in_new), collapse = ","), "...\n")
    cat("原数据中缺失的细胞前几位:", paste(head(missing_in_old), collapse = ","), "...\n")
    stop("提供的 Metadata 数据与对象细胞无法对齐，停止分析程序。")
  }else{
    cat("=> 元数据细胞 ID 完全匹配，已进行安全替换。\n")
  }
  # 使用与 data 对应的正确顺序写入 new_meta 
  new_meta <- new_meta[rownames(data@meta.data), , drop = FALSE]
  data@meta.data <- new_meta
}

#### 基础 QC 和指标计算 ####
cat("=> [3] 基础质量控制(QC)与指标标注...\n")
# 将逗号分隔字符串转为序列
format_genes <- function(gstr) {
  if (is.character(gstr) && trimws(gstr) != "") {
    vec <- strsplit(gstr, ",")[[1]]
    vec <- trimws(vec)
    return(vec[vec != ""])
  }
  return(NULL)
}

mt_genes_list <- format_genes(mt_genes)
hb_genes_list <- format_genes(hb_genes)

mt_genes <- get_mtgenes(data, species = species, genelist = mt_genes_list)
hb_genes <- get_hbgenes(data, species = species, genelist = hb_genes_list)

if(length(mt_genes) < 1){
  data[["percent.mt"]] <- 0
} else {
  data[["percent.mt"]] <- PercentageFeatureSet(data, features = mt_genes)
}

if(length(hb_genes) < 1){
  data[["percent.hb"]] <- 0
} else {
  data[["percent.hb"]] <- PercentageFeatureSet(data, features = hb_genes)
}

#### 1_Batch (基础整合和展示) ####
cat("=> [4] 运行批次降维可视化...\n")
data <- qc_cluster(data, col_sample = col_sample, col_group = col_group, nPCs = 30, runtsne = TRUE, outdir = file.path(outdir, "1_Batch"))

#### 2_QC Plots ###
cat("=> [5] 绘制基础质控(QC)小提琴图和散点映射...\n")
# plot
qc_1 <- FeatureStatPlot(data, stat.by = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), group.by = col_sample, ncol = 1, stack = TRUE) & NoLegend() & theme(axis.title.y = element_blank())
nsample <- length(unique(data@meta.data[, col_sample]))
wd <- max(6, nsample * 0.6)
ggsave(file.path(outdir, "2_QC/1_QC-vln.pdf"), plot = qc_1, width = wd, height = 6)
ggsave(file.path(outdir, "2_QC/1_QC-vln.png"), plot = qc_1, width = wd, height = 6)

FeatureDimPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 2, reduction = 'UMAP', theme_use = 'theme_blank', pt.size = 0.1) 
ggsave(file.path(outdir, "2_QC/2_QC-UMAP.pdf"), width = 10, height = 8)
ggsave(file.path(outdir, "2_QC/2_QC-UMAP.png"), width = 10, height = 8)

FeatureDimPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 2, reduction = 'TSNE', theme_use = 'theme_blank', pt.size = 0.1)
ggsave(file.path(outdir, "2_QC/2_QC-TSNE.pdf"), width = 10, height = 8)
ggsave(file.path(outdir, "2_QC/2_QC-TSNE.png"), width = 10, height = 8)

#### 3_Doublet (DoubletFinder 推断双细胞) ####
if (run_doublet) {
  cat("=> [6] 运行 DoubletFinder 双细胞预测模块...\n")
  run_doublet_prediction <- function(data) {
    samples <- unique(data$orig.ident)
    
    for(samp in samples) {
      # run_DF_for_sample 会内部自带 subset 机制并在源大对象上修改 meta
      data <- run_DF_for_sample(data, samp, expectedDoubletRate = doub_rate, sample_col = "orig.ident")
    }
    
    # 兼容并统一原来的命名列与绘图
    data$doublets_score <- data$PANN_score
    data$DF.classifications <- data$DF_result
    
    # 抽取以保存到 csv
    db_fd <- data@meta.data[, c('doublets_score', 'DF.classifications')]
    FeatureDimPlot(data,features='doublets_score', reduction = 'UMAP', pt.size = 0.1)
    ggsave(file.path(outdir, "3_Doublet/1_DoubletScore-UMAP.pdf"), width = 7, height = 5)
    ggsave(file.path(outdir, "3_Doublet/1_DoubletScore-UMAP.png"), width = 7, height = 5)
    
    CellDimPlot(data, group.by = 'DF.classifications', pt.size = 0.01, reduction = 'UMAP', palcolor = c('#b52b2b', '#4682b4')) 
    ggsave(file.path(outdir, "3_Doublet/2_DoubletFinder-UMAP.pdf"), width = 7, height = 5)
    ggsave(file.path(outdir, "3_Doublet/2_DoubletFinder-UMAP.png"), width = 7, height = 5)
    
    CellDimPlot(data, group.by = 'DF.classifications', pt.size = 0.01, reduction = 'TSNE', palcolor = c('#b52b2b', '#4682b4'))
    ggsave(file.path(outdir, "3_Doublet/3_DoubletFinder-TSNE.pdf"), width = 7, height = 5)
    ggsave(file.path(outdir, "3_Doublet/3_DoubletFinder-TSNE.png"), width = 7, height = 5)
    
    write.csv(db_fd, file = file.path(outdir, '3_Doublet/DoubletFinder_res.csv'), row.names = TRUE, quote = FALSE)
    
    # 统计各样本的双细胞情况
    if("DF.classifications" %in% names(data@meta.data)) {
        stat <- as.data.frame(table(data$orig.ident, data$DF.classifications))
        colnames(stat) <- c("Sample", "DF.classifications", "cell_counts")
        write.csv(stat, file.path(outdir, '3_Doublet/4_Doublet_stats.csv'), quote = FALSE, row.names = FALSE)
    }
    return(data)
  }
  
  data <- tryCatch({
    run_doublet_prediction(data)
  }, error = function(e){
    cat(sprintf("=> DoubletFinder过程报错：%s\n", e))
    return(data)
  })
} else {
  cat("=> [6] 跳过 DoubletFinder 双细胞预测模块...\n")
}

#### 4_Contamination (DecontX) ####
if (run_decontx) {
  cat("=> [7] 使用 DecontX 对游离 RNA 污染进行估算...\n")
  run_decontx_prediction <- function(data) {
    decontx_df <- data.frame()
    samples <- unique(data$orig.ident)
    for(samp in samples) {
      cat(sprintf("=> Processing DecontX for %s...\n", samp))
      object <- subset(data, orig.ident == samp)
      res <- get_decontx(object, prefix = file.path(outdir, 'tmp', samp))
      decontx_df <- rbind(decontx_df, res)
    }
    write.csv(decontx_df, file = file.path(outdir, '4_Contamination/DecontX_res.csv'), row.names = TRUE, quote = FALSE)
    data <- AddMetaData(data, decontx_df)
    
    pp <- ggplot(decontx_df, aes(x = Contamination, y = orig.ident, fill = orig.ident)) +
      geom_density_ridges_gradient() +
      labs(title = "Density Distribution of Contamination", x = "Contamination") +
      theme_ridges(font_size = 13, grid = TRUE) & NoLegend()
      
    ht <- 5 + ceiling((length(unique(decontx_df$orig.ident)) - 6) * 0.2)
    ggsave(file.path(outdir, '4_Contamination/1_Contamination-Density.pdf'), plot = pp, width = 5, height = ht)
    ggsave(file.path(outdir, '4_Contamination/1_Contamination-Density.png'), plot = pp, width = 5, height = ht)
    
    FeatureDimPlot(data, features = 'Contamination', pt.size = 0.01, theme_use = 'theme_blank', reduction = 'UMAP')
    ggsave(file.path(outdir, "4_Contamination/2_DecontX-UMAP.pdf"), width = 6, height = 5)
    ggsave(file.path(outdir, "4_Contamination/2_DecontX-UMAP.png"), width = 6, height = 5)
    
    FeatureDimPlot(data, features = 'Contamination', pt.size = 0.01, theme_use = 'theme_blank', reduction = 'TSNE')
    ggsave(file.path(outdir, "4_Contamination/2_DecontX-TSNE.pdf"), width = 6, height = 5)
    ggsave(file.path(outdir, "4_Contamination/2_DecontX-TSNE.png"), width = 6, height = 5)
    
    return(data)
  }
  
  data <- tryCatch({
    run_decontx_prediction(data)
  }, error = function(e){
    cat(sprintf("=> DecontX过程报错：%s\n", e))
    return(data)
  })
} else {
  cat("=> [7] 跳过 DecontX 游离 RNA 污染估算...\n")
}

# ================================
# 清理和导出部分
# ================================
cat("=> [8] 清理过载元胞属性分群遗留缓存并封存数据...\n")
rm.cols <- c(grep('RNA_snn', colnames(data@meta.data), value=TRUE), 'seurat_clusters')
keep.cols <- setdiff(colnames(data@meta.data), rm.cols)
data@meta.data <- data@meta.data[, keep.cols]

# 封存对象至本步质控最终输出点
cat("=> 导出 final_obj.rds 和 metadata \n")
saveRDS(data, file = file.path(outdir, "final_obj.rds"))
write.csv(data@meta.data, file = file.path(outdir, "final_metadata.csv"), row.names = TRUE, quote = FALSE)

cat("=> 质控分析流程结束，相关数据输出在当前 outdir 内。\n")
