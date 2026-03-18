#!/usr/bin/env Rscript

# ==========================================================
# 脚本名称：scCluster.R
# 脚本功能：单细胞数据去批次与多分辨率降维聚类分析
# ==========================================================

suppressMessages(library(argparse))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(SCP))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))
suppressMessages(library(tidyr))

parser <- argparse::ArgumentParser(description = '单细胞数据多策略整合去批次及多分辨率聚类')

parser$add_argument('--rds', dest='inputrds', required=TRUE, help='输入的 Seurat RDS 文件，如 02_filter 产出的 final_obj.rds')
parser$add_argument('--outdir', dest='outdir', required=TRUE, help='输出结果的总目录')
parser$add_argument('--subset_col', dest='subset_col', default="", help='取子集的列名（可选）', type='character')
parser$add_argument('--subset_val', dest='subset_val', default="", help='取子集的值，支持逗号分隔多个值（可选）', type='character')
parser$add_argument('--methods', dest='methods', default="harmony", help='去批次整合方法，支持逗号分隔执行多个作对比: rpca,harmony,cca,none', type='character')
parser$add_argument('--resolutions', dest='resolutions', default="0.2,0.4,0.6", help='聚类分辨率，支持逗号分隔多个', type='character')
parser$add_argument('--refmarker_file', dest='refmarker_file', default="", help='参考 Marker 基因集文件 (Tab 分隔的含有 CellType, Markers 列表格)', type='character')
parser$add_argument('--col_sample', dest='col_sample', default="Sample", help='样本列名称', type='character')
parser$add_argument('--col_group', dest='col_group', default="Group", help='分组列名称', type='character')
parser$add_argument('--integrate_by', dest='integrate_by', default="Sample", help='去批次整合时的依据列(默认Sample)', type='character')
parser$add_argument('--nfeatures', dest='nfeatures', default=2000, type='integer', help='高变基因数')
parser$add_argument('--npcs', dest='npcs', default=30, type='integer', help='PCA 主成分选取数量')

opt <- parser$parse_args()

# ===== 辅助函数 =====
plot_refmarkers_dot <- function(data, refmarker.file, outdir) {
  if (!file.exists(refmarker.file)) {
    warning(paste0("参考 marker 文件 ", refmarker.file, " 不存在，跳过 DotPlot 绘制。"))
    return()
  }
  mks_ref <- read.csv(refmarker.file, sep='\t', stringsAsFactors=FALSE)
  if(!("CellType" %in% colnames(mks_ref)) || !("Markers" %in% colnames(mks_ref))){
    warning("参考 marker 文件必须包含 CellType 和 Markers 列，跳过 DotPlot 绘制。")
    return()
  }
  
  df_long <- mks_ref[c('CellType','Markers')] %>% tidyr::separate_rows(Markers, sep = "[, ]+")
  df_unique <- df_long %>% unique() %>% as.data.frame()
  df_unique$Markers <- gsub(' ','', df_unique$Markers)
  valid_mks <- intersect(df_unique$Markers, rownames(data))
  if (length(valid_mks) == 0) {
    warning("提供的参考 marker 基因在数据中均未找到。")
    return()
  }
  df_unique <- df_unique[df_unique$Markers %in% valid_mks, ]
  
  cat("     - 正在绘制 SCP GroupHeatmap DotPlot ...\n")
  tryCatch({
    res_scp <- GroupHeatmap(data, features = df_unique$Markers,
                       feature_split = df_unique$CellType, assay = "RNA", slot = 'data',
                       group.by = "seurat_clusters", heatmap_palette = "RdBu",
                       show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 0, cluster_columns = TRUE,
                       anno_terms = FALSE, nlabel=0, add_dot = TRUE, add_bg = TRUE)
    
    ncluster <- length(unique(data$seurat_clusters))
    wd1 <- max(8, ncluster * 0.6)
    nmarkers <- nrow(df_unique)
    ht1 <- min(max(8, nmarkers * 0.5), 30)
    ggsave(file.path(outdir, "4_RefMarker_GroupHeatmap_DotPlot.pdf"), plot=res_scp$plot, width=wd1, height=ht1)
    ggsave(file.path(outdir, "4_RefMarker_GroupHeatmap_DotPlot.png"), plot=res_scp$plot, width=wd1, height=ht1)
  }, error = function(e){
    cat("       (SCP::GroupHeatmap 绘制失败: ", conditionMessage(e), ")\n")
  })
  
  cat("     - 正在绘制 Seurat DotPlot ...\n")
  tryCatch({
    p <- DotPlot(data, features = unique(df_unique$Markers), group.by="seurat_clusters") + 
      scale_colour_gradientn(colours = c("dodgerblue1","#44c1f0", "lightgoldenrod","#e20612",'#cc340c')) +
      geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2) +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), axis.title=element_blank())
      
    mks <- unique(df_unique$Markers)
    wd2 <- max(ceiling(length(mks)/3), 9)
    ht2 <- max(length(unique(data$seurat_clusters)) * 0.4, 5)
    ggsave(file.path(outdir, "4_RefMarker_Seurat_DotPlot.pdf"), plot=p, width=wd2, height=ht2)
    ggsave(file.path(outdir, "4_RefMarker_Seurat_DotPlot.png"), plot=p, width=wd2, height=ht2)
  }, error = function(e){
    cat("       (Seurat::DotPlot 绘制失败: ", conditionMessage(e), ")\n")
  })
}
# =====================

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

cat("=> [1] 读取输入 RDS 数据...\n")
if (!file.exists(opt$inputrds)) stop(paste0("文件 ", opt$inputrds, " 不存在！"))
data <- readRDS(opt$inputrds)

if (opt$subset_col != "" && opt$subset_val != "") {
  cat(sprintf("=> 根据参数取特征子集: %s 在 [%s] 中...\n", opt$subset_col, opt$subset_val))
  vals <- trimws(unlist(strsplit(opt$subset_val, ",")))
  if (opt$subset_col %in% colnames(data@meta.data)) {
    cells_keep <- rownames(data@meta.data)[data@meta.data[[opt$subset_col]] %in% vals]
    if (length(cells_keep) > 0) {
      data <- subset(data, cells = cells_keep)
      cat(sprintf("=> 取子集完成，剩余处理细胞数: %d\n", ncol(data)))
    } else {
      stop("取子集后细胞数为 0，请检查分类列名或数值是否拼写正确！")
    }
  } else {
    warning("输入的 subset_col 在 metadata 中不存在，分析已自动跳过取子集操作。")
  }
}

methods <- trimws(unlist(strsplit(opt$methods, ",")))
resolutions <- as.numeric(trimws(unlist(strsplit(opt$resolutions, ","))))

for (method in methods) {
  cat(paste0("\n========================================\n"))
  cat(sprintf("=> [2] 开始执行去批次整合，方法: %s\n", method))
  
  method_dir <- file.path(opt$outdir, method)
  dir.create(method_dir, showWarnings=FALSE, recursive=TRUE)
  
  DefaultAssay(data) <- "RNA"
  tmp <- data
  
  if (method == "rpca") {
    ifnb.list <- SplitObject(tmp, split.by = opt$integrate_by)
    ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = opt$nfeatures)
    })
    features <- SelectIntegrationFeatures(object.list = ifnb.list)
    ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
    })
    immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
    tmp <- IntegrateData(anchorset = immune.anchors)
    DefaultAssay(tmp) <- "integrated"
    tmp <- ScaleData(tmp, verbose = FALSE)
    tmp <- RunPCA(tmp, npcs = opt$npcs, verbose = FALSE)
    tmp <- RunTSNE(tmp, reduction = "pca", dims = 1:opt$npcs)
    tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:opt$npcs)
    tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:opt$npcs)
    
  } else if (method == "harmony") {
    suppressMessages(library(harmony))
    tmp <- NormalizeData(tmp) %>% FindVariableFeatures(nfeatures = opt$nfeatures) %>% ScaleData() %>% RunPCA(npcs = opt$npcs, verbose = FALSE)
    tmp <- RunHarmony(tmp, group.by.vars = opt$integrate_by)
    tmp <- RunTSNE(tmp, reduction = "harmony", dims = 1:opt$npcs)
    tmp <- RunUMAP(tmp, reduction = "harmony", dims = 1:opt$npcs)
    tmp <- FindNeighbors(tmp, reduction = "harmony", dims = 1:opt$npcs)
    
  } else if (method == "cca") {
    ifnb.list <- SplitObject(tmp, split.by = opt$integrate_by)
    ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = opt$nfeatures)
    })
    features <- SelectIntegrationFeatures(object.list = ifnb.list)
    immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
    tmp <- IntegrateData(anchorset = immune.anchors)
    DefaultAssay(tmp) <- "integrated"
    tmp <- ScaleData(tmp, verbose = FALSE)
    tmp <- RunPCA(tmp, npcs = opt$npcs, verbose = FALSE)
    tmp <- RunTSNE(tmp, reduction = "pca", dims = 1:opt$npcs)
    tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:opt$npcs)
    tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:opt$npcs)
    
  } else if (method == "none") {
    tmp <- NormalizeData(tmp) %>% FindVariableFeatures(nfeatures = opt$nfeatures) %>% ScaleData() %>% RunPCA(npcs = opt$npcs, verbose = FALSE)
    tmp <- RunTSNE(tmp, reduction = "pca", dims = 1:opt$npcs)
    tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:opt$npcs)
    tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:opt$npcs)
  } else {
    warning(paste0("未知方法 ", method, "，已跳过..."))
    next
  }
  
  # 基础降维聚类全貌图
  if(opt$col_sample %in% names(tmp@meta.data)) {
    p1 <- CellDimPlot(tmp, reduction='umap', group.by=opt$col_sample, pt.size=0.01, theme_use="theme_blank")
    ggsave(file.path(method_dir, paste0(method, "_Sample_UMAP.pdf")), p1, width=7, height=5)
  }
  if(opt$col_group %in% names(tmp@meta.data)) {
    p2 <- CellDimPlot(tmp, reduction='umap', group.by=opt$col_group, pt.size=0.01, theme_use="theme_blank")
    ggsave(file.path(method_dir, paste0(method, "_Group_UMAP.pdf")), p2, width=7, height=5)
  }
  
  cat("=> 导出基础整合后的基础 RDS ...\n")
  saveRDS(tmp, file.path(method_dir, paste0(method, "_integrate_base.rds")))
  
  cat("=> [3] 开展多分辨率尝试与细胞聚类...\n")
  DefaultAssay(tmp) <- "RNA"
  
  for (res in resolutions) {
    cat(sprintf("\n   - 正在处理分群 Resolution: %s\n", res))
    res_name <- paste0('res', res)
    res_dir <- file.path(method_dir, res_name)
    dir.create(res_dir, showWarnings=FALSE, recursive=TRUE)
    
    tmp <- FindClusters(tmp, resolution = res)
    # 重编号，避免混乱的因子级
    tmp@meta.data$seurat_clusters <- factor(as.numeric(as.character(tmp@meta.data$seurat_clusters)))
    tmp[[res_name]] <- tmp$seurat_clusters
    
    pal <- palette_scp(unique(tmp$seurat_clusters), palette='Paired')
    
    pumap <- CellDimPlot(tmp, reduction='umap', group.by='seurat_clusters', label=TRUE, palcolor=pal, theme_use="theme_blank")
    ggsave(file.path(res_dir, "1_Cluster_UMAP.pdf"), pumap, width=7, height=5)
    ptsne <- CellDimPlot(tmp, reduction='tsne', group.by='seurat_clusters', label=TRUE, palcolor=pal, theme_use="theme_blank")
    ggsave(file.path(res_dir, "1_Cluster_TSNE.pdf"), ptsne, width=7, height=5)
    
    # 按样本拆分降维绘图
    n_samp <- length(unique(tmp@meta.data[[opt$col_sample]]))
    w_split <- min(20, n_samp * 4 + 4)
    psplit <- CellDimPlot(tmp, reduction='umap', group.by='seurat_clusters', split.by=opt$col_sample, palcolor=pal, theme_use="theme_blank", ncol=min(4, n_samp))
    ggsave(file.path(res_dir, "1_Cluster_split_Sample_UMAP.pdf"), psplit, width=w_split, height=min(15, ceiling(n_samp/4)*4))
    
    # 比例堆叠条形图
    df_ratio <- as.data.frame(table(tmp@meta.data[[opt$col_sample]], tmp@meta.data$seurat_clusters))
    colnames(df_ratio) <- c("Sample", "Cluster", "Freq")
    p_bar <- ggplot(df_ratio, aes(x=Sample, y=Freq, fill=Cluster)) + 
             geom_bar(stat="identity", position="fill", width=0.7) + 
             theme_bw() + scale_fill_manual(values=pal) + 
             theme(axis.text.x=element_text(color="black", size=10, angle=45, hjust=1), legend.position="right") + 
             ylab("Fraction")
    ggsave(file.path(res_dir, "2_Cluster_Sample_CellRatio.pdf"), p_bar, width=max(6, n_samp*0.4 + 2), height=5)
    write.csv(df_ratio, file.path(res_dir, "2_Cluster_Sample_CellCount.csv"), quote=FALSE, row.names=FALSE)
    
    # 找 Marker 基因
    cat("     - 正在计算基于该分辨率下的亚群核心 Marker 基因 (为加速限定 max 300 cells) ...\n")
    tryCatch({
      # 这里不再过度依赖特殊包，使用公版的 FindAllMarkers 快速验证
      markers <- FindAllMarkers(tmp, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.25, max.cells.per.ident = 300)
      write.csv(markers, file.path(res_dir, "3_Cluster_all_markers.csv"), row.names=FALSE)
      
      top10 <- markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
      p_dot <- DotPlot(tmp, features = unique(top10$gene), group.by="seurat_clusters") + 
               theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8), panel.border = element_rect(fill=NA, color="black"))
      
      dot_width <- max(10, length(unique(top10$gene)) * 0.2 + 2)
      dot_height <- max(5, length(unique(tmp$seurat_clusters)) * 0.2 + 2)
      ggsave(file.path(res_dir, "3_Cluster_Top10_DotPlot.pdf"), p_dot, width=dot_width, height=dot_height)
    }, error = function(e) {
      cat("       (Marker 基因计算跳过或失败：", conditionMessage(e), ")\n")
    })
    
    # 绘制参考 Marker DotPlot (如果提供了标记列表)
    if (opt$refmarker_file != "") {
       plot_refmarkers_dot(tmp, opt$refmarker_file, res_dir)
    }
  }

  cat(sprintf("=> [%s] 综合导出最终对象文件 ...\n", method))
  saveRDS(tmp, file.path(method_dir, paste0(method, "_integrate_clustered.rds")))
  
}

cat("\n========================================\n")
cat("=> 03_cluster (去批次整合与多分辨率计算) 全部运行结束！结果存储于:", opt$outdir, "\n")
