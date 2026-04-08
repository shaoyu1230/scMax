#!/usr/bin/env Rscript

# ==========================================================
# 脚本名称：scSubCluster.R
# 脚本功能：亚群细分分析（重聚类、对小样本集特化的去批次与多测度解析）
# ==========================================================

suppressMessages(library(argparse))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(SCP))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))
suppressMessages(library(tidyr))

parser <- argparse::ArgumentParser(description = '单细胞数据亚群细分与专有去批次重聚类')

parser$add_argument('--rds', dest='inputrds', required=TRUE, help='输入的 Seurat RDS 文件，如上一层大群跑出来的 _integrate_clustered.rds')
parser$add_argument('--outdir', dest='outdir', required=TRUE, help='输出结果的总目录')
parser$add_argument('--subset_col', dest='subset_col', required=TRUE, help='必须提供取亚群的列名（比如 res0.4 或 CellType），这是决定提哪些细胞的关键')
parser$add_argument('--subset_val', dest='subset_val', required=TRUE, help='对应要提取的值（支持逗号分隔，比如提取大群中的 1,3,5 就填入 1,3,5）')
parser$add_argument('--methods', dest='methods', default="re_harmony,original_integrated", help='针对亚群特异的批次消除与评估策略: original_integrated, re_harmony, re_rpca, re_cca, re_none', type='character')
parser$add_argument('--resolutions', dest='resolutions', default="0.2,0.4", help='亚群再次聚类的分辨率尝试梯度', type='character')
parser$add_argument('--col_sample', dest='col_sample', default="Sample", help='样本列名称', type='character')
parser$add_argument('--col_group', dest='col_group', default="Group", help='分组列名称', type='character')
parser$add_argument('--integrate_by', dest='integrate_by', default="Sample", help='去批次整合时的依据列(默认Sample)', type='character')
parser$add_argument('--refmarker_file', dest='refmarker_file', default="", help='参考 Marker 基因集文件 (Tab 分隔的含有 CellType, Markers 列表格)', type='character')
parser$add_argument('--nfeatures', dest='nfeatures', default=1500, type='integer', help='亚群重新提取的高变基因数')
parser$add_argument('--npcs', dest='npcs', default=20, type='integer', help='亚群 PCA 降维参数（通常亚群不需要过大）')

parser$add_argument('--do_cluster', dest='do_cluster', action='store_true', help='执行高级群集特征表征(耗时长)')
parser$add_argument('--do_refmarker', dest='do_refmarker', action='store_true', help='执行参考靶点探查')
parser$add_argument('--orgdb', dest='orgdb', default='org.Hs.eg.db', type='character', help='富集分析的OrgDb')
parser$add_argument('--organism_kegg', dest='organism_kegg', default='hsa', type='character', help='富集分析的KEGG缩写')

opt <- parser$parse_args()

suppressMessages(library(this.path))
func_script <- file.path(dirname(this.path::this.dir()), "../05_celltype/script/func_scRNA_celltype_anno.R")
if (file.exists(func_script)) {
    source(func_script)
}

# ===== 辅助绘图函数 =====
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

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

cat("=> [1] 读取输入 RDS 数据并截取亚群...\n")
if (!file.exists(opt$inputrds)) stop(paste0("文件 ", opt$inputrds, " 不存在！"))
data <- readRDS(opt$inputrds)

cat(sprintf("=> 依据: %s 中提取亚群 [%s] ...\n", opt$subset_col, opt$subset_val))
vals <- trimws(unlist(strsplit(opt$subset_val, ",")))
if (opt$subset_col %in% colnames(data@meta.data)) {
  # 如果列是因子，转换为字符匹配
  cell_vals <- as.character(data@meta.data[[opt$subset_col]])
  cells_keep <- rownames(data@meta.data)[cell_vals %in% vals]
  if (length(cells_keep) > 0) {
    data <- subset(data, cells = cells_keep)
    cat(sprintf("=> 亚群提取完成，目标细胞数总计: %d\n", ncol(data)))
  } else {
    stop("提取后细胞数为 0，请检查列名与数值拼写或尝试换一个输入 RDS！")
  }
} else {
  stop(paste0("错误：列名 ", opt$subset_col, " 在提供的 RDS(metadata 中)不存在！没有数据用于子集划分。"))
}

methods <- trimws(unlist(strsplit(opt$methods, ",")))
resolutions <- as.numeric(trimws(unlist(strsplit(opt$resolutions, ","))))

# 统计分析批次情况，提供对如 CCA/RPCA 失败的预判与修正建议
batch_counts <- table(data@meta.data[[opt$integrate_by]])
cat("=> 各批次（整合基准列）拥有的样本数状态汇总:\n")
print(batch_counts)
min_batch_cell_num <- min(batch_counts)

for (method in methods) {
  cat(paste0("\n========================================\n"))
  cat(sprintf("=> [2] 开始执行亚群专属批次特征计算或整合，方法策略: %s\n", method))
  
  method_dir <- file.path(opt$outdir, method)
  dir.create(method_dir, showWarnings=FALSE, recursive=TRUE)
  
  tmp <- data
  is_success <- TRUE
  
  tryCatch({
    if (method == "original_integrated") {
      # 策略1：继承大群算出来的 Integrated 校正值而不重新找 Anchor，这对个别碎片细胞(数十个)的亚群最安全
      if ("integrated" %in% names(tmp@assays)) {
        cat("   - 检出继承的 integrated assay，安全剥离降维...\n")
        DefaultAssay(tmp) <- "integrated"
        tmp <- ScaleData(tmp, verbose = FALSE)
        tmp <- RunPCA(tmp, npcs = opt$npcs, verbose = FALSE)
        tmp <- RunTSNE(tmp, reduction = "pca", dims = 1:opt$npcs)
        tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:opt$npcs)
        tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:opt$npcs)
      } else {
        stop("没有在对象中发现 integrated 数据矩阵！这说明之前不是靠 RPCA/CCA 组装的结构，此策略在您数据中无法使用。")
      }
      
    } else if (method == "re_rpca" || method == "re_cca") {
      # 重新算，但是如果有的样本细胞数掉到极少值(比如<30)，CCA 和 RPCA 找靶向基因一定会失败。
      if(min_batch_cell_num < 30) {
         warning(sprintf("当前亚群存在细胞数等于 %d 的批次极小样本。RPCA/CCA 很有可能会由于难以产生稳定聚类描点而产生报错！", min_batch_cell_num))
      }
      DefaultAssay(tmp) <- "RNA"
      ifnb.list <- SplitObject(tmp, split.by = opt$integrate_by)
      ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = opt$nfeatures)
      })
      features <- SelectIntegrationFeatures(object.list = ifnb.list)
      
      if (method == "re_rpca") {
        ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
          x <- ScaleData(x, features = features, verbose = FALSE)
          x <- RunPCA(x, features = features, verbose = FALSE)
        })
        immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca", k.filter = max(10, min_batch_cell_num - 1))
      } else {
        immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, k.filter = max(10, min_batch_cell_num - 1))
      }
      
      tmp <- IntegrateData(anchorset = immune.anchors)
      DefaultAssay(tmp) <- "integrated"
      tmp <- ScaleData(tmp, verbose = FALSE)
      tmp <- RunPCA(tmp, npcs = opt$npcs, verbose = FALSE)
      tmp <- RunTSNE(tmp, reduction = "pca", dims = 1:opt$npcs)
      tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:opt$npcs)
      tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:opt$npcs)
      
    } else if (method == "re_harmony") {
      # Harmony 是最适合细分细胞的，不仅快而且能够应对样本间细胞数量不均的问题
      DefaultAssay(tmp) <- "RNA"
      suppressMessages(library(harmony))
      tmp <- NormalizeData(tmp) %>% FindVariableFeatures(nfeatures = opt$nfeatures) %>% ScaleData() %>% RunPCA(npcs = opt$npcs, verbose = FALSE)
      tmp <- RunHarmony(tmp, group.by.vars = opt$integrate_by)
      tmp <- RunTSNE(tmp, reduction = "harmony", dims = 1:opt$npcs)
      tmp <- RunUMAP(tmp, reduction = "harmony", dims = 1:opt$npcs)
      tmp <- FindNeighbors(tmp, reduction = "harmony", dims = 1:opt$npcs)
      
    } else if (method == "re_none") {
      DefaultAssay(tmp) <- "RNA"
      tmp <- NormalizeData(tmp) %>% FindVariableFeatures(nfeatures = opt$nfeatures) %>% ScaleData() %>% RunPCA(npcs = opt$npcs, verbose = FALSE)
      tmp <- RunTSNE(tmp, reduction = "pca", dims = 1:opt$npcs)
      tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:opt$npcs)
      tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:opt$npcs)
    } else {
      stop(paste0("未知的处理方法: ", method))
    }
  }, error = function(e){
    cat(sprintf("\n   !!! 警告: 策略 %s 在细分阶段因为细胞量或数值问题执行失败。报错信息: %s\n", method, conditionMessage(e)))
    cat("   建议：若属于 RPCA/CCA 报错，通常是某个批次的该亚群细胞数过少不足以建模锚点；建议切为 re_harmony(强大稳定) 或 original_integrated(保留原大图相对位置) \n")
    is_success <<- FALSE
  })
  
  if(!is_success) next # 失败则换下一个策略不跑了
  
  if(opt$col_sample %in% names(tmp@meta.data)) {
    p1 <- CellDimPlot(tmp, reduction='umap', group.by=opt$col_sample, pt.size=0.1, theme_use="theme_blank")
    ggsave(file.path(method_dir, paste0(method, "_Sample_UMAP.pdf")), p1, width=7, height=5)
  }
  if(opt$col_group %in% names(tmp@meta.data)) {
    p2 <- CellDimPlot(tmp, reduction='umap', group.by=opt$col_group, pt.size=0.1, theme_use="theme_blank")
    ggsave(file.path(method_dir, paste0(method, "_Group_UMAP.pdf")), p2, width=7, height=5)
  }
  
  cat("=> 导出针对子集整合基础 RDS ...\n")
  saveRDS(tmp, file.path(method_dir, paste0(method, "_subcluster_base.rds")))
  
  cat("=> [3] 开展亚群细分聚类与差异探针...\n")
  
  # 捕获刚刚去批次整合生成的特定 graph
  active_graph <- paste0(DefaultAssay(tmp), "_snn")
  
  DefaultAssay(tmp) <- "RNA"
  
  for (res in resolutions) {
    cat(sprintf("\n   - 亚群分辨解析 Resolution: %s\n", res))
    res_name <- paste0('res', res)
    res_dir <- file.path(method_dir, res_name)
    dir.create(res_dir, showWarnings=FALSE, recursive=TRUE)
    
    tmp <- FindClusters(tmp, resolution = res, graph.name = active_graph)
    tmp@meta.data$seurat_clusters <- factor(as.numeric(as.character(tmp@meta.data$seurat_clusters)))
    tmp[[res_name]] <- tmp$seurat_clusters
    
    pal <- palette_scp(unique(tmp$seurat_clusters), palette='Paired')
    
    pumap <- CellDimPlot(tmp, reduction='umap', group.by='seurat_clusters', label=TRUE, palcolor=pal, theme_use="theme_blank")
    ggsave(file.path(res_dir, "1_SubCluster_UMAP.pdf"), pumap, width=7, height=5)
    ptsne <- CellDimPlot(tmp, reduction='tsne', group.by='seurat_clusters', label=TRUE, palcolor=pal, theme_use="theme_blank")
    ggsave(file.path(res_dir, "1_SubCluster_TSNE.pdf"), ptsne, width=7, height=5)
    
    n_samp <- length(unique(tmp@meta.data[[opt$col_sample]]))
    w_split <- min(20, n_samp * 4 + 4)
    psplit <- CellDimPlot(tmp, reduction='umap', group.by='seurat_clusters', split.by=opt$col_sample, palcolor=pal, theme_use="theme_blank", ncol=min(4, n_samp))
    ggsave(file.path(res_dir, "1_SubCluster_split_Sample_UMAP.pdf"), psplit, width=w_split, height=min(15, ceiling(n_samp/4)*4))
    
    df_ratio <- as.data.frame(table(tmp@meta.data[[opt$col_sample]], tmp@meta.data$seurat_clusters))
    colnames(df_ratio) <- c("Sample", "SubCluster", "Freq")
    p_bar <- ggplot(df_ratio, aes(x=Sample, y=Freq, fill=SubCluster)) + 
             geom_bar(stat="identity", position="fill", width=0.7) + 
             theme_bw() + scale_fill_manual(values=pal) + 
             theme(axis.text.x=element_text(color="black", size=10, angle=45, hjust=1), legend.position="right") + 
             ylab("Fraction")
    ggsave(file.path(res_dir, "2_SubCluster_Sample_CellRatio.pdf"), p_bar, width=max(6, n_samp*0.4 + 2), height=5)
    write.csv(df_ratio, file.path(res_dir, "2_SubCluster_Sample_CellCount.csv"), quote=FALSE, row.names=FALSE)
    
    cat("     - 提取这支特殊亚群在分辨率下的精细 Marker ...\n")
    tryCatch({
      # 亚群由于细胞总数相对少，取消细胞取样上限保证挖掘精度
      markers <- FindAllMarkers(tmp, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.25)
      write.csv(markers, file.path(res_dir, "3_SubCluster_all_markers.csv"), row.names=FALSE)
      
      top10 <- markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
      p_dot <- DotPlot(tmp, features = unique(top10$gene), group.by="seurat_clusters") + 
               theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8), panel.border = element_rect(fill=NA, color="black"))
      
      dot_width <- max(10, length(unique(top10$gene)) * 0.2 + 2)
      dot_height <- max(5, length(unique(tmp$seurat_clusters)) * 0.2 + 2)
      ggsave(file.path(res_dir, "3_SubCluster_Top10_DotPlot.pdf"), p_dot, width=dot_width, height=dot_height)
    }, error = function(e) {
      cat("       (精细 Marker 基因计算由于某些类细胞严重不足导致报错跳过：", conditionMessage(e), ")\n")
    })
    
    # 绘制参考 Marker DotPlot (如果提供了标记列表)
    if (opt$refmarker_file != "") {
       plot_refmarkers_dot(tmp, opt$refmarker_file, res_dir)
    }
    
    # +++ 高级特征分析 (与 05_celltype 平齐，便于跳步软链) +++
    if (isTRUE(opt$do_cluster)) {
       cat("     - [高级选项] 正在执行全方位亚群特征表征 (cluster_characterization) ...\n")
       c_outdir <- file.path(res_dir, "cluster_characterization")
       dir.create(c_outdir, showWarnings=FALSE, recursive=TRUE)
       tryCatch({
         if (exists("cluster_plots")) {
             cluster_plots(tmp, c_outdir, groups=unique(c(opt$col_sample, opt$col_group)), palcolor=pal)
             cluster_de(tmp, rdsdir=file.path(c_outdir, "Rdata"), outdir=file.path(c_outdir,"4_Cluster.DE"))
             cluster_enrich(derds=file.path(c_outdir, "Rdata", "ClusterDEGs.rds"), orgdb=opt$orgdb, organism_kegg=opt$organism_kegg, enrichdir=file.path(c_outdir,"5_Enrich"))
         } else {
             warning("func_scRNA_celltype_anno.R 加载失败，跳过高级表征！")
         }
       }, error = function(e){
         cat("       (SubCluster 高级表征执行报错失败，将跳过当前组：", conditionMessage(e), ")\n")
       })
    }
    
    if (isTRUE(opt$do_refmarker) && opt$refmarker_file != "") {
       cat("     - [高级选项] 正在执行全套参考 Marker 映射 (marker_expression) ...\n")
       m_outdir <- file.path(res_dir, "marker_expression")
       dir.create(m_outdir, showWarnings=FALSE, recursive=TRUE)
       tryCatch({
         if (exists("plot_refmarkers")) {
            plot_refmarkers(tmp, opt$refmarker_file, outdir=m_outdir)
         } else {
            warning("func_scRNA_celltype_anno.R 加载失败，跳过全套 Reference Marker 探查！")
         }
       }, error = function(e){
         cat("       (高级参考 Marker 全套绘制失败：", conditionMessage(e), ")\n")
       })
    }
  }

  cat(sprintf("=> [%s] 综合导出子集细分最终对象文件 ...\n", method))
  saveRDS(tmp, file.path(method_dir, paste0(method, "_subcluster_final.rds")))
  
}

cat("\n========================================\n")
cat("=> 04_subcluster (高级亚群细分探针) 全部运行结束！结果存储于:", opt$outdir, "\n")
