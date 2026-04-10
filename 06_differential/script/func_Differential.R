# ==========================================================
# 库函数名称：func_Differential.R
# 功能：差异分析核心函数库 (全标签差异 vs 精细组间两两比较)
# ==========================================================

suppressMessages(library(Seurat))
suppressMessages(library(openxlsx))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(SCP))

# --- [模式 A] 全标签组差异核心函数 ---
# input: Seurat对象, 比较列, 输出路径
cluster_de_all <- function(data, assay='RNA', groupby='seurat_clusters', outdir){
  cat("=> [DE] 启动模式 A: 全标签组差异挖掘 (Global Markers)...\n")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # 设置标识
  if (!(groupby %in% colnames(data@meta.data))) {
    stop(sprintf("错误: 比较列 %s 不在 Metadata 中!", groupby))
  }
  
  Idents(data) <- data@meta.data[[groupby]]
  DefaultAssay(data) <- assay
  
  cat("   - 寻找差异基因 (FindAllMarkers)... \n")
  # 优先尝试执行快速差异分析
  degs <- RunPrestoAll(data, group.by = groupby, only.pos = FALSE)
  
  # 排序与判定
  degs <- as.data.frame(degs)
  degs <- setorderv(degs, c("cluster", "avg_log2FC", "pct.1", "pct.2"), c(1, -1, 1, -1))
  degs$UpDown <- 'not.sig'
  degs[degs$p_val_adj < 0.05 & degs$avg_log2FC >= 0.25 & degs$pct.1 >= 0.1, ]$UpDown <- 'Up'
  degs[degs$p_val_adj < 0.05 & degs$avg_log2FC <= -0.25 & degs$pct.1 >= 0.1, ]$UpDown <- 'Down'
  
  # 归档数据
  saveRDS(degs, file.path(outdir, 'DEGs.rds'))
  diffgene <- degs[degs$UpDown != 'not.sig', ]
  
  # 导出 Excel
  if(nrow(diffgene) > 0){
    cat("   - 导出显著基因明细至 Excel... \n")
    wb_list <- split(diffgene, diffgene$cluster)
    write.xlsx(wb_list, file.path(outdir, '1_DEGs_FDR0.05.xlsx'))
  }
  
  # --- 基础绘图 ---
  cat("   - 正在绘制可视化大图... \n")
  # Top10 DotPlot
  top10df <- diffgene %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  if (nrow(top10df) > 0) {
    wd <- max(10, length(unique(top10df$gene)) * 0.2 + 2)
    ht <- max(6, length(unique(degs$cluster)) * 0.5 + 2)
    p_dot <- DotPlot(data, features = unique(top10df$gene)) +
      scale_colour_gradientn(colours = c("dodgerblue1","#44c1f0", "lightgoldenrod","#e20612",'#cc340c')) +
      geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.2) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    ggsave(file.path(outdir, '2_Top10DEGs.DotPlot.pdf'), p_dot, width = wd, height = ht)
    ggsave(file.path(outdir, '2_Top10DEGs.DotPlot.png'), p_dot, width = wd, height = ht)
  }

  # Top5 FeaturePlot
  fe_dir <- file.path(outdir, "3_FeatureExprUMAP")
  dir.create(fe_dir, showWarnings = FALSE)
  top5df <- degs %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  for (g in head(unique(top5df$gene), 50)){ # 最多画 50 个，避免IO崩溃
    p_fe <- FeatureDimPlot(data, features = g, pt.size = 0.01, theme_use = "theme_blank")
    ggsave(file.path(fe_dir, paste0(g, ".png")), p_fe, width = 5, height = 4)
  }
  
  return(degs)
}

# --- [模式 B] 亚群内组间两两比较核心函数 ---
# input: Seurat对象, 亚群拆分列(如CellType), 比较列(如Group), 比较列表文件, 输出路径
group_de_by_subtype <- function(data, split_by, condition, cmp_file, outdir, assay='RNA'){
  cat("=> [DE] 启动模式 B: 亚群内组间精细比较 (Subtype Pairwise DE)...\n")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  min_cells_per_group <- 3
  
  if (!file.exists(cmp_file)) stop("错误: 找不到比较列表文件 (cmp_file)!")
  cmp_list <- read.table(cmp_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
  
  DefaultAssay(data) <- assay
  
  for (i in 1:nrow(cmp_list)){
    g1 <- as.character(cmp_list[i, 1])
    g2 <- as.character(cmp_list[i, 2])
    cat(sprintf("\n>>> 正在处理比较对: %s vs %s <<<\n", g1, g2))
    
    comp_name <- paste0(g1, ".vs.", g2)
    pair_dir <- file.path(outdir, comp_name)
    dir.create(pair_dir, showWarnings = FALSE)
    
    # 建立项目级的 Excel
    debook <- createWorkbook()
    
    # 遍历亚群
    subtypes <- sort(unique(as.character(data@meta.data[[split_by]])))
    for (st in subtypes){
      cat(sprintf("   - 挖掘子类: %s ", st))
      st_safe <- gsub('[ |+|:|/]', '_', st)
      
      # 提子集并初步检查
      cells_idx <- which(data@meta.data[[split_by]] == st)
      if (length(cells_idx) < 10) { cat("(过少, 跳过)\n"); next }
      
      sub_obj <- subset(data, cells = rownames(data@meta.data)[cells_idx])
      existing_groups <- unique(as.character(sub_obj@meta.data[[condition]]))
      
      if (!(g1 %in% existing_groups && g2 %in% existing_groups)){
        cat("(缺失组别数据, 跳过)\n")
        next
      }

      group_counts <- table(as.character(sub_obj@meta.data[[condition]]))
      n1 <- unname(group_counts[g1])
      n2 <- unname(group_counts[g2])
      n1 <- ifelse(length(n1) == 0 || is.na(n1), 0, n1)
      n2 <- ifelse(length(n2) == 0 || is.na(n2), 0, n2)

      if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
        cat(sprintf("(组细胞数不足, 跳过: %s=%d, %s=%d; 至少需要 %d)\n",
                    g1, n1, g2, n2, min_cells_per_group))
        next
      }
      
      # 开始差异分析
      cat("... 正在计算 \n")
      Idents(sub_obj) <- sub_obj@meta.data[[condition]]
      mks <- tryCatch(
        FindMarkers(sub_obj, ident.1 = g1, ident.2 = g2, min.pct = 0.1, logfc.threshold = 0),
        error = function(e) {
          cat(sprintf("(差异分析失败, 跳过: %s)\n", conditionMessage(e)))
          return(NULL)
        }
      )
      if (is.null(mks)) {
        next
      }
      mks$compare <- comp_name
      mks$gene <- rownames(mks)
      mks$subtype <- st
      
      # 判定
      mks$UpDown <- 'Not.sig'
      mks[mks$avg_log2FC >= 0.25 & mks$p_val_adj < 0.05, ]$UpDown <- 'Up'
      mks[mks$avg_log2FC <= -0.25 & mks$p_val_adj < 0.05, ]$UpDown <- 'Down'
      mks <- setorderv(mks, c("avg_log2FC", "pct.1", "pct.2"), c(-1, -1, 1))
      
      # 写出明细
      st_out_dir <- file.path(pair_dir, paste0(comp_name, "_", st_safe))
      dir.create(st_out_dir, showWarnings = FALSE)
      write.csv(mks, file.path(st_out_dir, "1_DEGs_all.csv"), row.names = FALSE, quote = FALSE)
      
      # 写入 Sheet
      sig_df <- mks[mks$UpDown != 'Not.sig', ]
      if (nrow(sig_df) > 0) {
        # 截断名字以符合Excel限制 (31字符)
        sheet_n <- paste0(substr(st_safe, 1, 25), "_", i)
        addWorksheet(debook, sheetName = sheet_n)
        writeData(debook, sheet = sheet_n, x = sig_df)
      }
    }
    saveWorkbook(debook, file.path(pair_dir, paste0(comp_name, "_DEGs_FDR0.05.xlsx")), overwrite = TRUE)
  }
}

# --- 常规富集分析控制函数 ---
cluster_enrich_auto <- function(derds, orgdb, organism_kegg, outdir){
  cat("=> [Enrich] 正在执行全自动化功能富集程序...\n")
  degs <- readRDS(derds)
  if (!("cluster" %in% colnames(degs)) && "subtype" %in% colnames(degs)) {
    # 如果是模式 B 后续由于输出结构分散，暂时由主脚本外接循环调用，内部不做多层循环处理
    return()
  }
  
  # 复用已定义的逻辑
  if (exists("cluster_enrich")) {
     cluster_enrich(derds = derds, orgdb = orgdb, organism_kegg = organism_kegg, 
                    enrichdir = file.path(outdir, "4_Enrich"))
  } else {
     warning("无法定位到基础富集函数 cluster_enrich，请检查 func_scRNA_celltype_anno.R 是否成功加载。")
  }
}
