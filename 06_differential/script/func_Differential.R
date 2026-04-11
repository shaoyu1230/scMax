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
suppressMessages(library(SeuratWrappers))

plotgo <- function(Markers_GO_enrich, outdir = "./", prefix = "",
                   titlei = "GO Enrichment") {
  CPCOLS <- c("#6495ED", "#8FBC8F", "#F4A460")
  mtx <- Markers_GO_enrich@result
  mtx$plt <- -log10(mtx$pvalue)
  dmtx <- dplyr::group_by(mtx, ONTOLOGY) %>% dplyr::top_n(10, plt)
  dorder <- factor(rev(as.integer(rownames(dmtx))), labels = rev(dmtx$Description))

  pbar <- ggplot(dmtx, aes(x = Description, y = plt, fill = ONTOLOGY)) +
    geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.5, aes(x = dorder)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 70)) +
    coord_flip() +
    scale_y_log10(breaks = c(1, 10, 100, 1000)) +
    theme(panel.background = element_rect(fill = "transparent", colour = NA)) +
    xlab("GO_Term") + ylab(expression(-"log"["10"] * "(PValue)")) + theme_bw() +
    scale_fill_manual(values = CPCOLS) + labs(title = titlei)

  pdot <- ggplot(dmtx, aes(x = plt, y = dorder)) +
    geom_point(aes(size = Count, color = -1 * log(pvalue), shape = ONTOLOGY)) +
    scale_color_gradient(low = "steelblue", high = "indianred") +
    scale_y_discrete(labels = function(y) stringr::str_wrap(y, width = 70)) +
    labs(color = expression(-log[10](pvalue)), size = "Count",
         x = expression(-"log"["10"] * "(PValue)"), y = "Go_term", title = titlei) +
    theme_bw()

  ggsave(paste0(outdir, "/", prefix, "GO_bar.pdf"), pbar, width = 8, height = 8)
  ggsave(paste0(outdir, "/", prefix, "GO_bar.png"), pbar, width = 8, height = 8)
  ggsave(paste0(outdir, "/", prefix, "GO_dot.png"), pdot, width = 8, height = 8)
  ggsave(paste0(outdir, "/", prefix, "GO_dot.pdf"), pdot, width = 8, height = 8)
}

plotkegg <- function(keggresultfile, outdir = "./", prefix = "", titlei = "") {
  mtx <- openxlsx::read.xlsx(keggresultfile)
  mtx$plt <- -log10(mtx$pvalue)
  if (nrow(mtx) >= 30) {
    dmtx <- dplyr::top_n(mtx, 30, plt)
  } else {
    dmtx <- mtx
  }
  dorder <- factor(rev(as.integer(rownames(dmtx))), labels = rev(dmtx$Description))

  pbar <- ggplot(dmtx, aes(x = Description, y = plt, fill = plt)) +
    geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.5, aes(x = dorder)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 60)) +
    coord_flip() +
    theme(panel.background = element_rect(fill = "transparent", colour = NA)) +
    xlab("Term") + ylab(expression(-"log"["10"] * "(PValue)")) +
    theme_bw() +
    scale_fill_gradientn(colours = c("steelblue", "yellow", "red")) +
    labs(title = titlei, fill = expression(-log[10](pvalue)))

  pdot <- ggplot(dmtx, aes(x = plt, y = dorder)) +
    geom_point(aes(size = Count, color = -1 * log(pvalue))) +
    scale_color_gradient(low = "steelblue", high = "indianred") +
    scale_y_discrete(labels = function(y) stringr::str_wrap(y, width = 70)) +
    labs(color = expression(-log[10](pvalue)), size = "Count",
         x = expression(-"log"["10"] * "(PValue)"), y = "KEGG_term", title = titlei) +
    theme_bw()

  ggsave(paste0(outdir, "/", prefix, "KEGG_bar.pdf"), pbar, width = 8, height = 8)
  ggsave(paste0(outdir, "/", prefix, "KEGG_bar.png"), pbar, width = 8, height = 8)
  ggsave(paste0(outdir, "/", prefix, "KEGG_dot.pdf"), pdot, width = 8, height = 8)
  ggsave(paste0(outdir, "/", prefix, "KEGG_dot.png"), pdot, width = 8, height = 8)
}

GeneList2GOKEGG <- function(infile, orgdb = "org.Rn.eg.db", organism_kegg = "rno",
                            outdir, prefix = "", runGO = TRUE, runKEGG = TRUE, title = "") {
  suppressMessages(library(clusterProfiler))
  suppressMessages(library(KEGG.db))

  intersectGenes <- read.table(infile)[, 1]

  if (runGO) {
    Markers_GO_enrich <- enrichGO(
      intersectGenes,
      OrgDb = orgdb,
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      keyType = "SYMBOL",
      readable = TRUE
    )
    if (!is.null(Markers_GO_enrich) && nrow(Markers_GO_enrich@result) > 0) {
      openxlsx::write.xlsx(Markers_GO_enrich@result, paste0(outdir, "/", prefix, "GO.xlsx"))
      plotgo(Markers_GO_enrich, outdir = outdir, prefix = prefix, titlei = paste0(title, "GO Enrichment"))
    }
  }

  if (runKEGG) {
    eg <- bitr(intersectGenes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = orgdb)
    genelist <- eg$ENTREZID[!is.na(eg$ENTREZID)]
    keggresult <- enrichKEGG(
      gene = genelist,
      use_internal_data = TRUE,
      organism = organism_kegg,
      keyType = "kegg",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    if (!is.null(keggresult)) {
      keggresult <- setReadable(keggresult, orgdb, keyType = "ENTREZID")
      openxlsx::write.xlsx(keggresult@result, paste0(outdir, "/", prefix, "KEGG.xlsx"))
      plotkegg(
        keggresultfile = paste0(outdir, "/", prefix, "KEGG.xlsx"),
        outdir = outdir,
        prefix = prefix,
        titlei = paste0(title, "KEGG Enrichment")
      )
    }
  }
}

plot_group_volcano <- function(deg_df, g1, g2, out_prefix, top_n = 8) {
  df <- as.data.frame(deg_df)
  if (nrow(df) == 0) {
    return(invisible(NULL))
  }

  df$gene_label <- ""
  df$neglog10_padj <- -log10(pmax(df$p_val_adj, 1e-300))
  df$neglog10_padj[!is.finite(df$neglog10_padj)] <- max(df$neglog10_padj[is.finite(df$neglog10_padj)], na.rm = TRUE)
  if (!all(is.finite(df$neglog10_padj))) {
    df$neglog10_padj[!is.finite(df$neglog10_padj)] <- 0
  }

  sig_up <- df[df$UpDown == "Up", , drop = FALSE]
  sig_down <- df[df$UpDown == "Down", , drop = FALSE]
  label_df <- rbind(
    head(sig_up[order(-sig_up$avg_log2FC, -sig_up$neglog10_padj), , drop = FALSE], top_n),
    head(sig_down[order(sig_down$avg_log2FC, -sig_down$neglog10_padj), , drop = FALSE], top_n)
  )
  label_genes <- unique(as.character(label_df$gene))
  df$gene_label[df$gene %in% label_genes] <- df$gene[df$gene %in% label_genes]

  p <- ggplot(df, aes(x = avg_log2FC, y = neglog10_padj, color = UpDown)) +
    geom_point(size = 1.1, alpha = 0.75) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
    scale_color_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "Not.sig" = "grey75")) +
    labs(
      title = paste0(g1, " vs ", g2),
      x = "avg_log2FC",
      y = expression(-log[10]("adj.P"))
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())

  if (any(df$gene_label != "")) {
    p <- p + geom_text(
      data = df[df$gene_label != "", , drop = FALSE],
      aes(label = gene_label),
      size = 2.8,
      check_overlap = TRUE,
      vjust = -0.2
    )
  }

  ggsave(paste0(out_prefix, ".pdf"), p, width = 6, height = 5)
  ggsave(paste0(out_prefix, ".png"), p, width = 6, height = 5)
  invisible(p)
}

run_group_enrich <- function(deg_df, g1, g2, subtype_name, outdir, orgdb, organism_kegg) {
  enrich_dir <- file.path(outdir, "3_Enrich")
  dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)

  up_genes <- unique(as.character(deg_df$gene[deg_df$UpDown == "Up"]))
  down_genes <- unique(as.character(deg_df$gene[deg_df$UpDown == "Down"]))
  gene_sets <- list(Up = up_genes, Down = down_genes)

  for (direction in names(gene_sets)) {
    genes <- gene_sets[[direction]]
    if (length(genes) < 5) {
      message(sprintf("      富集跳过 [%s - %s]: 基因数过少 (%d)", subtype_name, direction, length(genes)))
      next
    }

    gene_file <- file.path(enrich_dir, paste0(direction, ".genes.txt"))
    write.table(genes, gene_file, row.names = FALSE, col.names = FALSE, quote = FALSE)

    prefix <- paste0(direction, ".")
    title <- paste0(subtype_name, " ", g1, " vs ", g2, " ", direction, " ")
    tryCatch(
      GeneList2GOKEGG(
        infile = gene_file,
        orgdb = orgdb,
        organism_kegg = organism_kegg,
        outdir = enrich_dir,
        prefix = prefix,
        runGO = TRUE,
        runKEGG = TRUE,
        title = title
      ),
      error = function(e) {
        warning(sprintf("富集分析失败 [%s - %s]: %s", subtype_name, direction, conditionMessage(e)))
      }
    )
  }

  invisible(TRUE)
}

cluster_enrich <- function(derds, orgdb = "org.Hs.eg.db", organism_kegg = "hsa",
                           enrichdir = "Enrich") {
  suppressMessages(library(KEGG.db))
  if (!dir.exists(enrichdir)) {
    dir.create(enrichdir, recursive = TRUE)
  }

  degs <- readRDS(derds)
  clusters <- unique(degs$cluster)
  degs.up <- degs[degs$p_val_adj < 0.05 & degs$avg_log2FC >= 0.25 & degs$pct.1 >= 0.1, ]
  tmpdir <- file.path(enrichdir, "..", "..", "tmp")
  if (!dir.exists(tmpdir)) {
    dir.create(tmpdir, recursive = TRUE)
  }

  for (c in clusters) {
    genes <- degs.up[degs.up$cluster == c, ]$gene
    write.table(genes, file.path(tmpdir, paste0(c, ".up.txt")),
                row.names = FALSE, quote = FALSE, col.names = FALSE)
  }

  print(paste0(date(), " - Enrichment preparation finish!"))
  file.list <- list.files(tmpdir, pattern = "up.txt$")
  print(paste0(date(), " - Enrichment analysis starts!"))
  for (file in file.list) {
    prefix <- gsub(".txt$", "", file)
    infile <- file.path(tmpdir, file)
    GeneList2GOKEGG(
      infile = infile,
      orgdb = orgdb,
      organism_kegg = organism_kegg,
      outdir = enrichdir,
      prefix = paste0(prefix, "."),
      runGO = TRUE,
      runKEGG = TRUE,
      title = paste0(prefix, " ")
    )
  }
  print(paste0(date(), " - Enrichment analysis finish!"))
}

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
group_de_by_subtype <- function(data, split_by, condition, cmp_file, outdir, assay='RNA',
                                do_enrich = FALSE, orgdb = "org.Hs.eg.db", organism_kegg = "hsa"){
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
      mks <- tryCatch(
        RunPresto(
          sub_obj,
          group.by = condition,
          ident.1 = g1,
          ident.2 = g2,
          logfc.threshold = 0,
          min.pct = 0.1,
          min.cells.group = min_cells_per_group
        ),
        error = function(e) {
          cat(sprintf("(差异分析失败, 跳过: %s)\n", conditionMessage(e)))
          return(NULL)
        }
      )
      if (is.null(mks)) {
        next
      }
      mks <- as.data.frame(mks)
      if (nrow(mks) == 0) {
        cat("(未检出差异基因, 跳过)\n")
        next
      }
      mks$compare <- comp_name
      if (!("gene" %in% colnames(mks))) {
        mks$gene <- rownames(mks)
      }
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
      sig_df <- mks[mks$UpDown != 'Not.sig', ]
      write.csv(sig_df, file.path(st_out_dir, "1_DEGs_FDR0.05.csv"), row.names = FALSE, quote = FALSE)
      plot_group_volcano(mks, g1 = g1, g2 = g2, out_prefix = file.path(st_out_dir, "2_Volcano"))
      
      if (isTRUE(do_enrich)) {
        run_group_enrich(
          deg_df = mks,
          g1 = g1,
          g2 = g2,
          subtype_name = st,
          outdir = st_out_dir,
          orgdb = orgdb,
          organism_kegg = organism_kegg
        )
      }

      # 写入 Sheet
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
  
  cluster_enrich(
    derds = derds,
    orgdb = orgdb,
    organism_kegg = organism_kegg,
    enrichdir = file.path(outdir, "4_Enrich")
  )
}
