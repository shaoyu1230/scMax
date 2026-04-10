suppressMessages(library(Seurat))
suppressMessages(library(SCP))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

assert_meta_col <- function(obj, colname, obj_name = "Seurat object") {
  if (is.null(colname) || colname == "" || !(colname %in% colnames(obj@meta.data))) {
    stop(sprintf("%s metadata 中缺少列: %s", obj_name, colname))
  }
}

split_config_values <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || x == "") {
    return(character(0))
  }
  unique(trimws(unlist(strsplit(as.character(x), "\\s*,\\s*"))))
}

safe_factor <- function(x, numeric_sort = FALSE) {
  x <- as.character(x)
  lvls <- unique(x)
  if (numeric_sort) {
    lvls <- stringr::str_sort(lvls, numeric = TRUE)
  }
  factor(x, levels = lvls)
}

merge_sub_labels <- function(main_obj,
                             sub_list,
                             major_col = "MajorType",
                             sub_label_col = "CellType",
                             subtype_col = "SubType") {
  assert_meta_col(main_obj, major_col, "main_obj")

  main_cells <- colnames(main_obj)
  main_major <- as.character(main_obj@meta.data[[major_col]])
  names(main_major) <- main_cells

  subtype_vec <- main_major
  names(subtype_vec) <- main_cells

  if (length(sub_list) == 0) {
    main_obj@meta.data[[subtype_col]] <- safe_factor(subtype_vec)
    summary_df <- data.frame(
      major_name = character(0),
      main_major_cells = integer(0),
      sub_rds_cells = integer(0),
      kept_cells = integer(0),
      removed_cells = integer(0),
      extra_cells_not_in_main = integer(0),
      stringsAsFactors = FALSE
    )
    return(list(
      object = main_obj,
      kept_cells = main_cells,
      removed_cells = character(0),
      summary = summary_df
    ))
  }

  if (is.null(names(sub_list)) || any(names(sub_list) == "")) {
    stop("sub_list 必须是具名列表，名称为大群名，值为对应亚群 RDS 路径。")
  }

  refined_majors <- unique(names(sub_list))
  untouched_cells <- main_cells[!main_major %in% refined_majors]
  kept_cells <- untouched_cells
  removed_cells <- character(0)
  summary_list <- list()

  message(">>> 开始处理亚群回填...")

  for (major_name in refined_majors) {
    rds_path <- sub_list[[major_name]]
    major_cells <- main_cells[main_major == major_name]

    if (length(major_cells) == 0) {
      warning(sprintf("主对象中不存在大群 [%s]，跳过。", major_name))
      next
    }
    if (is.null(rds_path) || rds_path == "" || !file.exists(rds_path)) {
      stop(sprintf("大群 [%s] 的亚群 RDS 不存在: %s", major_name, rds_path))
    }

    message(sprintf("--- 处理 [%s] <- %s", major_name, rds_path))
    sub_obj <- readRDS(rds_path)
    assert_meta_col(sub_obj, sub_label_col, sprintf("sub_obj (%s)", major_name))

    sub_cells <- colnames(sub_obj)
    retained_cells <- intersect(major_cells, sub_cells)
    missing_cells <- setdiff(major_cells, retained_cells)
    extra_cells <- setdiff(sub_cells, main_cells)

    if (length(retained_cells) == 0) {
      warning(sprintf("大群 [%s] 在亚群 RDS 中没有匹配到任何细胞，将全部视为删除。", major_name))
    } else {
      subtype_vec[retained_cells] <- as.character(sub_obj@meta.data[retained_cells, sub_label_col])
      kept_cells <- c(kept_cells, retained_cells)
    }

    removed_cells <- c(removed_cells, missing_cells)

    summary_list[[major_name]] <- data.frame(
      major_name = major_name,
      main_major_cells = length(major_cells),
      sub_rds_cells = length(sub_cells),
      kept_cells = length(retained_cells),
      removed_cells = length(missing_cells),
      extra_cells_not_in_main = length(extra_cells),
      stringsAsFactors = FALSE
    )

    if (length(extra_cells) > 0) {
      warning(sprintf("亚群 RDS [%s] 中有 %d 个细胞不在主对象中，已忽略。", major_name, length(extra_cells)))
    }
  }

  kept_cells <- unique(kept_cells)
  removed_cells <- unique(removed_cells)

  main_obj@meta.data[[subtype_col]] <- safe_factor(subtype_vec)
  message(sprintf(">>> 原始细胞数: %d", ncol(main_obj)))
  main_obj <- subset(main_obj, cells = kept_cells)
  main_obj@meta.data[[major_col]] <- safe_factor(main_obj@meta.data[[major_col]])
  main_obj@meta.data[[subtype_col]] <- safe_factor(main_obj@meta.data[[subtype_col]])
  message(sprintf(">>> 合并后保留细胞数: %d", ncol(main_obj)))
  message(sprintf(">>> 因亚群分析被删除的细胞数: %d", length(removed_cells)))

  summary_df <- dplyr::bind_rows(summary_list)

  list(
    object = main_obj,
    kept_cells = kept_cells,
    removed_cells = removed_cells,
    summary = summary_df
  )
}

save_merge_summary <- function(merge_res, outdir) {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  if (!is.null(merge_res$summary) && nrow(merge_res$summary) > 0) {
    write.csv(merge_res$summary,
              file.path(outdir, "merge_sub.summary.csv"),
              row.names = FALSE)
  }

  removed_df <- data.frame(cell = merge_res$removed_cells, stringsAsFactors = FALSE)
  write.csv(removed_df,
            file.path(outdir, "merge_sub.removed_cells.csv"),
            row.names = FALSE)

  invisible(TRUE)
}

prepare_plot_groups <- function(obj, sample_col = "Sample", group_col = "Group") {
  plot_groups <- character(0)

  assert_meta_col(obj, sample_col, "merged_obj")
  if (sample_col != "Sample") {
    obj$Sample <- obj@meta.data[[sample_col]]
  }
  plot_groups <- c(plot_groups, "Sample")

  group_cols <- split_config_values(group_col)
  if (length(group_cols) > 0) {
    for (gcol in group_cols) {
      assert_meta_col(obj, gcol, "merged_obj")
      if (gcol != "Group" && !"Group" %in% plot_groups) {
        obj$Group <- obj@meta.data[[gcol]]
      }
      plot_groups <- c(plot_groups, gcol)
    }
  }

  list(
    object = obj,
    groups = unique(plot_groups)
  )
}

plot_hierarchy_umaps <- function(obj,
                                 outdir,
                                 major_col = "MajorType",
                                 subtype_col = "SubType",
                                 groups = c("Sample", "Group")) {
  annotation_dir <- file.path(outdir, "annotation")
  major_dir <- file.path(annotation_dir, "MajorType")
  subtype_dir <- file.path(annotation_dir, "SubType")
  dir.create(major_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(subtype_dir, recursive = TRUE, showWarnings = FALSE)

  major_colors <- palette_scp(unique(as.character(obj@meta.data[[major_col]])),
                              palette = "Paired", NA_keep = TRUE)
  subtype_colors <- palette_scp(unique(as.character(obj@meta.data[[subtype_col]])),
                                palette = "Paired", NA_keep = TRUE)

  celltype_umap_plots(obj,
                      outdir = major_dir,
                      celltype_col = major_col,
                      groups = groups,
                      palcolor = major_colors)

  celltype_umap_plots(obj,
                      outdir = subtype_dir,
                      celltype_col = subtype_col,
                      groups = groups,
                      palcolor = subtype_colors)

  invisible(list(major_colors = major_colors, subtype_colors = subtype_colors))
}

plot_subtype_statistics <- function(obj,
                                    outdir,
                                    subtype_col = "SubType",
                                    sample_col = "Sample",
                                    group_col = "Group",
                                    subtype_colors = NULL) {
  stat_dir <- file.path(outdir, "subtype_fraction")
  dir.create(stat_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(subtype_colors) || length(subtype_colors) < length(unique(obj@meta.data[[subtype_col]]))) {
    subtype_colors <- palette_scp(unique(as.character(obj@meta.data[[subtype_col]])),
                                  palette = "Paired", NA_keep = TRUE)
  }

  stat_groups <- unique(c(sample_col, split_config_values(group_col)))
  stat_groups <- stat_groups[stat_groups %in% colnames(obj@meta.data)]

  if (length(stat_groups) == 0) {
    warning("未找到可用于统计的样本/分组列，跳过 subtype_fraction 绘图。")
    return(invisible(NULL))
  }

  cell_fraction_plots(obj,
                      out_prefix = file.path(stat_dir, "1_SubType.in."),
                      celltype_col = subtype_col,
                      groups = stat_groups,
                      do.boxplot = FALSE,
                      palcolor = subtype_colors)

  group_cols <- split_config_values(group_col)
  for (gcol in group_cols) {
    if (!(gcol %in% colnames(obj@meta.data))) {
      next
    }
    if (length(unique(obj@meta.data[[gcol]])) <= 1) {
      message(sprintf(">>> 跳过 Roe 分析 [%s]: 只有 1 个分组。", gcol))
      next
    }
    if (!("Sample" %in% colnames(obj@meta.data)) || length(unique(obj@meta.data$Sample)) == length(unique(obj@meta.data[[gcol]]))) {
      message(sprintf(">>> 跳过 Roe 分析 [%s]: 缺少样本重复。", gcol))
      next
    }

    roe_plot <- distribution_Roe(
      meta_data = obj@meta.data,
      celltype_column = subtype_col,
      celltype_level = levels(safe_factor(obj@meta.data[[subtype_col]])),
      condition_column = gcol,
      out_prefix = file.path(stat_dir, paste0("2_SubType.in.", gcol, "-"))
    )
    ggsave(file.path(stat_dir, paste0("2_SubType.in.", gcol, "-Roe.pdf")), roe_plot, width = 5, height = 5)
    ggsave(file.path(stat_dir, paste0("2_SubType.in.", gcol, "-Roe.png")), roe_plot, width = 5, height = 5)
  }

  invisible(TRUE)
}
