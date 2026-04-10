suppressMessages(library(yaml))
suppressMessages(library(argparse))

script_dir <- this.path::this.dir()
source(file.path(script_dir, "func_scRNA_celltype_anno.R"))
source(file.path(script_dir, "func_MergeSub.R"))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

parser <- argparse::ArgumentParser(
  description = "scMax 06_merge_sub: 将亚群注释结果回填到大群对象，并输出层级 UMAP 与样本/分组统计"
)
parser$add_argument("-c", "--config", dest = "config", default = "", help = "YAML 配置文件路径")
parser$add_argument("-m", "--main_rds", dest = "main_rds", default = "", help = "主大群 RDS 路径")
parser$add_argument("-s", "--sub_map", dest = "sub_map", default = "", help = "亚群映射 (MajorA:pathA,MajorB:pathB)")
parser$add_argument("-o", "--outdir", dest = "outdir", default = "", help = "输出目录")
parser$add_argument("--major_col", dest = "major_col", default = "", help = "大群列名")
parser$add_argument("--sub_col", dest = "sub_col", default = "", help = "亚群 RDS 中的注释列名")
parser$add_argument("--subtype_col", dest = "subtype_col", default = "", help = "回填到主对象的新列名")
parser$add_argument("--sample_col", dest = "sample_col", default = "", help = "样本列名")
parser$add_argument("--group_col", dest = "group_col", default = "", help = "分组列名，多个用逗号分隔")
opt <- parser$parse_args()

conf_root <- list()
if (!is.null(opt$config) && opt$config != "") {
  conf_root <- yaml::read_yaml(opt$config)
}
merge_conf <- conf_root$`06_merge_sub` %||% conf_root$merge_sub %||% conf_root$`07_merge_sub` %||% list()

main_rds <- if (opt$main_rds != "") opt$main_rds else merge_conf$main_rds
outdir <- if (opt$outdir != "") opt$outdir else merge_conf$outdir %||% "./06_merge_sub_out"
major_col <- if (opt$major_col != "") opt$major_col else merge_conf$major_col %||% "MajorType"
sub_col <- if (opt$sub_col != "") opt$sub_col else merge_conf$sub_col %||% "CellType"
subtype_col <- if (opt$subtype_col != "") opt$subtype_col else merge_conf$subtype_col %||% "SubType"
sample_col <- if (opt$sample_col != "") opt$sample_col else merge_conf$sample_col %||% "Sample"
group_col <- if (opt$group_col != "") opt$group_col else merge_conf$group_col %||% "Group"

sub_list <- list()
if (!is.null(opt$sub_map) && opt$sub_map != "") {
  pairs <- strsplit(opt$sub_map, ",")[[1]]
  for (p in pairs) {
    kv <- strsplit(p, ":", fixed = TRUE)[[1]]
    if (length(kv) < 2) {
      stop(sprintf("sub_map 格式错误: %s", p))
    }
    sub_list[[trimws(kv[1])]] <- trimws(paste(kv[-1], collapse = ":"))
  }
} else if (!is.null(merge_conf$sub_rds_map)) {
  sub_list <- merge_conf$sub_rds_map
}

if (is.null(main_rds) || main_rds == "") {
  stop("必须指定 main_rds。")
}
if (!file.exists(main_rds)) {
  stop(sprintf("main_rds 不存在: %s", main_rds))
}
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

dir.create(file.path(outdir, "Rdata"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "summary"), recursive = TRUE, showWarnings = FALSE)

message(">>> 加载主对象: ", main_rds)
main_obj <- readRDS(main_rds)

message(">>> 开始合并亚群注释结果...")
merge_res <- merge_sub_labels(
  main_obj = main_obj,
  sub_list = sub_list,
  major_col = major_col,
  sub_label_col = sub_col,
  subtype_col = subtype_col
)
merged_obj <- merge_res$object

save_merge_summary(merge_res, file.path(outdir, "summary"))

plot_prep <- prepare_plot_groups(
  obj = merged_obj,
  sample_col = sample_col,
  group_col = group_col
)
merged_obj <- plot_prep$object
plot_groups <- plot_prep$groups

message(">>> 绘制大群 / 亚群层级 UMAP ...")
plot_info <- plot_hierarchy_umaps(
  obj = merged_obj,
  outdir = outdir,
  major_col = major_col,
  subtype_col = subtype_col,
  groups = plot_groups
)

message(">>> 统计亚群在不同样本 / 分组中的细胞数与比例 ...")
plot_subtype_statistics(
  obj = merged_obj,
  outdir = outdir,
  subtype_col = subtype_col,
  sample_col = sample_col,
  group_col = group_col,
  subtype_colors = plot_info$subtype_colors
)

message(">>> 保存整合后的 RDS ...")
saveRDS(merged_obj, file.path(outdir, "Rdata", "Data-Merged_SubType.rds"))

message(">>> 06_merge_sub 完成，结果目录: ", outdir)
