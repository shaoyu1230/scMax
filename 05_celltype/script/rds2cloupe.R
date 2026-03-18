#/annogene/data2/share/software/SCV/SC_tools/Miniforge3/envs/R4.2_seurat/bin/Rscript

# 加载必要的库
library(optparse)
library(loupeR)
default_executable_path = "/annogene/data2/share/software/SCV/SC_tools/Miniforge3/envs/R4.2_seurat/lib/R/library/loupeR/exec/louper-linux-x64"
# 创建命令行选项
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "输入Seurat对象的RDS文件路径", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./",
              help = "输出目录 [默认=%default]", metavar = "character"),
  make_option(c("-n", "--output_name"), type = "character", default = "LoupeOutput",
              help = "输出文件名称 [默认=%default]", metavar = "character"),
  make_option(c("-l", "--louper_path"), type = "character", default = "",
              help = "louper 可执行文件路径（可选）", metavar = "character"),
  make_option(c("-f", "--force"), type = "logical", default = TRUE,
              help = "是否覆盖已存在的文件 [默认=%default]", metavar = "logical")
)

# 解析参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("必须提供输入文件路径", call. = FALSE)
}

# 读取Seurat对象
cat("正在读取Seurat对象...\n")
obj <- readRDS(opt$input)

# 选择 louper 可执行路径：参数 > 环境变量 > 默认
executable_path <- opt$louper_path
if (is.null(executable_path) || executable_path == "") {
  env_path <- Sys.getenv("LOUPER_PATH", "")
  executable_path <- ifelse(env_path == "", default_executable_path, env_path)
}

if (!file.exists(executable_path)) {
  cat("警告：louper 可执行文件不存在，仍尝试调用 loupeR 默认机制：", executable_path, "\n")
}

# 创建Loupe可视化文件
cat("正在创建Loupe可视化文件...\n")
create_loupe_from_seurat(
  obj,
  output_dir = opt$output_dir,
  output_name = opt$output_name,
  force = opt$force,
  executable_path = executable_path
)

cat(paste0("完成! 输出文件已保存至: ", opt$output_dir, "/", opt$output_name, "\n"))
