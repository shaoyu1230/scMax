# /annoroad/data1/bioinfo/PMO/shaoyusong/01.software/installed/miniconda3/envs/SCP1/bin/Rscript
suppressMessages(library(Seurat))
suppressMessages(library(SCP))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))
suppressMessages(library(Cairo))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(tidyr))
suppressMessages(library(ggsci))
suppressMessages(library(gridExtra))
suppressMessages(library(openxlsx))
suppressMessages(library(getopt))
suppressMessages(library(clusterProfiler))
suppressMessages(library(configr))
bindir <- this.path::this.dir()
source(file.path(bindir,'func_scRNA_celltype_anno.R'))
parser = argparse::ArgumentParser(description = 'Script for cluster Plots, DE and GO KEGG enrichment analysis.')
parser$add_argument('-c', '--config', dest = 'config', help = 'config file')

opt = parser$parse_args()
config = read.config(opt$config)
config = config$Para
outdir = config$outdir
inputrds = config$inputrds
metacsv = config$metacsv
sample_col = config$sample_col
sample_order = config$sample_order
group_col = config$group_col
if(!is.null(group_col)){
  group_col <- strsplit(group_col,split=',')[[1]]
}
group_order = config$group_order

clustercol = config$cluster_col #分辨率所在列名
cluster_colors = config$cluster_colors
# celltypecolors = config$celltype_colors
refmarker.file = config$refmarker.file
orgdb = config$orgdb
organism_kegg= config$organism_kegg

########## 0_preparations ################
# 创建文件输出目录
if(!dir.exists(paste0(outdir,'/Rdata'))){dir.create(paste0(outdir,'/Rdata'),recursive = T)}
if(!dir.exists(paste0(outdir,'/tmp'))){dir.create(paste0(outdir,'/tmp'),recursive = T)}
if(!dir.exists(paste0(outdir,'/1_cluster_characterization'))){dir.create(paste0(outdir,'/1_cluster_characterization'))}

# 判断rds文件是否存在，不存在报错退出
if(!file.exists(inputrds)){
  stop(paste0(inputrds,' file not exists!'))
}

# 读入文件，记录时间
print(Sys.time())
data=readRDS(inputrds)
DefaultAssay(data)="RNA"
if(is.null(metacsv)){
  print('No meta file! use default metadata!')
}else if(file.exists(metacsv)){
  meta <- read.csv(metacsv,stringsAsFactors = F,row.names=1)
  data <- AddMetaData(data,meta)  #更新metadata
}else{
  stop(paste0(metacsv,' file not exists!'))
}

# 如果提供了sample_order，更改data中相应列为factor，并添加level
if(!is.null(sample_order)){
  sample_order=strsplit(sample_order,split=',')[[1]]
  data@meta.data[,sample_col] <- as.character(data@meta.data[,sample_col])
  data@meta.data[,sample_col]=factor(data@meta.data[,sample_col],levels=sample_order)
}
if (sample_col!='Sample'){
  data$Sample <- data@meta.data[,sample_col] #报告展示Sample
}

# 如果提供了group_order，更改data中相应列为factor，并添加level
if(!is.null(group_order) & !is.null(group_col)){
  gcol <- group_col[1]
  group_order=strsplit(group_order,split=',')[[1]]
  if (length(group_order[group_order %in% unique(data@meta.data[,gcol])]) != length(unique(data@meta.data[,gcol]))){
    stop(paste0('group_order not same as ',gcol,' in meta.data'))
  }
  data@meta.data[,gcol] <- as.character(data@meta.data[,gcol])
  data@meta.data[,gcol]=factor(data@meta.data[,gcol],levels=group_order)
}


if (!is.null(group_col)){
  if (group_col[1] != 'Group'){
    data$Group <- data@meta.data[,group_col[1]] #报告展示Group，不能有Sample
  }
  plot_groups = unique(c('Sample',group_col)) #sample和group列
}else{
  print("Group col is not set!" )
  plot_groups = 'Sample'
}
# 检查提供的sample_col和group_col是否在数据中存在
columns <- colnames(data@meta.data)
columns <- columns[columns %in% plot_groups]
print(paste0('plot_groups: ',plot_groups))
print(paste0('plot_groups in metadata: ',columns))
if(length(plot_groups)>length(columns)){
  stop(paste0('Not all plot_groups in the data!'))
}


if(!is.null(clustercol)){
  if(clustercol %in% colnames(data@meta.data)){
    data$seurat_clusters <- data@meta.data[,clustercol]
  }else{
    stop(paste0(clustercol, 'Not in the metadata!'))
  }
}
Idents(data) <- data$seurat_clusters
# 颜色
ncluster <- length(unique(data$seurat_clusters))
if(length(cluster_colors) < ncluster){
  print(paste0('cluster_colors length is not enough! use default colors!'))
  cluster_colors <- NULL
}
print(paste0(Sys.time(),' - data read in done!'))
##################### 1_cluster_characterization-DE,Enrich ##########################
## DEGs analysis
print(paste0(Sys.time(),' ############## - start 3_DEGs'))

# 默认计算seurat_cluster各群差异基因
library(SeuratWrappers)
diffgene <- cluster_de(data,rdsdir=file.path(outdir,'/Rdata'),
           outdir=file.path(outdir,'1_cluster_characterization/4_Cluster.DE'))
print(paste0(Sys.time(),' ############## - 3_DEGs done!'))
print(paste0("DEGs results in: ",outdir,'/1_cluster_characterization/4_Cluster.DE'))
## enrichment
print(paste0(Sys.time(),' ############## - start 4_Enrich'))
cluster_enrich(derds=file.path(outdir,'/Rdata/ClusterDEGs.rds'),
               orgdb=orgdb,organism_kegg=organism_kegg,
               enrichdir=paste0(outdir,'/1_cluster_characterization/5_Enrich'))
print(paste0(Sys.time(),' ############## - 4_Enrich done!'))
print(paste0("DEGs results in: ",outdir,'/1_cluster_characterization/5_Enrich'))






