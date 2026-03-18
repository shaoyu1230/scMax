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
#source(file.path(bindir,'func_distribution_Roe.R'))
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
group_order = config$group_order
if(!is.null(group_col)){
  group_col <- strsplit(group_col,split=',')[[1]]
}
roe_groups = group_col
plot_groups = unique(c(sample_col,group_col)) 
celltype.levels = config$celltype.levels
if(!is.null(celltype.levels)){
  celltype.levels = strsplit(celltype.levels, "\\s*,\\s*")[[1]]#以逗号拆分字符串，忽略逗号前后的空格
  print("celltype levels:")
  print(celltype.levels)
}
celltype.keep = config$celltype.keep
celltype.keep = strsplit(celltype.keep, "\\s*,\\s*")[[1]]
rdsprefix = config$rdsprefix
annofile = config$annofile # celltype annotation file
clustercol = config$cluster_col #分辨率所在列名
celltype_col = config$celltype_col # 细胞类型注释列名
celltype_colors = config$celltype_colors # 细胞类型注释颜色
if(!is.null(celltype_colors)){
	celltype_colors = strsplit(celltype_colors,"\\s*,\\s*")[[1]]
}
orgdb = config$orgdb
organism_kegg= config$organism_kegg

########## 0_preparations ################
# 创建文件输出目录
if(!dir.exists(paste0(outdir,'/Rdata'))){dir.create(paste0(outdir,'/Rdata'),recursive = T)}
if(!dir.exists(paste0(outdir,'/tmp'))){dir.create(paste0(outdir,'/tmp'),recursive = T)}
if(!dir.exists(paste0(outdir,'/3_annotation'))){dir.create(paste0(outdir,'/3_annotation'))}
if(!dir.exists(paste0(outdir,'/4_celltype_fraction'))){dir.create(paste0(outdir,'/4_celltype_fraction'))}
if(!dir.exists(paste0(outdir,'/5_celltype_characterization'))){dir.create(paste0(outdir,'/5_celltype_characterization'))}

# 判断文件是否存在，不存在报错退出
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
  meta <- read.csv(metacsv,stringsAsFactors = F, row.names=1)
  data <- AddMetaData(data,meta)  #更新metadata
}else{
  stop(paste0(metacsv,' file not exists!'))
}

# 保存一份原始metadata
metaraw <- data@meta.data

# 如果提供了sample_order，更改data中相应列为factor，并添加level
if(!is.null(sample_order)){
  sample_order=strsplit(sample_order,split=',')[[1]]
  if (length(sample_order[sample_order %in% unique(data@meta.data[,sample_col])]) != length(unique(data@meta.data[,sample_col]))){
    stop(paste0('sample_order not same as samples!'))
  }
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
print(paste0(Sys.time(),' - data read in done!'))

if(is.null(rdsprefix)){
  rdsprefix <- 'Annotation'
}

print(paste0(Sys.time(),' ############## - Part2: CellType Annotation !'))

##################### 3_annotation ##########################
print(paste0(Sys.time(),' ############## - start 3_annotation'))
## add annotation
data <- add_anno(data,outdir=paste0(outdir,'/3_annotation'),annofile,
                        celltype.levels,celltype_col=celltype_col)
metaraw[,celltype_col]=data@meta.data[,celltype_col] #原meta添加细胞类型注释
# 颜色
ncluster <- length(unique(data@meta.data[,celltype_col]))
if(length(celltype_colors)>=ncluster){
  celltype_colors <- celltype_colors
}else{
  print(paste0('celltype_colors length is not enough! use default colors!'))
  celltype_colors <- palette_scp(unique(data@meta.data[,celltype_col]),palette='Paired',NA_keep=TRUE)
}
celltype_umap_plots(data,outdir=paste0(outdir,'/3_annotation/CellType_All'),celltype_col=celltype_col,groups=plot_groups,palcolor=celltype_colors)
library(gridExtra)
df <- read.table(annofile,sep='\t',header=TRUE)
df <- df[,c('Cluster','CellType','Markers')]
png(paste0(outdir,'/3_annotation/CellType_All/0_Cluster-CellType.anno.png'), width = 2000, height = 1600, res = 150)
grid.table(df)
dev.off()
file.copy(annofile,paste0(outdir,'/3_annotation/CellType_All/0_Cluster-CellType.anno.xls'))
## feature plot
anno.df <- read.csv(annofile,sep='\t')
if('Markers' %in% colnames(anno.df)){
  anno_mks_plot1(data,annofile,celltype_col=celltype_col,celltype.levels,outdir=paste0(outdir,'/3_annotation/CellType_All'))
}

## subset celltype
keep.cells <- rownames(data@meta.data[data@meta.data[,celltype_col] %in% celltype.keep,])
data.f <- subset(data,cells=keep.cells)
data.f@meta.data[,celltype_col] <- droplevels(data.f@meta.data[,celltype_col])
metaraw.f <- metaraw[keep.cells,]
celltype_umap_plots(data.f,outdir=paste0(outdir,'/3_annotation/CellType_Keep'),celltype_col=celltype_col,groups=plot_groups,palcolor=celltype_colors) #1,2

## feature plot
if('Markers' %in% colnames(anno.df)){
  anno_mks_plot2(data.f,annofile,celltype_col=celltype_col,celltype.keep,outdir=paste0(outdir,'/3_annotation/CellType_Keep'))
}
saveRDS(data,paste0(outdir,'/Rdata/Data-',rdsprefix,'_CellType.rds'))
saveRDS(data.f,paste0(outdir,'/Rdata/Data-',rdsprefix,'_CellType.Keep.rds'))
write.csv(metaraw.f,paste0(outdir,'/Rdata/Meta-',rdsprefix,'_CellType.Keep.csv'),row.names = F)
write.csv(metaraw,paste0(outdir,'/Rdata/Meta-',rdsprefix,'_CellType.csv'),row.names = F)
print(paste0(Sys.time(),' ############## - 3_annotation done!'))

##################### 4_celltype_fraction ##########################
cell_fraction_plots(data,out_prefix=paste0(outdir,'/4_celltype_fraction/1_CellType_All.in.'),
                    celltype_col=celltype_col,groups=plot_groups,do.boxplot=FALSE,palcolor=celltype_colors)
cell_fraction_plots(data.f,out_prefix=paste0(outdir,'/4_celltype_fraction/2_CellType_Keep.in.'),
                    celltype_col=celltype_col,groups=plot_groups,do.boxplot=FALSE,palcolor=celltype_colors)
## Ro/e plots
roe_groups <- intersect(roe_groups,colnames(data.f@meta.data))
if (!is.null(roe_groups)){
  for(group in roe_groups){
    Roe <- distribution_Roe(data.f@meta.data,celltype_column = celltype_col,celltype_level = celltype.keep,
                            condition_column = group,out_prefix=paste0(outdir,'/4_celltype_fraction/3_CellType_Keep.in.',group,'-'))
    ggsave(paste0(outdir,'/4_celltype_fraction/3_CellType_Keep.in.',group,'-Roe.pdf'),width=5,height=5)
    ggsave(paste0(outdir,'/4_celltype_fraction/3_CellType_Keep.in.',group,'-Roe.png'),width=5,height=5)
  }
}
print(paste0(Sys.time(),' ############## - 4_celltype_fraction done!'))

##################### 5_celltype_characterization ##########################
print(paste0(Sys.time(),' ############## - start 5_celltype_characterization'))
library(SeuratWrappers)
diffgene <- celltype_de(data.f,rdsdir=file.path(outdir,'/Rdata'),celltype_col=celltype_col,outdir=paste0(outdir,'/5_celltype_characterization/1_CellType_DE'))
files <- list.files(paste0(outdir,'/5_celltype_characterization/1_CellType_DE'),pattern='4_FeatureExpr')
example.file <- files[grep('png',files)][1]
file.copy(paste0(outdir,'/5_celltype_characterization/1_CellType_DE/',example.file),paste0(outdir,'/5_celltype_characterization/1_CellType_DE/4_FeatureExprUMAP.example.png'),)

celltype_enrich(derds=paste0(outdir,'/Rdata/CellTypeDEGs.rds'),
                orgdb=orgdb,organism_kegg=organism_kegg,
                enrichdir=paste0(outdir,'/5_celltype_characterization/2_CellType_Enrich'))
print(paste0(Sys.time(),' ############## - CellType Annotation Part2 Done!'))




