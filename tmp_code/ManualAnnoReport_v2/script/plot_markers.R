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

######### plot_refmarkers ############
# 功能：根据参考marker文件绘制DotPlot和FeaturePlot
# 参考marker文件格式：第一列为CellType，第二列为Markers，Markers之间用逗号或空格分隔，列名用Tab分隔
plot_refmarkers_dot <- function(data,refmarker.file,outdir=paste0(outdir,'/2_marker_expression')){
  if (!dir.exists(outdir)){dir.create(outdir,recursive=T)}
  ####### read markers ########
  mks_ref <- read.csv(refmarker.file,sep='\t')
  df_long <- mks_ref[c('CellType','Markers')] %>%  separate_rows(Markers, sep = "[, ]+")
  df_unique <- df_long %>% unique() %>% data.frame() # 转换为数据框
  df_unique$Markers <- gsub(' ','',df_unique$Markers)
  ####### 1_DotPlot ########
  ht <- GroupHeatmap(data,features = df_unique$Markers,
                     feature_split= df_unique$CellType,assay = "RNA",slot = 'data',
                     group.by = "seurat_clusters",heatmap_palette = "RdBu",
                     show_row_names = TRUE, show_column_names = T,column_names_rot = 0,cluster_columns = T,
                     anno_terms = FALSE,nlabel=0,add_dot = TRUE,add_bg = TRUE)
  ht$plot
  ncluster=length(unique(data$seurat_clusters))
  wd=max(8,ncluster*0.6)
  nmarkers=nrow(df_unique)
  ht=min(max(8,nmarkers*0.5),30)
  ggsave(file.path(outdir,"dotplot1.pdf"),width=wd,height=ht)
  ggsave(file.path(outdir,"dotplot1.png"),width=wd,height=ht)

  # dotplot
  DotPlot(data,features =unique(df_unique$Markers),group.by="seurat_clusters")+ #coord_flip()+
    scale_colour_gradientn(colours = c("dodgerblue1","#44c1f0", "lightgoldenrod","#e20612",'#cc340c'))+
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2)+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.title=element_blank())
  mks <- unique(df_unique$Markers)
  wd <- max(ceiling(length(mks)/3),9)
  ht <- max(length(unique(data@meta.data[,celltype_col]))*0.4,5)
  ggsave(file.path(outdir,"dotplot2.png"),width=wd,height=ht)
  ggsave(file.path(outdir,"dotplot2.pdf"),width=wd,height=ht)
}
plot_refmarkers_umap <- function(data,refmarker.file,outdir=paste0(outdir,'/2_marker_expression')){
  if (!dir.exists(outdir)){dir.create(outdir,recursive=T)}
  ####### read markers ########
  mks_ref <- read.csv(refmarker.file,sep='\t')
  df_long <- mks_ref[c('CellType','Markers')] %>%  separate_rows(Markers, sep = "[, ]+")
  df_unique <- df_long %>% unique() %>% data.frame() # 转换为数据框
  df_unique$Markers <- gsub(' ','',df_unique$Markers)
  ####### FeaturePlot ######## 
  if (length(unique(df_unique$Markers)>30)){
    for (i in mks_ref$CellType){
      mks=df_long[df_long$CellType==i,]$Markers
      ct1=gsub('[ |+|/]','_',i)
      file1=paste0(outdir,'/FeaturePlot-',ct1)
      mks1 = intersect(unique(mks),rownames(data))
      if(length(mks1)>0){
      	FeatureDimPlot(data,features=mks1,theme_use='theme_blank',pt.size=0.01,title=i)
      	ggsave(paste0(file1,'.png'),width=12,height=12)
      	ggsave(paste0(file1,'.pdf'),width=12,height=12)
      }
    }
  }else{
    FeatureDimPlot(data,features = unique(df_unique$Markers),theme_use = 'theme_blank',pt.size=0.01)
    ggsave(file.path(outdir,"FeaturePlot.pdf"),width=20,height=18)
    ggsave(file.path(outdir,"FeaturePlot.png"),width=20,height=18)
  }
}



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

##################### 2_marker_expression ##########################
# 判断refmarker.file文件是否存在，不存在跳过
if(is.null(refmarker.file)){refmarker.file='NA'}
if(file.exists(refmarker.file)){
  print(paste0(Sys.time(),' - ############## start 2_marker_expression'))
  if(!dir.exists(paste0(outdir,'/2_marker_expression'))){dir.create(paste0(outdir,'/2_marker_expression'))}
  plot_refmarkers(data,refmarker.file,outdir=paste0(outdir,'/2_marker_expression'))
  files <- list.files(paste0(outdir,'/2_marker_expression'),pattern='2_FeaturePlot')
  example.file <- files[grep('png',files)][1]
  file.copy(paste0(outdir,'/2_marker_expression/',example.file),paste0(outdir,'/2_marker_expression/2_FeaturePlot.example.png'),)
  print(paste0(Sys.time(),' - ############## 2_marker_expression done!'))
}else{
  print(paste0(refmarker.file,' file not exists!'))
}






