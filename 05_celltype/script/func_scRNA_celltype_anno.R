library(tidyr)  
library(dplyr)

get_frac <- function(meta,fill.bar,x){
  plotdatalist <- list()
  for (i in x) {
    plotdata <- meta %>% group_by_(i, fill.bar) %>% dplyr::summarise(count = n()) %>%
      dplyr::mutate(count/sum(count)) %>% dplyr::arrange_(fill.bar)
    plotdata <- as.data.frame(plotdata)
    plotdata$facet <- i
    colnames(plotdata) <- c(x, fill.bar, "count", "frac","facet")
    plotdatalist[[i]] <- as.data.frame(plotdata)
  }
  p_multi <- do.call(rbind, plotdatalist)
  outfile <- p_multi[,1:4]
  return(outfile)
}

## Usage：
# cluster_plots(data,outputDir,groups=c('seurat_clusters','Sample','Group'))
cluster_plots <- function(data,outputDir,groups=c('Sample','Group'),palcolor=NULL){
  ############ 1_QC ###############
  print(paste0(date(),"- QC plots started."))
  FeatureStatPlot(data,stat.by=c("nCount_RNA","nFeature_RNA","percent.mt"),group.by='seurat_clusters',add_box=TRUE,stack=TRUE,palcolor=palcolor)
  ggsave(paste0(outputDir,'/0_QC.pdf'),width=10,height=7)
  ggsave(paste0(outputDir,'/0_QC.png'),width=10,height=7)
  print(paste0(date(),"- QC plots finished"))
  
  ############ 2_UMAP plots ###########
  print(paste0(date(),"- UMAP plots started."))
  celltype_umap_plots(data,outputDir,celltype_col = 'seurat_clusters',groups = groups,palcolor=palcolor)
  print(paste0(date(),"- UMAP plots finished"))
  #################### 3_Fraction plots ###################
  print(paste0(date(),"- Fraction plots started."))
  for (groupby in groups[groups!='seurat_clusters']){
    # get frac
    df=get_frac(data@meta.data,'seurat_clusters',groupby)
    write.csv(df,paste0(outputDir,'/3_Frac-seurat_clusters.in.',groupby,'.csv'),row.names=FALSE)
    # stat heatmap
    wd <- max(ceiling(length(unique(df$seurat_clusters))/2),6)
    ht <- max(ceiling(length(unique(df[,groupby]))/2),4)
    plot_stat_heatmap(df,prefix=paste0(outputDir,'/3_Frac-seurat_clusters.in.',groupby),x='seurat_clusters',y=groupby,fill='frac',width = wd,height = ht)
    plot_stat_heatmap(df,prefix=paste0(outputDir,'/3_Count-seurat_clusters.in.',groupby),x='seurat_clusters',y=groupby,fill='count',width = wd,height = ht)
    # frac barplot, trend
    plot_stat_bar(data,statby='seurat_clusters',groupby=groupby,prefix=paste0(outputDir,'/3_Frac-'),palcolor=palcolor)
  }
  print(paste0(date(),"- Fraction plots finished"))
}

plot_stat_bar <- function(data,statby='seurat_clusters',groupby='Sample',prefix='./',palcolor=NULL){
  statby1 <- gsub('[ |+]','.',statby)
  # frac barplot
  n=length(unique(data@meta.data[,groupby]))
  wd=max(5,ceiling(n/3))
  print(paste0('width:',wd))
  CellStatPlot(data, stat.by = statby, group.by = groupby, label = FALSE,palcolor=palcolor)
  ggsave(paste0(prefix,statby1,".in.",groupby,"-bar.pdf"),width =wd, height = 4)
  ggsave(paste0(prefix,statby1,".in.",groupby,"-bar.png"),width =wd, height = 4)
  # frac trend
  CellStatPlot(data, stat.by = statby, group.by = groupby, plot_type='trend', label = FALSE,palcolor=palcolor)
  ggsave(paste0(prefix,statby1,".in.",groupby,"-trend.pdf"),width =wd, height = 4)
  ggsave(paste0(prefix,statby1,".in.",groupby,"-trend.png"),width = wd, height = 4)
}
#### cell count/frac heatmap ####
plot_stat_heatmap <- function(data, prefix = "./CellType.percentage.heatmap",
                              x='Sample',y='CellType',fill=c('count','frac')[1],
                                    width = 6, height = 4, 
                                    low_color = "white", high_color = "red",
                                    label_format = "%.2f",font_size = 12) {
  if (fill == "count") {
    data$fill <- data$count
    label_format <- "%d"
    title <- "Cell count"
  }else if (fill == "frac") {
    data$fill <- data$frac * 100
    label_format <- "%.2f"
    title <- "Cell fraction (%)"
  }
  data$x=data[,x];data$y=data[,y]
  # 生成热图
  p <- ggplot(data,aes(x = x, y = y, fill = fill)) + 
    geom_tile(color = "white", linewidth = 0.5) + #添加热图
    geom_text(aes(label = sprintf(label_format, fill)),color = "black",size = 3) + #添加文本
    scale_fill_gradient(low = low_color,high = high_color,name ="") +
    theme_minimal(base_size = font_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_blank(),panel.grid = element_blank()) +
    labs(title = title)
  # 保存文件
  ggsave(filename = paste0(prefix, "-heatmap.pdf"),plot = p,width = width,height = height)
  ggsave(filename = paste0(prefix, "-heatmap.png"),plot = p,width = width,height = height)
  # 返回绘图对象（可选）
  invisible(p)
}

######### Roe plots #########
distribution_Roe = function(
    meta_data=metainfo,
    celltype_column = "majorCluster",
    celltype_level = NULL,
    condition_column = "tissue",
    condition_level = NULL,
    max_threshold = 2,min_threshold = 0.01,
    add_label = "number",#可选NA, "number", "sign"
    condition_label_angle = 60,condition_label_hjust = 1,
    celltype_color = NULL,relative_width = NULL,#这两个参数是共存的；都为空或都不为空。
    out_prefix='./', tile_color = NA,
    tile_fill = c("#f6f8e6","#eb632e")#"viridis"等颜色标题或者c("#f6f8e6","#ec6725")
){
  library(tidyverse)
  colnames(meta_data)[which(colnames(meta_data) == celltype_column)] = "celltypE"
  colnames(meta_data)[which(colnames(meta_data) == condition_column)] = "conditioN"
  if(is.null(celltype_level)){
    meta_data$celltypE = as.character(meta_data$celltypE)
    meta_data$celltypE = factor(meta_data$celltypE,levels = sort(unique(meta_data$celltypE)))
  } else {
    meta_data$celltypE = factor(meta_data$celltypE,levels = celltype_level)
  }
  
  if(is.null(condition_level)) {
    meta_data$conditioN = as.character(meta_data$conditioN)
    meta_data$conditioN = factor(meta_data$conditioN,levels = sort(unique(meta_data$conditioN)))
  } else {
    meta_data$conditioN = factor(meta_data$conditioN,levels = condition_level)
  }
  
  ### 卡方独立性检验
  #仅仅是用到了其中的chisq$expected，独立性检验的结论并不关心
  #Note that, Chi-square test should only be applied when the expected frequency of any cell is at least 5.
  contengency_table = xtabs(~celltypE+conditioN,data = meta_data)
  chisq <- chisq.test(as.matrix(contengency_table))
  Roe = chisq$observed / chisq$expected
  Roe = as.matrix(Roe) %>% as.data.frame()
  colnames(Roe)[3] = "value"
  write.csv(Roe,paste0(out_prefix,'Roe.csv'),quote=FALSE,row.names = FALSE)
  Roe$old_value = Roe$value
  Roe$value[Roe$value > max_threshold] = max_threshold
  Roe$value[Roe$value < min_threshold] = 0
  ### 画图
  Roe$sign = ""
  Roe$sign[Roe$value > 1] = "+++"
  Roe$sign[Roe$value > 0.8 & Roe$value <= 1] = "++"
  Roe$sign[Roe$value >= 0.2 & Roe$value <= 0.8] = "+"
  Roe$sign[Roe$value > 0 & Roe$value < 0.2] = "+/−"
  Roe$sign[Roe$value == 0] = "-"
  
  pb = Roe %>% ggplot(aes(x=conditioN,y=celltypE))+
    geom_tile(aes(fill = value),color = tile_color)+
    scale_y_discrete(expand = c(0,0),position = "right")+
    scale_x_discrete(expand = c(0,0))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y.right = element_text(size = 12,colour = "black"),
      axis.text.x.bottom = element_text(size = 12,colour = "black",angle = condition_label_angle,hjust = condition_label_hjust)
    )
  # 配色
  if (length(tile_fill) == 1) {
    pb=pb+scale_fill_viridis_c("Ro/e",option = tile_fill)
  } else if (length(tile_fill) == 2) {
    pb=pb+scale_fill_gradient("Ro/e",low = tile_fill[1],high = tile_fill[2])
  }
  # 添加符号
  if(is.na(add_label)) {
    pb=pb
  } else if (add_label == "number") {
    pb=pb+geom_text(aes(label=round(old_value,2)))
  } else if (add_label == "sign") {
    pb=pb+geom_text(aes(label=sign))
  }
  
  ### 返回值
  if (is.null(celltype_color) & is.null(relative_width)) {
    return(pb)
  } else if(!is.null(celltype_color) & !is.null(relative_width)) {
    library(patchwork)
    tmpdf = factor(levels(Roe$celltypE),levels = levels(Roe$celltypE)) %>% as.data.frame()
    colnames(tmpdf) = "celltype"
    
    pa = ggplot(data = tmpdf,aes(x=0,y=celltype))+
      geom_text(aes(label = celltype,color = celltype),hjust = 0,size = 5)+
      scale_color_manual(values = celltype_color)+
      scale_x_continuous(expand = c(0,0),limits = c(0,1))+
      theme_void()+theme(legend.position = "none")
    
    pb=pb+theme(legend.position = "top",axis.text.y.right = element_blank())
    pc = pb+pa+plot_layout(widths = c(1,relative_width))
    return(pc)
  } else {
    return(print("celltype_color and relative_width do not match!"))
  }
}

######### plot_refmarkers ############
# 功能：根据参考marker文件绘制DotPlot和FeaturePlot
# 参考marker文件格式：第一列为CellType，第二列为Markers，Markers之间用逗号或空格分隔，列名用Tab分隔
plot_refmarkers <- function(data,refmarker.file,outdir=paste0(outdir,'/marker_expression')){
  if (!dir.exists(outdir)){dir.create(outdir)}
  ####### read markers ########
  mks_ref <- read.csv(refmarker.file,sep='\t')
  library(tidyr)  
  library(dplyr) 
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
  ggsave(file.path(outdir,"1_DotPlot-RefMarkers.pdf"),width=wd,height=ht)
  ggsave(file.path(outdir,"1_DotPlot-RefMarkers.png"),width=wd,height=ht)
  ####### 2_FeaturePlot ######## 
  if (length(unique(df_unique$Markers)>30)){
    for (i in mks_ref$CellType){
      mks=df_long[df_long$CellType==i,]$Markers
      ct1=gsub('[ |+|/]','_',i)
      file1=paste0(outdir,'/2_FeaturePlot-',ct1)
      mks1 = intersect(unique(mks),rownames(data))
      if(length(mks1)>0){
      	FeatureDimPlot(data,features=mks1,theme_use='theme_blank',pt.size=0.01,title=i)
      	ggsave(paste0(file1,'.png'),width=12,height=12)
      	ggsave(paste0(file1,'.pdf'),width=12,height=12)
      }
    }
  }else{
    FeatureDimPlot(data,features = unique(df_unique$Markers),theme_use = 'theme_blank',pt.size=0.01)
    ggsave(file.path(outdir,"1_FeaturePlot-RefMarkers.pdf"),width=20,height=18)
    ggsave(file.path(outdir,"1_FeaturePlot-RefMarkers.png"),width=20,height=18)
  }
}

cluster_de <- function(data,rdsdir,outdir='cluster_characterization/4_Cluster.DE'){
  if(!dir.exists(outdir)){dir.create(outdir)}
  if(!dir.exists(rdsdir)){dir.create(rdsdir)}
  tmpdir <- paste0(rdsdir,'/../tmp')
  if(!dir.exists(tmpdir)){dir.create(tmpdir)}
  degs <- RunPrestoAll(data,group.by='seurat_clusters',only.pos=T) #比FindAllMarkers快
  degs <- data.table::setorderv(as.data.frame(degs),c("cluster","avg_log2FC","pct.1","pct.2"),c(1,-1,1,-1))
  degs$UpDown <- 'not.sig'
  degs[degs$p_val_adj<0.05 & degs$avg_log2FC>=0.25 & degs$pct.1>=0.1,]$UpDown <- 'Up'
  # degs[degs$p_val_adj<0.05 & degs$avg_log2FC <= -0.25 & degs$pct.1>=0.1,]$UpDown <- 'Down'
  saveRDS(degs,file.path(rdsdir,'ClusterDEGs.rds'))
  write.csv(degs,file.path(tmpdir,'ClusterDEGs.csv'),row.names=F)
  diffgene <- degs[degs$UpDown != 'not.sig',]
  diffgene %>% split(diffgene$cluster) %>% openxlsx::write.xlsx(file.path(outdir,'1_Cluster.DEGs.xlsx'))
  #### plots #####
  # top5 heatmap用于快速看cluster聚类
  ncluster <- length(unique(diffgene$cluster))
  topgene.df <- diffgene %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)
  ht <- GroupHeatmap(data,features = topgene.df$gene,group.by = "seurat_clusters",heatmap_palette = "RdBu",
                     show_row_names = TRUE, show_column_names = T,column_names_rot = 0,cluster_columns = T,
                     feature_split= topgene.df$cluster,assay = "RNA",anno_terms = FALSE,nlabel=0,
                     add_dot = TRUE,add_bg = T)
  ht$plot
  ggsave(paste0(tmpdir,'/ClusterTopDEGs.clustered.pdf'),width=20,height=25)
  ggsave(paste0(tmpdir,'/ClusterTopDEGs.clustered.png'),width=20,height=25)
  # top 10 DotPlot
  top10df <- degs %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
  DotPlot(data,features = unique(top10df$gene))+
    scale_colour_gradientn(colours = c("dodgerblue1","#44c1f0", "lightgoldenrod","#e20612",'#cc340c'))+
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2)+
    theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
  ggsave(file.path(outdir,'2_ClusterTop10DEGs.pdf'),width=28,height=6)
  ggsave(file.path(outdir,'2_ClusterTop10DEGs.png'),width=28,height=6)
  # top 10 FeaturePlot
  if(!dir.exists(file.path(outdir,'3_FeatureExprUMAP/'))){dir.create(file.path(outdir,'3_FeatureExprUMAP/'))}
  top10df <- degs %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
  # 绘制每个cluster高表达基因的UMAP表达图
  for (c in unique(top10df$cluster)){
    FeatureDimPlot(data,features=top10df[top10df$cluster==c,'gene'],pt.size=0.01,theme_use='theme_blank',ncol=5,title=paste0("Cluster ",c)) #2行5列
    ggsave(paste0(outdir,'/3_FeatureExprUMAP/Cluster_',c,'.png'),width=15,height=7)
    ggsave(paste0(outdir,'/3_FeatureExprUMAP/Cluster_',c,'.pdf'),width=15,height=7)
  }
  return(diffgene)
}

#### cluster 富集分析 ####
## Usage
# cluster_enrich(derds=file.path(outdir,'/Rdata/ClusterDEGs.rds'),
#                orgdb='org.Mmu.eg.db',organism_kegg='mcc',
#                outdir,
#                enrichdir=paste0(outdir,'Enrich'))
# 需要在node026运行kegg分析（联网），已改为本地数据库20250801
cluster_enrich <- function(derds,orgdb='org.Mmu.eg.db',organism_kegg='mcc',
                           enrichdir='Enrich'){
  library(KEGG.db)
  if(!dir.exists(enrichdir)){dir.create(enrichdir)}
  degs <- readRDS(derds)
  clusters <- unique(degs$cluster)
  degs.up <- degs[degs$p_val_adj<0.05 & degs$avg_log2FC>=0.25 & degs$pct.1>=0.1,]
  tmpdir <- paste0(enrichdir,'/../../tmp')
  if(!dir.exists(tmpdir)){dir.create(tmpdir)}
  for (c in clusters){
    genes <- degs.up[degs.up$cluster==c,]$gene
    write.table(genes,paste0(tmpdir,'/',c,'.up.txt'),row.names=F,quote=F,col.names=F)
  }
  print(paste0(date(),' - Enrichment preparation finish!'))
  
  if(!dir.exists(enrichdir)){dir.create(enrichdir)}
  #source('/annogene/data2/bioinfo/PMO/shaoyusong/03.scripts/02.tidy/Enrich/func_GeneList2GOKEGG.r')
  file.list <- list.files(tmpdir,pattern='up.txt')
  print(paste0(date(),' - Enrichment analysis starts!'))
  for(file in file.list){
    prefix=gsub('.txt|.txt','',file)
    file=file.path(tmpdir,file)
    GeneList2GOKEGG(infile=file,orgdb=orgdb,organism_kegg=organism_kegg,
                    outdir=enrichdir,
                    prefix=paste0(prefix,'.'),runGO=TRUE,runKEGG=TRUE,title=paste0(prefix,' '))
  }
  print(paste0(date(),' - Enrichment analysis finish!'))
}
#### celltype annotation ####
## 规则：样本分组列名为Group，细胞注释列名为CellType
#celltype_plots(data,'3_cluster_annotation/CellType_Keep/',groups=groupcols)
celltype_umap_plots <- function(data,outdir,celltype_col="CellType",groups=c('Sample','Group'),palcolor=NULL){
  if(!dir.exists(outdir)){dir.create(outdir,recursive = TRUE)}
  if(celltype_col!='seurat_clusters'){
    data$CellType <- data@meta.data[,celltype_col]
    celltype_col1 <- 'CellType'
  }else if(celltype_col=='seurat_clusters'){
    celltype_col1 <- 'seurat_clusters'
  }
  # celltype umap
  ncts <- length(unique(data@meta.data[,celltype_col]))
  wd=ifelse(ncts>14,8,6)
  CellDimPlot(srt = data, group.by = celltype_col, label = TRUE,label_insitu=TRUE,label.size=3,label_repel=TRUE,pt.size=0.01,palcolor=palcolor)
  ggsave(paste0(outdir,'/1_',celltype_col1,'.UMAP.pdf'),width=wd,height=5)
  ggsave(paste0(outdir,'/1_',celltype_col1,'.UMAP.png'),width=wd,height=5)
  CellDimPlot(srt = data, group.by = celltype_col, label = TRUE,label_insitu=TRUE,label.size=3,label_repel=TRUE,pt.size=0.01,theme_use='theme_blank',palcolor=palcolor)
  ggsave(paste0(outdir,'/1_',celltype_col1,'.UMAP-blank.pdf'),width=wd,height=5)
  ggsave(paste0(outdir,'/1_',celltype_col1,'.UMAP-blank.png'),width=wd,height=5)
  CellDimPlot(srt = data, group.by = celltype_col, split.by=celltype_col,label = FALSE,pt.size=0.01,theme_use='theme_blank',ncol=5,palcolor=palcolor)&NoLegend()
  ht <- ceiling(length(unique(data@meta.data[,celltype_col]))/5)*3
  ggsave(paste0(outdir,'/1_',celltype_col1,'.UMAP-split.pdf'),width=18,height = ht)
  ggsave(paste0(outdir,'/1_',celltype_col1,'.UMAP-split.png'),width=18,height = ht)
  # celltype umap-split by group
  for (groupby in groups){
    print(paste0("Plot group: ",groupby))
    if (class(data@meta.data[,groupby])=='integer'|class(data@meta.data[,groupby])=='numeric'){
      data@meta.data[,groupby] <- as.character(data@meta.data[,groupby])
    }
    CellDimPlot(srt = data, group.by = groupby, label = FALSE,pt.size=0.01)
    ggsave(paste0(outdir,'/1_',groupby,'.UMAP.pdf'),width=6,height=5)
    ggsave(paste0(outdir,'/1_',groupby,'.UMAP.png'),width=6,height=5)
    ngroup <- length(unique(data@meta.data[,groupby]))
    if (ngroup>=3){
      ncol=3
    }else if(ngroup==2){
      ncol=2
    }else{
      ncol=1
    }
    print(ncol)
    ncts <- length(unique(data@meta.data[,celltype_col]))
    h=ifelse(ncts>14,6,4)
#细胞种类多的话就加高umap图
    w=ifelse(ncts>14,8,5)
    CellDimPlot(srt = data, group.by = celltype_col, split.by=groupby,label = FALSE,pt.size=0.01,ncol=ncol, palcolor=palcolor)
    ht <- min(40,ceiling(length(unique(data@meta.data[,groupby]))/ncol)*h)
    wd <- min(40,ncol*w)
    ggsave(paste0(outdir,'/2_',groupby,'.UMAP-split.pdf'),
           width=wd,height = ht,limitsize = FALSE)
    ggsave(paste0(outdir,'/2_',groupby,'.UMAP-split.png'),
           width=wd,height=ht,limitsize = FALSE)
    
  }
}
###### fraction plots ######
# 绘制细胞比例柱状图，输出细胞数和比例统计表
cell_fraction_plots <- function(data,out_prefix="celltype_fraction/1_CellType_All.in.",celltype_col='CellType',
                                groups=c('Sample','Group'),do.boxplot=FALSE,palcolor=NULL){
  #################### fraction plots ###################
  print(paste0(date(),"- Fraction plots started."))
  # for (groupby in groups){
  #   out_prefix <- paste0(out_prefix,groupby)
  #   # barplot
  #   CellStatPlot(data, stat.by = celltype_col, group.by = groupby, label = FALSE)
  #   ggsave(paste0(out_prefix,"-bar.pdf"),width =5, height = 4)
  #   ggsave(paste0(out_prefix,"-bar.png"),width = 5, height = 4)
  #   # frac table
  #   df=get_frac(data@meta.data,celltype_col,groupby)
  #   write.csv(df,paste0(out_prefix,'.frac.csv'),row.names=FALSE)
  # }
  for (groupby in groups){
    # get frac
    df=get_frac(data@meta.data,celltype_col,groupby)
    write.csv(df,paste0(out_prefix,groupby,'.csv'),row.names=FALSE)
    # stat heatmap
    wd <- max(length(unique(df[,celltype_col]))/2,6)
    ht <- max(length(unique(df[,groupby]))/2,4)
    plot_stat_heatmap(df,prefix=paste0(out_prefix,groupby,'-frac'),x=celltype_col,y=groupby,fill='frac',width = wd,height = ht)
    plot_stat_heatmap(df,prefix=paste0(out_prefix,groupby,'-count'),x=celltype_col,y=groupby,fill='count',width = wd,height = ht)
    # frac barplot, trend
    CellStatPlot(data, stat.by = celltype_col, group.by = groupby, label = FALSE,palcolor=palcolor)
    wd <- max(length(unique(data@meta.data[,groupby]))/2,5)
    ht <- ifelse(length(unique(data@meta.data[,celltype_col]))>14,6,4)
    ggsave(paste0(out_prefix,groupby,"-bar.pdf"),width =wd, height = ht)
    ggsave(paste0(out_prefix,groupby,"-bar.png"),width = wd, height = ht)
    CellStatPlot(data, stat.by = celltype_col, group.by = groupby, label = FALSE,plot_type='trend',palcolor=palcolor)
    ggsave(paste0(out_prefix,groupby,"-trend.pdf"),width =wd, height = ht) 
    ggsave(paste0(out_prefix,groupby,"-trend.png"),width = wd, height = ht)
  }
  # boxplot
  # if(do.boxplot){
  #   df <- read.csv(paste0(outputDir,'/3_CellType.in.Group.frac.csv'))
  #   groups <- unique(df$Group)
  #   cmps <- combn(groups, 2, simplify = FALSE)
  #   plot_fracbox(df,celltype='CellType',x='Group',y='frac',color='Group',
  #                method='t .test',comparisons=cmps,
  #                col=NULL,
  #                outfile_prefix="/3_CellType.Fraction-box",
  #                wd=18,ht=9)
  # }
  print(paste0(date(),"- Fraction plots finished"))
}


plot_fracbox <- function(frac_df,celltype,x='Response.2',y='frac',color='Phase',
                         method='t.test',comparisons=list(c("PCR", "MPR"),
                                                          c("PCR", "PR"),
                                                          c("MPR", "PR")),
                         col=NULL,
                         outfile_prefix,wd=18,ht=9){
  plist <- list()
  frac_df$CellType <- frac_df[,celltype]
  celltypes <- unique(frac_df$CellType)
  if (is.null(col)){
    col <- palette_scp(unique(frac_df[,color]), palette = "Paired",  NA_keep = TRUE)
  }
  for(ct in celltypes){
    print(ct)
    df_tmp <- frac_df[frac_df$CellType==ct,]
    pb1 <- ggplot(df_tmp)+
      geom_boxplot(aes_string(x=x,y=y))+
      geom_point(size=3,aes_string(x=x,y=y,color=color))+scale_color_manual(values=col)+
      ggtitle(ct)+theme_classic()+
      theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
            axis.title.x=element_blank(),
            # legend.position = "NA",
            panel.background = element_blank())+
      stat_compare_means(method = method, hide.ns = FALSE, 
                         comparisons = comparisons,#label = "p.signif",
                         vjust=0.02,bracket.size=0.6) # 执行显著性测试，使用 t /wilcox检验，并添加不同组别之间的比较
    plist[[ct]] <- pb1   
  }
  #plot_grid(plotlist=plist, ncol=5)#只需要指定列数
  # ggarange可以合并图例，要指定几行几列
  ggarrange(plotlist=plist, ncol=5,nrow=3,legend = "right",common.legend = T)
  ggsave(paste0(outfile_prefix,'.pdf'),width=wd,height=ht)
  ggsave(paste0(outfile_prefix,'.png'),width=wd,height=ht)
}

## add_anno功能：添加细胞类型注释
add_anno <- function(data,outdir,annofile,celltype.levels,celltype_col="CellType"){
  #annofile: Cluster,CellType,Markers
  ##### Add annotation #######
  anno.df <- read.csv(annofile,sep='\t')
  # 检查anno.df的cluster数量是否和data一致
  if(length(unique(anno.df$Cluster)) != length(unique(data$seurat_clusters))){
    stop('anno.df的cluster数量和data不一致，请检查！')
  }
  anno.majortype <- anno.df$CellType
  names(anno.majortype) <- anno.df$Cluster
  Idents(data) <- data$seurat_clusters# 默认data的Idents为seurat_clusters列
  data <- RenameIdents(data, anno.majortype)
  data@meta.data[,celltype_col] <- Idents(data)
  if(!is.null(celltype.levels) & length(celltype.levels)>0){
    if(length(celltype.levels)==length(unique(data@meta.data[,celltype_col]))){
      data@meta.data[,celltype_col] <- factor(data@meta.data[,celltype_col],levels=c(celltype.levels) )
      Idents(data) <- data@meta.data[,celltype_col]
    }else{
      print('CellType levels not equal to Celltypes, CellTypes: ')
      print(unique(data@meta.data[,celltype_col]))
      print('CellType levels: ')
      print(celltype.levels)
      # 提示差异
      print('Diff celltypes: ', setdiff(unique(data@meta.data[,celltype_col]),celltype.levels))
      print('Diff celltype levels: ', setdiff(celltype.levels,unique(data@meta.data[,celltype_col])))
      stop('Please check the celltype levels and celltypes!')
    }
  }
  # saveRDS(data,paste0(outdir,'/Data-',rdsprefix,'_CellType.rds'))
  return(data)
}

anno_mks_plot1 <- function(data,annofile,celltype_col="CellType",celltype.levels,outdir){
  # feature dotplot
  library(tidyr)  
  library(dplyr)
  #### all cells ####
  anno.df <- read.csv(annofile,sep='\t')
  df_long <- anno.df[,c('CellType','Markers')] %>%  separate_rows(Markers, sep = "[, ]+")
  df_unique <- df_long %>% unique() %>%  data.frame()  
  if(is.null(celltype.levels)|length(celltype.levels)<length(unique(anno.df$CellType))){
    celltype.levels <- unique(df_unique$CellType)
  }
  df_unique$CellType<- factor(df_unique$CellType,levels=celltype.levels)
  df_unique<- df_unique[order(df_unique$CellType),]
  # celltype-cluster heatmap用于检查注释是否准确
  ht <- GroupHeatmap(data,features = df_unique$Markers,
                     feature_split= df_unique$CellType,
                     group.by = celltype_col,split.by="seurat_clusters",
                     heatmap_palette = "RdBu",
                     show_row_names = TRUE, show_column_names = T,column_names_rot = 30,cluster_columns = F,
                     assay = "RNA",anno_terms = FALSE,nlabel=0,
                     add_dot = TRUE,add_bg = TRUE)
  ht$plot
  wd <- max(9,length(unique(data$seurat_clusters))/1.4)
  ht <- max(length(unique(df_unique$Markers))*0.3,10)
  ggsave(paste0(outdir,"/3_CellType.dotplot.clusterAnno.pdf"),width=wd,height=ht)
  ggsave(paste0(outdir,"/3_CellType.dotplot.clusterAnno.png"),width=wd,height=ht)
  # dotplot, by seurat_clusters
  DotPlot(data,features =unique(df_unique$Markers),group.by="seurat_clusters")+ #coord_flip()+
    scale_colour_gradientn(colours = c("dodgerblue1","#44c1f0", "lightgoldenrod","#e20612",'#cc340c'))+
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2)+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.title=element_blank())
  wd <- max(ceiling(length(mks)/3),9)
  ht <- max(length(unique(data@meta.data[,celltype_col]))*0.4,5)
  ggsave(paste0(outdir,"/4_cluster.dotplot.png"),width=wd,height=ht)
  ggsave(paste0(outdir,"/4_cluster.dotplot.pdf"),width=wd,height=ht)
  # dotplot
  DotPlot(data,features =unique(df_unique$Markers),group.by=celltype_col)+ #coord_flip()+
    scale_colour_gradientn(colours = c("dodgerblue1","#44c1f0", "lightgoldenrod","#e20612",'#cc340c'))+
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2)+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.title=element_blank())
  mks <- unique(df_unique$Markers)
  wd <- max(ceiling(length(mks)/3),9)
  ht <- max(length(unique(data@meta.data[,celltype_col]))*0.4,5)
  ggsave(paste0(outdir,"/4_CellType.dotplot.png"),width=wd,height=ht)
  ggsave(paste0(outdir,"/4_CellType.dotplot.pdf"),width=wd,height=ht)
  # FeaturePlot
  for(ct in unique(df_unique$CellType)){
    mks=unique(df_unique[df_unique$CellType==ct,]$Markers)
    mks <- intersect(mks,rownames(data))
    if(length(mks)>=1){
      ncol=min(length(mks),5)
      FeatureDimPlot(data,features = mks,theme_use = 'theme_blank',ncol=ncol,pt.size=0.01)
      ht <- ceiling(length(mks)/ncol)*3
      wd <- ncol*3
      ggsave(paste0(outdir,"/6_FeaturePlot-",ct,".pdf"),width=wd,height=ht,limitsize=FALSE)
      ggsave(paste0(outdir,"/6_FeaturePlot-",ct,".png"),width=wd,height=ht,limitsize=FALSE)
    }
  }
  file.copy(paste0(outdir,"/6_FeaturePlot-",df_unique$CellType[1],".png"),paste0(outdir,"/6_FeaturePlot-example.png"))

#  mks=unique(intersect(rownames(data),df_unique$Markers))
#  ht <- ceiling(length(mks)/6)
#  FeatureDimPlot(data,features = mks,theme_use = 'theme_blank',ncol=6,pt.size=0.01)
#  ggsave(paste0(outdir,"/5_FeaturePlot.pdf"),width=18,height=ht*3)
#  ggsave(paste0(outdir,"/5_FeaturePlot.png"),width=18,height=ht*3)
}


#### CellType_Keep dotplot ####
anno_mks_plot2 <- function(data.f,annofile,celltype_col="CellType",celltype.keep,outdir){
  anno.df <- read.csv(annofile,sep='\t')
  # celltype heatmap
  anno.df1 <- unique(anno.df[anno.df$CellType %in% celltype.keep,c('CellType','Markers')])
  df_unique <- anno.df1 %>%  separate_rows(Markers, sep = "[, ]+") %>% unique() %>% data.frame()
  df_unique$CellType<- factor(df_unique$CellType,levels=celltype.keep)
  df_unique <- df_unique[order(df_unique$CellType),]
  ht <- GroupHeatmap(data.f,features = df_unique$Markers,#feature_split= df_unique$CellType,
                     row_names_side = 'right',
                     column_names_side = 'bottom',reticle_color = NA,
                     group.by = celltype_col,
                     heatmap_palette = "RdBu",
                     show_row_names = TRUE, show_column_names = T,
                     column_names_rot = 60,cluster_columns = F,
                     assay = "RNA",anno_terms = FALSE,nlabel=0,add_bg = FALSE)
  ht$plot
  wd <- max(9,length(unique(data@meta.data[,celltype_col]))*0.9)
  ht <- max(length(unique(df_unique$Markers))*0.3,8)
  ggsave(paste0(outdir,"/3_CellType.heatmap.png"),width=wd,height=ht)
  ggsave(paste0(outdir,"/3_CellType.heatmap.pdf"),width=wd,height=ht)
  DotPlot(data.f,features =unique(df_unique$Markers),group.by=celltype_col)+ #coord_flip()+
    scale_colour_gradientn(colours = c("dodgerblue1","#44c1f0", "lightgoldenrod","#e20612",'#cc340c'))+
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2)+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.title=element_blank())
  ht <- max(ceiling(length(unique(df_unique$CellType))/3)*0.4,5)
  wd <- max(ceiling(length(unique(df_unique$Markers))/3),10)
  ggsave(paste0(outdir,"/4_CellType.dotplot.png"),width=wd,height=ht)
  ggsave(paste0(outdir,"/4_CellType.dotplot.pdf"),width=wd,height=ht)
  ht <- max(ceiling(length(unique(df_unique$Markers))/3),10)
  FeatureStatPlot(data.f,stat.by=df_unique$Markers,group.by=celltype_col,stack=TRUE)
  ggsave(paste0(outdir,"/5_CellType.violin.png"),width=10,height=ht)
  ggsave(paste0(outdir,"/5_CellType.violin.pdf"),width=10,height=ht)
  for(ct in unique(df_unique$CellType)){
    mks=unique(df_unique[df_unique$CellType==ct,]$Markers)
    mks <- intersect(mks,rownames(data.f))
    if(length(mks)>=1){
      ncol=min(length(mks),5)
      FeatureDimPlot(data.f,features = mks,theme_use = 'theme_blank',ncol=ncol,pt.size=0.01)
      ht <- ceiling(length(mks)/ncol)*3
      wd <- ncol*3
      ggsave(paste0(outdir,"/6_FeaturePlot-",ct,".pdf"),width=wd,height=ht,limitsize=FALSE)
      ggsave(paste0(outdir,"/6_FeaturePlot-",ct,".png"),width=wd,height=ht,limitsize=FALSE)
    }
  }
  file.copy(paste0(outdir,"/6_FeaturePlot-",df_unique$CellType[1],".png"),paste0(outdir,"/6_FeaturePlot-example.png"))
}
## 细胞类型差异分析
## Usage
## celltype_de(data.f,rdsdir,outdir='1_cluster_properties/4_')
celltype_de <- function(data,degs,rdsdir,celltype_col='CellType',outdir='cluster_characterization/4_Cluster.DE'){
  if(!dir.exists(outdir)){dir.create(outdir)}
  if(!dir.exists(rdsdir)){dir.create(rdsdir)}
  tmpdir <- paste0(outdir,'/../../tmp')
  if(!dir.exists(tmpdir)){dir.create(tmpdir)}
  degs <- RunPrestoAll(data,group.by=celltype_col,only.pos=T)
  degs<-setorderv(as.data.frame(degs),c("cluster","avg_log2FC","pct.1","pct.2"),c(1,-1,1,-1))
  degs$UpDown <- 'not.sig'
  degs[degs$p_val_adj<0.05 & degs$avg_log2FC>=0.25 & degs$pct.1>=0.1,]$UpDown <- 'Up'
  # degs[degs$p_val_adj<0.05 & degs$avg_log2FC <= -0.25 & degs$pct.1>=0.1,]$UpDown <- 'Down'
  saveRDS(degs,paste0(rdsdir,'/CellTypeDEGs.rds'))
  write.csv(degs,paste0(rdsdir,'/../tmp/CellTypeDEGs.csv'),row.names=F)
  diffgene <- degs[degs$UpDown != 'not.sig',]
  diffgene %>% split(diffgene$cluster) %>% openxlsx::write.xlsx(paste0(outdir,'/1_CellType.DEGs.xlsx'))
  #### plots #####
  # top5 heatmap用于快速看cluster聚类
  ncluster <- length(unique(diffgene$cluster))
  topgene.df <- diffgene %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)
  ht <- GroupHeatmap(data,features = topgene.df$gene,group.by = celltype_col,heatmap_palette = "RdBu",
                     show_row_names = TRUE, show_column_names = T,column_names_rot = 45,cluster_columns = FALSE,
                     assay = "RNA",anno_terms = FALSE,nlabel=0,
                     add_dot = FALSE,add_bg = FALSE)
  ht$plot
  wd <- max(9,length(unique(data@meta.data[,celltype_col]))/1.9)
  ht <- max(length(unique(topgene.df$gene))*0.2,12)
  ggsave(paste0(outdir,'/2_CellType.Top5DEGs-heatmap.pdf'),width=wd,height=ht)
  ggsave(paste0(outdir,'/2_CellType.Top5DEGs-heatmap.png'),width=wd,height=ht)
  # top 10 DotPlot
  #if(ncluster>10){topn=5}
  topn=5
  top10df <- degs %>% group_by(cluster) %>% top_n(n=topn,wt=avg_log2FC)
  DotPlot(data,features = unique(top10df$gene))+
    scale_colour_gradientn(colours = c("dodgerblue1","#44c1f0", "lightgoldenrod","#e20612",'#cc340c'))+
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.2)+
    theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
  wd <- max(ceiling(length(unique(top10df$gene))/3),15)
  ht <- max(ceiling(length(unique(data@meta.data[,celltype_col]))*0.6),6)
  ggsave(paste0(outdir,'/3_CellType.Top5DEGs.pdf'),width=wd,height=ht)
  ggsave(paste0(outdir,'/3_CellType.Top5DEGs.png'),width=wd,height=ht)
  dev.off()
  # top 3 FeaturePlot
  # if(!dir.exists(file.path(outdir,'4_FeatureExprUMAP/'))){dir.create(file.path(outdir,'4_FeatureExprUMAP/'))}
  # top10df <- degs %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
  # 绘制每个cluster高表达基因的UMAP表达图
  top5.df <- diffgene %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)
  for (c in unique(diffgene$cluster)){
    FeatureDimPlot(data,features=top5.df[top5.df$cluster==c,'gene'],pt.size=0.01,theme_use='theme_blank',ncol=5,title=c) #1行5列
    c <- gsub('[ |+|-|/]','.',c)#去除文件名中的特殊字符
    ggsave(paste0(outdir,'/4_FeatureExprUMAP-',c,'.png'),width=15,height=7)
    ggsave(paste0(outdir,'/4_FeatureExprUMAP-',c,'.pdf'),width=15,height=7)
  }

#  top5.df <- diffgene %>% group_by(cluster) %>% top_n(n=3,wt=avg_log2FC)
#  FeatureDimPlot(data,features=unique(top5.df$gene),pt.size=0.01,theme_use='theme_blank',ncol=5) #5列
#  ht <- ceiling(length(unique(top5.df$gene))/5)*3
#  ggsave(paste0(outdir,'/4_FeatureExprUMAP.png'),width=15,height=ht,limitsize = FALSE)
#  ggsave(paste0(outdir,'/4_FeatureExprUMAP.pdf'),width=15,height=ht,limitsize = FALSE)
  
  return(diffgene)
}

# top100 DEGs enrichment
library(KEGG.db)
celltype_enrich <- function(derds='CellTypeDEGs.rds',orgdb='org.Hs.eg.db',organism_kegg='hsa',
                            enrichdir='celltype_characterization/2_CellType_Enrich'){
  if(!dir.exists(file.path(enrichdir,'../../Rdata/'))){dir.create(file.path(enrichdir,'../../Rdata/'),recursive = TRUE)}
  diffgene <- readRDS(derds)
  top100df <- diffgene[diffgene$UpDown=='Up',] %>% group_by(cluster) %>% top_n(n=100,wt=avg_log2FC)
  go <- compareCluster(gene~cluster,data=top100df,fun='enrichGO',
                       pvalueCutoff = 0.05,qvalueCutoff = 0.05,ont='ALL',
                       OrgDb=orgdb,keyType='SYMBOL')
  write.csv(go@compareClusterResult,file.path(enrichdir,'1_CellType.GO.csv'),quote=F,row.names=F)
  saveRDS(go,file.path(enrichdir,'../../Rdata/CellType.GO.rds'))
  go@compareClusterResult <- subset(go@compareClusterResult,ONTOLOGY=="BP")
  dotplot(go,showCategory=3,color="p.adjust",font.size=10,label_format=35,
          title='GO terms of top100 upregulated genes')+
    theme(panel.grid = element_blank(),axis.text.x=element_text(angle=30,hjust=1),
          axis.title.x=element_blank())+
    scale_fill_gradientn(colors=c('darkred',"red","pink"))
  cts.count <- length(unique(top100df$cluster))
  wd <- max(cts.count,8)
  ht <- max(cts.count,8)
  ggsave(file.path(enrichdir,'2_CellType.GO.pdf'),width=wd,height=ht)
  ggsave(file.path(enrichdir,'2_CellType.GO.png'),width=wd,height=ht)
  
  eg <- bitr(top100df$gene,fromType="SYMBOL", toType=c("ENTREZID"),OrgDb=orgdb)
  top100df=left_join(top100df,eg,by=c('gene'="SYMBOL"))
  top100df <- top100df[!is.na(top100df$ENTREZID),]
  kegg <- compareCluster(ENTREZID~cluster,data=top100df,fun = "enrichKEGG",  
                         keyType = 'kegg',  #KEGG 富集
			 use_internal_data = TRUE, #使用本地KEGG数据库，ssy于20250715更新
                         organism=organism_kegg,pvalueCutoff = 0.05,qvalueCutoff = 0.05 )
  write.csv(kegg@compareClusterResult,file.path(enrichdir,'3_CellType.KEGG.csv'),quote=F,row.names=F)
  saveRDS(kegg,file.path(enrichdir,'../../Rdata/CellType.KEGG.rds'))
  dotplot(kegg,showCategory=3,color="p.adjust",font.size=10,label_format=35,
          title='KEGG terms of top100 upregulated genes')+
    theme(panel.grid = element_blank(),axis.text.x=element_text(angle=30,hjust=1),
          axis.title.x=element_blank())+
    scale_fill_gradientn(colors=c('darkred',"red","pink"))
  ggsave(file.path(enrichdir,'4_CellType.KEGG.pdf'),width=wd,height=ht)
  ggsave(file.path(enrichdir,'4_CellType.KEGG.png'),width=wd,height=ht)
}



GeneList2GOKEGG <- function(infile,orgdb='org.Rn.eg.db',organism_kegg='rno',
                            outdir,prefix='',runGO=TRUE,runKEGG=TRUE,title=''){
  library(stringr)
  library(clusterProfiler)
  library(openxlsx)
  library(dplyr)
  library(ggplot2)
  #### readin genelist
  intersectGenes <- read.table(infile)
  intersectGenes <- intersectGenes[,1]
  #### run GO
  if (runGO){
    Markers_GO_enrich <- enrichGO(intersectGenes, OrgDb = orgdb, ont = "ALL",
                                  pAdjustMethod = 'BH',
                                  pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                                  keyType = 'SYMBOL',readable = TRUE)
    if (!is.null(Markers_GO_enrich) & nrow(Markers_GO_enrich@result)>0){
		write.xlsx(Markers_GO_enrich@result,paste0(outdir,'/',prefix,'GO.xlsx'))
		plotgo(Markers_GO_enrich,outdir=outdir,prefix=prefix,titlei=paste0(title,'GO Enrichment'))
	}
  }
  
  #### run KEGG
  if (runKEGG){
    library(KEGG.db)
    eg <- bitr(intersectGenes, fromType="SYMBOL", toType=c("ENTREZID"),							                OrgDb=orgdb)
    genelist <- eg$ENTREZID[!is.na(eg$ENTREZID)]
    keggresult <- enrichKEGG(gene = genelist,use_internal_data = TRUE,
                               organism = organism_kegg,
                             keyType = "kegg",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2)
    if (!is.null(keggresult)){
		keggresult <- setReadable(keggresult,orgdb,keyType="ENTREZID")
		write.xlsx(keggresult@result,paste0(outdir,'/',prefix,'KEGG.xlsx'))
		plotkegg(keggresultfile=paste0(outdir,'/',prefix,'KEGG.xlsx'),
             outdir=outdir,prefix=prefix,titlei=paste0(title,'KEGG Enrichment'))
	}
  }
}

# plot top10 GO terms of each Oncology
plotgo <- function(Markers_GO_enrich,outdir='./',prefix='',
                   titlei = 'Pro-Lac TimePoint Common DEGs GO Enrichment'){
  CPCOLS<- c("#6495ED", "#8FBC8F", "#F4A460")
  mtx <- Markers_GO_enrich@result
  mtx$plt<- -log10(mtx$pvalue)
  mtx %>% group_by(ONTOLOGY) %>% top_n(10,plt) -> dmtx
  dorder = factor(rev(as.integer(rownames(dmtx))),labels=rev(dmtx$Description))
  x=dmtx$plt
  y=dorder
  pbar<-ggplot(dmtx,aes(x=Description,y=plt,fill=ONTOLOGY))+
    geom_bar(stat="identity",position=position_dodge(0.7),width=0.5,aes(x=dorder))+
	scale_x_discrete(labels=function(x)str_wrap(x,width = 70))+
    coord_flip() +scale_y_log10(breaks = c(1,10,100,1000))+
    theme(panel.background = element_rect(fill = "transparent",colour = NA))+
    xlab("GO_Term") +ylab(expression(-"log"["10"]*"(PValue)"))+theme_bw()+
    scale_fill_manual(values = CPCOLS)+labs(title=titlei)
  pdot <- ggplot(dmtx,aes(x,y)) +
    geom_point(aes(size=Count,color=-1*log(pvalue),shape=ONTOLOGY,))+
    scale_color_gradient(low = "steelblue", high = "indianred")+
	scale_y_discrete(labels=function(y)str_wrap(y,width = 70))+
    labs(color=expression(-log[10](pvalue)),size="Count",
         x=expression(-"log"["10"]*"(PValue)"),y="Go_term", title=titlei)+theme_bw()
  print(pbar)
  ggsave(paste0(outdir,'/',prefix,'GO_bar.pdf'),width=8,height=8)
  ggsave(paste0(outdir,'/',prefix,'GO_bar.png'),width=8,height=8)
  print(pdot)
  ggsave(paste0(outdir,'/',prefix,'GO_dot.png'),width=8,height=8)
  ggsave(paste0(outdir,'/',prefix,'GO_dot.pdf'),width=8,height=8)
}

plotkegg <- function(keggresultfile,outdir='./',prefix='',
                     titlei = ''){
  # CPCOLS<- c("#6495ED", "#8FBC8F", "#F4A460")
  mtx <- read.xlsx(keggresultfile)
  mtx$plt<- -log10(mtx$pvalue)
  if(nrow(mtx)>=30){
    mtx %>% top_n(30,plt) -> dmtx
  }else{
    dmtx <- mtx
  }
  dorder = factor(rev(as.integer(rownames(dmtx))),labels=rev(dmtx$Description))
  x=dmtx$plt
  y=dorder
  pbar <- ggplot(dmtx,aes(x=Description,y=plt,fill=plt))+
    geom_bar(stat="identity",position=position_dodge(0.7),width=0.5,aes(x=dorder))+
	scale_x_discrete(labels=function(x)str_wrap(x,width = 60))+
    coord_flip()+
    theme(panel.background = element_rect(fill = "transparent",colour = NA))+
    xlab("Term") +ylab(expression(-"log"["10"]*"(PValue)"))+
    theme_bw()+
    scale_fill_gradientn(colours=c("steelblue","yellow","red"))+
    labs(title=titlei,fill=expression(-log[10](pvalue)))
  pdot <- ggplot(dmtx,aes(x,y)) +
    geom_point(aes(size=Count,color=-1*log(pvalue),y=dorder))+
    scale_color_gradient(low = "steelblue", high = "indianred")+
	scale_y_discrete(labels=function(y)str_wrap(y,width = 70))+
    labs(color=expression(-log[10](pvalue)),size="Count",
         x=expression(-"log"["10"]*"(PValue)"),y="KEGG_term", title=titlei)+
    theme_bw()
  print(pbar)
  ggsave(paste0(outdir,'/',prefix,'KEGG_bar.pdf'),width=8,height = 8)
  ggsave(paste0(outdir,'/',prefix,'KEGG_bar.png'),width=8,height = 8)
  print(pdot)
  ggsave(paste0(outdir,'/',prefix,'KEGG_dot.pdf'),width=8,height = 8)
  ggsave(paste0(outdir,'/',prefix,'KEGG_dot.png'),width=8,height = 8)
}
