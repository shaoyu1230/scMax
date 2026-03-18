library(tidyr)  
library(dplyr)
library(ggplot2)
library(SCP)
library(Seurat)

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
#  mks=unique(df_unique$Markers)
#  FeatureDimPlot(data.f,features = mks,theme_use = 'theme_blank',ncol=6,pt.size=0.01)
#  ht <- ceiling(length(mks)/6)*3
#  ggsave(paste0(outdir,"/6_FeaturePlot.pdf"),width=18,height=ht,limitsize=FALSE)
#  ggsave(paste0(outdir,"/6_FeaturePlot.png"),width=18,height=ht,limitsize=FALSE)
}
## 细胞类型差异分析
## Usage
## celltype_de(data.f,rdsdir,outdir='1_cluster_properties/4_')
celltype_de <- function(data,degs,rdsdir,celltype_col='CellType',outdir='1_cluster_properties/4_Cluster.DE'){
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

