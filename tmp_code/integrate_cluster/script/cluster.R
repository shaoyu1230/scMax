library('getopt')
para<- matrix(c(
        'help',         'h',    0,      "logical",
        'rds',          'r',    1,      "character",
        'config',       'c',    1,      "character",
        'outdir',       'o',    1,      "character",
		'meta',       'm',    2,      "character"
),byrow=TRUE,ncol=4)
# 2是表示可选的。
#===========================================================
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
        cat(getopt(para,usage=TRUE))
        cat("
        ==================================================================================================================
		按照不同的分辨率分群并绘制umap,细胞占比图,top marker基因图
        ==================================================================================================================
        Usage example:
        Rscript this.r -r rds -c config.ini  -o outdir -m meta.datafile
        Options:
        --help          h       NULL            get this help
        --rds           r       character       rds file for analysis by seurat [forced]
        --config        c       character       ini file for analysis by seurat [forced]
        --outdir        o       character       The resurt of out dir for analysis [forced]
		--meta        m       character       meta.data file for rdsfile [option]
        \n")
        q(status=1)
}
#===========================================================
getAbsolutePath <- function(relative_path) {
  return(normalizePath(relative_path))
}

if ( !is.null(opt$help) )       { print_usage(para) }
if ( is.null(opt$rds) )         { cat("Please input the rds data ...\n\n") ; print_usage(para)}
if ( is.null(opt$config) )      { cat("Please give the congfig file for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$outdir) )      { cat("Please give the outdir for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$meta) )      { opt$meta <- 'None'}


 # if (!is.null(opt$rename) & file.exists(opt$rename)){
   # rename <- getAbsolutePath(opt$rename)
 # }else{
   # rename <- NULL}


# 读取配置文件
library(configr)
rdsfile <-opt$rds
outdir<-opt$outdir
metafile<-opt$meta


ini.list <- read.config(file = opt$config)
#参数
resolution     			 <-ini.list$cluster_ana$resolution
resolutionlist<-unique(unlist(strsplit(resolution,split = ",",fixed=T)))

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
#suppressMessages(library(openxlsx))
suppressMessages(library(getopt))
library(harmony)
library(pheatmap)
library(cowplot) #plot_grid
library(SeuratWrappers) #RunPrestoAll

marker_gene <- function(immune.combined,min.pct = 0.1, logfc.threshold = 0.25,outdir=getwd(),pref='10x',test.use= "wilcox"){
        DefaultAssay(immune.combined) <- "RNA"
        cluster_num<-sort(unique(Idents(immune.combined)))
        print("Starting FindMarkers")
        #all.markers <- FindAllMarkers(immune.combined, only.pos = FALSE, min.pct =min.pct, logfc.threshold = logfc.threshold,test.use=test.use)
        all.markers <- RunPrestoAll(immune.combined,only.pos = FALSE, min.pct =min.pct, logfc.threshold = logfc.threshold,test.use=test.use)
        head(all.markers, n = 2)
        #marker_result <- subset(all.markers, select=c(7,1,2,3,4,5,6))
        marker_result <- subset(all.markers, select=c(8,1,3,4,5,6,7))
        names(marker_result)[names(marker_result)=='gene'] <- 'gene_name'
        #write.csv(all.markers,paste(outdir,'3_marker',paste(pref,'all.markers.csv',sep='_'),sep='/'),quote=F)
        write.csv(marker_result,paste(outdir,'3_marker',paste(pref,'all.markers.csv',sep='_'),sep='/'),quote=F,row.names=F)
        print("Finished FindMarkers")
        return (all.markers)
}


mkdirs <- function(outdir,fp) {
        if(!file.exists(file.path(outdir,fp))) {
#               mkdirs(dirname(fp))
                dir.create(file.path(outdir,fp))
        }else{
                        print(paste(fp,"Dir already exists!",sep="     "))
                        unlink(file.path(outdir,fp), recursive=TRUE)
                        dir.create(file.path(outdir,fp))
                }
}

join_c <- function(ts,Connector=''){
        paste(ts,collapse=Connector)
}

Parse_abspath_c <- function(input_abspath){
  tmp <- c( dirname(input_abspath), basename(input_abspath))
  return(tmp)
}

#读取rds文件
print("###开始读取rds文件........................................")
print(Sys.time())
print("")
rds0<-readRDS(rdsfile)
rds<-rds0

if (metafile!='None'){
meta<-read.csv(metafile,header=T,row.names=1)
rds@meta.data<-meta[rownames(rds0@meta.data),]
}


tmp<-rds
#根据指定的分辨率分群
res_list<-c()
for (res in resolutionlist){
	res<-as.numeric(res)
	tmp <- FindClusters(tmp,resolution =res)
	current.cluster.ids <- levels(Idents(tmp))
	#new.cluster.ids <- as.numeric(current.cluster.ids)+1
	new.cluster.ids <- as.numeric(current.cluster.ids)
	Idents(tmp) <- plyr::mapvalues(x = Idents(tmp), from = current.cluster.ids, to = new.cluster.ids)
	tmp@meta.data$seurat_clusters<-Idents(tmp)[rownames(tmp@meta.data)]
	table(tmp@meta.data$seurat_clusters)
	tmp[[paste0("res", res)]] <- tmp$seurat_clusters
	res_list<-c(res_list,paste0("res", res))
}


##------------------------
#开始绘图，需要分开，因为DefaultAssay(tmp) <- "RNA"，会影响FindClusters
DefaultAssay(tmp) <- "RNA"
#存数据
tmp0<-tmp
meta<-tmp0@meta.data[,setdiff(colnames(tmp0@meta.data),res_list)]
tmp0@meta.data<-meta
saveRDS(tmp0,file=paste0(outdir,'/../../final_obj.rds'))
write.csv(tmp0@meta.data,file=paste0(outdir,'/../../final_metadata.csv'),quote=F,row.names=T)

##############--------------------------
for (res in resolutionlist){
	res_name<-paste0('res',res)
	tmp$seurat_clusters  <-tmp[[res_name]]
	mkdirs(outdir,res_name)
	outdir_pre<-paste0(outdir,'/',res_name)
	#setwd(outdir_pre)
	##===========================================================
	#1_QC #绘制质控图
	#percent.mt','percent.HB','nCount_RNA','nFeature_RNA',"S.Score","G2M.Score"
	#qc_list<-c('percent.mt','percent.HB','nCount_RNA','nFeature_RNA',"S.Score","G2M.Score")
	#qc_list<-unique(unlist(strsplit(ini.list$Subclusters$qc_list,split = ",",fixed=T)))
	#qc_list<-c('percent.mt','percent.HB','nCount_RNA','nFeature_RNA')
	#VlnPlot
	celltype_colors <- palette_scp(unique(tmp@meta.data$seurat_clusters),palette='Paired',NA_keep=TRUE)
	qc_list0<-c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb")
	qc_list<-intersect(qc_list0,colnames(tmp@meta.data))
	p<-FeatureStatPlot(tmp, stat.by =qc_list, group.by = 'seurat_clusters',ncol=1,stack = TRUE,legend.position='none',ylab='',palcolor=celltype_colors,add_box=TRUE)
	pdf(paste0(outdir_pre,'/','1_qc_clusters_VlnPlot.pdf'),w=5+0.2*length(unique(tmp$seurat_clusters)),h=3+0.3*length(qc_list))
	print(p)
	dev.off()
	##===========================================================	
	#2_clusters #绘制分群结果
	#umap
	p1<-CellDimPlot(srt = tmp, group.by = 'seurat_clusters', label = TRUE,label_insitu=T,label.fg = "black",label.bg = NA,label.size=4,label_repel=F,pt.size=0.01,theme_use="theme_blank",palcolor=celltype_colors)
	p2<-CellDimPlot(srt = tmp, group.by = 'Sample', label = F,label_insitu=T,label.fg = "black",label.bg = NA,label.size=4,label_repel=F,pt.size=0.01,theme_use="theme_blank",palcolor=celltype_colors)
	p<-plot_grid(p1,p2,ncol = 2)
	pdf(paste0(outdir_pre,'/','2_clusters_sample_UMAP.pdf'),w=15,h=5)
	print(p)
	dev.off()
	pdf(paste0(outdir_pre,'/','2_clusters_UMAP.pdf'),w=7,h=5)
	print(p1)
	dev.off()
	#tsne	
	p1<-CellDimPlot(srt = tmp, reduction='tsne', group.by = 'seurat_clusters', label = TRUE,label_insitu=T,label.fg = "black",label.bg = NA,label.size=4,label_repel=F,pt.size=0.01,theme_use="theme_blank",palcolor=celltype_colors)
	p2<-CellDimPlot(srt = tmp, reduction='tsne', group.by = 'Sample', label = F,label_insitu=T,label.fg = "black",label.bg = NA,label.size=4,label_repel=F,pt.size=0.01,theme_use="theme_blank",palcolor=celltype_colors)
	p<-plot_grid(p1,p2,ncol = 2)
	pdf(paste0(outdir_pre,'/','2_clusters_sample_TSNE.pdf'),w=15,h=5)
	print(p)
	dev.off()
	pdf(paste0(outdir_pre,'/','2_clusters_TSNE.pdf'),w=7,h=5)
	print(p1)
	dev.off()	
	#按照样本分隔绘图
	n1<-length(unique(tmp$Sample))
	if (n1<=4){a=n1;b=1}else{a=4;b=ceiling(n1/4)}
	#umap
	CellDimPlot(srt = tmp, reduction='umap', group.by = 'seurat_clusters',split.by='Sample', label = TRUE,label_insitu=T,label.fg = "black",label.bg = NA,label.size=4,label_repel=F,pt.size=0.01,palcolor=celltype_colors,theme_use="theme_blank",ncol=a)
	ggsave(paste0(outdir_pre,'/','2_clusters_splitbySample_UMAP.pdf'),width=6*a,height=5*b)
	#tsne
	CellDimPlot(srt = tmp, reduction='tsne', group.by = 'seurat_clusters',split.by='Sample', label = TRUE,label_insitu=T,label.fg = "black",label.bg = NA,label.size=4,label_repel=F,pt.size=0.01,palcolor=celltype_colors,theme_use="theme_blank",ncol=a)
	ggsave(paste0(outdir_pre,'/','2_clusters_splitbySample_TSNE.pdf'),width=6*a,height=5*b)
	#按照组分隔绘图
	n1<-length(unique(tmp$Group))
	if (n1<=4){a=n1;b=1}else{a=4;b=ceiling(n1/4)}
	#umap
	CellDimPlot(srt = tmp, reduction='umap', group.by = 'seurat_clusters',split.by='Group', label = TRUE,label_insitu=T,label.fg = "black",label.bg = NA,label.size=4,label_repel=F,pt.size=0.01,palcolor=celltype_colors,theme_use="theme_blank",ncol=a)
	ggsave(paste0(outdir_pre,'/','2_clusters_splitbyGroup_UMAP.pdf'),width=6*a,height=5*b)
	#tsne
	CellDimPlot(srt = tmp, reduction='tsne', group.by = 'seurat_clusters',split.by='Group', label = TRUE,label_insitu=T,label.fg = "black",label.bg = NA,label.size=4,label_repel=F,pt.size=0.01,palcolor=celltype_colors,theme_use="theme_blank",ncol=a)
	ggsave(paste0(outdir_pre,'/','2_clusters_splitbyGroup_TSNE.pdf'),width=6*a,height=5*b)	
	#细胞占比分析
	#按照样本绘图
	a<-table(tmp$Sample,tmp$seurat_clusters)
	n<-length(unique(tmp$Sample))
	write.csv(a,file=paste0(outdir_pre,'/','2_cluster_sample_cellnum.csv'),quote=F)
	a<-data.frame(seurat_clusters=tmp@meta.data[,'seurat_clusters'],Samples=tmp@meta.data[,'Sample'])
	p0<-theme(legend.position='right',panel.grid=element_blank(), legend.background = element_rect(colour = NA),
			legend.title = element_blank(),legend.text = element_text(face="plain", color="black",size=8),
			axis.text.x = element_text(color="black",size=10,angle=30, hjust = 1),
			axis.text.y = element_text(color="black",size=10),
			axis.title.x = element_text(face="plain", color="black",size=12),
			axis.title.y = element_text(face="plain", color="black",size=12))
	pdf(paste0(outdir_pre,'/','2_clusters_sample_cellratio.pdf'),w=5+0.3*n,h=5)
	p<-ggplot(a, aes(Samples)) + geom_bar(aes(fill=seurat_clusters), position='fill',width=0.6)+labs(x=" ", y = "Cell type distribution",fill= "seurat_clusters")+theme_bw()+p0+scale_fill_manual(values=celltype_colors)+ guides(fill = guide_legend(ncol = 2))
	print(p)
	dev.off()
	#按照组绘图
	a<-table(tmp$Group,tmp$seurat_clusters)
	n<-length(unique(tmp$Group))
	write.csv(a,file=paste0(outdir_pre,'/','2_cluster_group_cellnum.csv'),quote=F)
	a<-data.frame(seurat_clusters=tmp@meta.data[,'seurat_clusters'],Samples=tmp@meta.data[,'Group'])
	p0<-theme(legend.position='right',panel.grid=element_blank(), legend.background = element_rect(colour = NA),
			legend.title = element_blank(),legend.text = element_text(face="plain", color="black",size=8),
			axis.text.x = element_text(color="black",size=10,angle=30, hjust = 1),
			axis.text.y = element_text(color="black",size=10),
			axis.title.x = element_text(face="plain", color="black",size=12),
			axis.title.y = element_text(face="plain", color="black",size=12))
	pdf(paste0(outdir_pre,'/','2_clusters_group_cellratio.pdf'),w=5+0.3*n,h=5)
	p<-ggplot(a, aes(Samples)) + geom_bar(aes(fill=seurat_clusters), position='fill',width=0.6)+labs(x=" ", y = "Cell type distribution",fill= "seurat_clusters")+theme_bw()+p0+scale_fill_manual(values=celltype_colors)+ guides(fill = guide_legend(ncol = 2))
	print(p)
	dev.off()
	##===========================================================	
	#3_markers #绘制分群结果
	##===========================================================
	#marker基因
	Idents(tmp)<-'seurat_clusters'
	cluster.averages <- AverageExpression(object = tmp, assays ='RNA',return.seurat = F)
	write.csv(cluster.averages$RNA,file=paste0(outdir_pre,'/','3_clusters_allgenes_averageExpression.csv'),quote=F,row.names=T)
	average.genes<-cluster.averages$RNA
	#all.markers<-FindAllMarkers(tmp, only.pos = T, min.pct =0.1, logfc.threshold = 0.25,test.use='wilcox')
	# all.markers<-marker_gene(immune.combined,min.pct =as.numeric(ini.list$Para$marker_gene_min.pct), logfc.threshold = as.numeric(ini.list$Para$marker_gene_logfc.threshold),outdir=outdir,test.use=ini.list$Para$marker_gene_test.use,pref=prefix)
	#SeuratWrappers::RunPrestoAll 是 Seurat::FindAllMarkers 基于 Presto 的实现，功能参数和输出与 Seurat::FindAllMarkers 基本一致，由于 Presto 底层是基于 Rcpp 重编译的，可以非常迅速的进行 Wilcoxon 秩和检验。
	#remotes::install_github('satijalab/seurat-wrappers')
	#/annogene/data2/bioinfo/SCV/Repositories/Pipeline/Stable/SinCell_10X/SinCell_10X_v4.9.5/bin/Integration/script
	all.markers <- RunPrestoAll(tmp,only.pos = T, min.pct =0.1, logfc.threshold = 0.25,test.use='wilcox')
	marker_result <- subset(all.markers, select=c(8,1,3,4,5,6,7))
	markers<-subset(marker_result,p_val_adj<0.05)
	da<-arrange(markers,cluster,desc(avg_log2FC))
	write.csv(da,file=paste0(outdir_pre,'/','3_clusters_all.markers_FDR0.05.csv'),quote=F,row.names=T)
	#da1<-da[-grep("^mt-",da$gene),]
	topn<-10
	top_gene <- da %>% group_by(cluster) %>% top_n(topn,avg_log2FC)
	filename<-paste0(outdir_pre,'/','3_cluster_top10_markers_genes_pheatmap.pdf')
	genes<-unique(as.vector(top_gene$gene))
	data<-average.genes
	da0<-log2(data[genes,]+1)
	n1<-length(unique(tmp$seurat_clusters))
	colorlist<-c("navy","white","firebrick3")
	breaksList = seq(-3, 3, by = 0.01)
	color = colorRampPalette(colorlist)(length(breaksList))
	pheatmap::pheatmap(da0,filename=filename,width=8+0.3*n1,height=10+0.5*n1,cluster_rows=F,cluster_cols=F,display_numbers = F,number_format = "%.0f",fontsize = 8,fontsize_col = 10,border_color=NA,angle_col = "45",scale='row',color=color)

	#top10
	da %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top5
	genes<-unique(as.vector(top5$gene))
	#options (repr.plot.width=28, repr.plot.height=5)
	p0=DotPlot(tmp,features =genes,group.by = "seurat_clusters")#,'FTH1','FTL'
	p3<-p0+theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+xlab('')+ylab('')+
	  guides(color = guide_colorbar(title = 'Scaled expression',order = 1),size = guide_legend("Percent expressed"),order = 0)+
	  scale_colour_gradientn(colours = c("dodgerblue1","#44c1f0", "lightgoldenrod","#e20612",'#cc340c'))+
	  theme(legend.position = "top")+
	  theme(panel.border = element_rect(fill=NA,color="black", size=0.3, linetype="solid"))+
	  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
	pdf(paste0(outdir_pre,'/','3_cluster_top10_markers_genes_DotPlot.pdf'), w=10+n1,h=3+0.2*n1)
	print(p3)
	dev.off()
	#绘制top5 FeaturePlot
	# mkdirs(paste0(outdir_pre,'/3_markers/'),'FeaturePlot_top5')
	# da %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
	# genes<-unique(as.vector(top5$gene))
	# for (g in genes){
	# p1<-FeaturePlot(tmp,features=g,cols = c("lightgrey", "red"),label=F,pt.size=0,raster=F,label.size=3)
	# pdf(paste0('FeaturePlot_top5/',g,'_FeaturePlot.pdf'),w=5,h=4)
	# print(p1)
	# dev.off()
	# }
	##===========================================================
	print("###........................................")
	print("完成所有分析......")
	print(Sys.time())
}
