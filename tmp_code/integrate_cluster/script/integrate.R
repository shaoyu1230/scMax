library('getopt')
para<- matrix(c(
        'help',         'h',    0,      "logical",
        'rds',          'r',    2,      "character",
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
		去批次的方法有：rpca,harmony,cca,none
		推荐rpca或者harmony; cca整合分析主要用于大群，none是指不进行整合分群，也就是不去批次效应亚群细分（有些整合效果不好也可以尝试，一般不推荐）
		#nfeatures<-2000;npc<-30;dims<-20 #默认一般参数。
        ==================================================================================================================
        Usage example:
        Rscript this.r -r rds -c config.ini  -o outdir -m meta.datafile
        Options:
        --help          h       NULL            get this help
        --rds           r       character       rds file for analysis by seurat [option]
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
if ( is.null(opt$rds) )         { opt$rds <- 'None'}
if ( is.null(opt$config) )      { cat("Please give the congfig file for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$outdir) )      { cat("Please give the outdir for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$meta) )      { opt$meta <- 'None'}


 # if (!is.null(opt$rename) & file.exists(opt$rename)){
   # rename <- getAbsolutePath(opt$rename)
 # }else{
   # rename <- NULL}


# 读取配置文件
library(configr)
outdir<-opt$outdir


ini.list <- read.config(file = opt$config)
#参数
nvariablefeatures     			 <-as.numeric(ini.list$integrate_ana$nvariablefeatures)
reduction_npc_num     			 <-as.numeric(ini.list$integrate_ana$reduction_npc_num)
debatch_method     			 	 <-ini.list$integrate_ana$debatch_method
debatch_object     				 <-ini.list$integrate_ana$debatch_object

dims<-reduction_npc_num
npc<-dims+10


#文件获取
#  <-opt$rds
# metafile<-opt$meta
rdsfile<-ini.list$data$f01_rdsfile
metafile<-ini.list$data$f02_metadatafile
#if ( is.null(metafile) )      {metafile  <- 'None'}
if ( metafile=='' )      {metafile  <- 'None'}
metafile


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


#开始进行亚群细分
#nfeatures<-2000;npc<-30;dims<-20 #默认一般参数
if ('seurat_clusters' %in% colnames(rds@meta.data)){
    rds@meta.data$seurat_clusters_raw <- rds@meta.data$seurat_clusters
}
DefaultAssay(rds) <- "RNA"
#https://satijalab.org/seurat/articles/integration_rpca.html
if (debatch_method == 'rpca'){
	ifnb.list <- SplitObject(rds, split.by = debatch_object)
	# normalize and identify variable features for each dataset independently
	ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
		x <- NormalizeData(x)
		x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nvariablefeatures)
	})
	features <- SelectIntegrationFeatures(object.list = ifnb.list)
	ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
		x <- ScaleData(x, features = features, verbose = FALSE)
		x <- RunPCA(x, features = features, verbose = FALSE)
	})
	immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
	# this command creates an 'integrated' data assay
	immune.combined <- IntegrateData(anchorset = immune.anchors)
	# specify that we will perform downstream analysis on the corrected data note that the
	# original unmodified data still resides in the 'RNA' assay
	DefaultAssay(immune.combined) <- "integrated"
	# Run the standard workflow for visualization and clustering
	immune.combined <- ScaleData(immune.combined, verbose = FALSE)
	immune.combined <- RunPCA(immune.combined, npcs = npc, verbose = FALSE)
	# p <- ElbowPlot(immune.combined)
	# ggsave(paste0(prefix,'_PCA_ElbowPlot.pdf'), w=6,h=5)
	immune.combined <- RunTSNE(immune.combined, reduction = "pca",dims = 1:dims)
	immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:dims)
	tmp <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:dims)
	
}else if (debatch_method == 'harmony'){
	#nfeatures<-2000;npc<-30;dims<-20 #默认一般参数
	tmp <- NormalizeData(rds) %>% FindVariableFeatures(nfeatures =nvariablefeatures) %>% ScaleData() %>% RunPCA(npcs = as.numeric(npc),verbose = FALSE)
	tmp <- RunHarmony(tmp, group.by.vars = debatch_object)
	tmp <- RunTSNE(tmp, reduction = "harmony", dims = 1:dims)
	tmp <- RunUMAP(tmp, reduction = "harmony", dims = 1:dims)
	tmp <- FindNeighbors(tmp, reduction = "harmony", dims = 1:dims)#harmony
#CCA方法，亚群细分时一般不建议
}else if (debatch_method == 'cca'){
	ifnb.list <- SplitObject(rds, split.by = debatch_object)
	# normalize and identify variable features for each dataset independently
	ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
		x <- NormalizeData(x)
		x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nvariablefeatures)
	})
	features <- SelectIntegrationFeatures(object.list = ifnb.list)
	# ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
		# x <- ScaleData(x, features = features, verbose = FALSE)
		# x <- RunPCA(x, features = features, verbose = FALSE)
	# })
	immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
	# this command creates an 'integrated' data assay
	immune.combined <- IntegrateData(anchorset = immune.anchors)
	# specify that we will perform downstream analysis on the corrected data note that the
	# original unmodified data still resides in the 'RNA' assay
	DefaultAssay(immune.combined) <- "integrated"
	# Run the standard workflow for visualization and clustering
	immune.combined <- ScaleData(immune.combined, verbose = FALSE)
	immune.combined <- RunPCA(immune.combined, npcs = npc, verbose = FALSE)
	# p <- ElbowPlot(immune.combined)
	# ggsave(paste0(prefix,'_PCA_ElbowPlot.pdf'), w=6,h=5)
	immune.combined <- RunTSNE(immune.combined, reduction = "pca",dims = 1:dims)
	immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:dims)
	tmp <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:dims)

}else if (debatch_method == 'none'){
	tmp <- NormalizeData(rds) %>% FindVariableFeatures(nfeatures =nvariablefeatures) %>% ScaleData() %>% RunPCA(npcs = as.numeric(npc),verbose = FALSE)
	#tmp <- RunHarmony(tmp, group.by.vars = "stim")
	tmp <- RunTSNE(tmp, reduction = "pca", dims = 1:dims)
	tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:dims)
	tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:dims)#harmony
}



#setwd(outdir)

#按照样本绘图
CellDimPlot(srt = tmp, reduction='umap',group.by = 'Sample', label = F,label_insitu=T,label.fg = "black",label.bg = NA,label.size=4,label_repel=F,pt.size=0.01,theme_use="theme_blank")
ggsave(paste0(outdir,'/',debatch_method,'_Sample_UMAP.pdf'),width=6.5,height=5)
CellDimPlot(srt = tmp, reduction='tsne',group.by = 'Sample', label = F,label_insitu=T,label.fg = "black",label.bg = NA,label.size=4,label_repel=F,pt.size=0.01,theme_use="theme_blank")
ggsave(paste0(outdir,'/',debatch_method,'_Sample_TSNE.pdf'),width=6.5,height=5)

#按照组绘图
CellDimPlot(srt = tmp, reduction='umap',group.by = 'Group',label = F,label_insitu=T,label.fg = "black",label.bg = NA,label.size=4,label_repel=F,pt.size=0.01,theme_use="theme_blank")
ggsave(paste0(outdir,'/',debatch_method,'_Group_UMAP.pdf'),width=6.5,height=5)
CellDimPlot(srt = tmp, reduction='tsne',group.by = 'Group',label = F,label_insitu=T,label.fg = "black",label.bg = NA,label.size=4,label_repel=F,pt.size=0.01,theme_use="theme_blank")
ggsave(paste0(outdir,'/',debatch_method,'_Group_TSNE.pdf'),width=6.5,height=5)

saveRDS(tmp,file=paste0(outdir,'/','integrate.rds'))
write.csv(tmp@meta.data,file=paste0(outdir,'/','integrate_metadata.csv'),quote=F,row.names=T)