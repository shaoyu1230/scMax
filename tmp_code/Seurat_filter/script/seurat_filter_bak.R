library('getopt')
library(Seurat)
library(dplyr)
library(Matrix)
library(magrittr)
library(tibble)
library(readr)
library(patchwork)  
library(ggplot2)
library(SCP)

para<-matrix(c(
	'help',	'h',	0,	"logical",
	'datadir',	'd',	1,	"character",
	'metadata', 'f',     1,  "character",
	'outdir',	'o',	1,	"character",
	'minc',	'c',	2,	"integer",
	'maxc',	'e',	2,	"integer",
	'min',	'm',	2,	"integer",
	'max',	'a',	2,	"integer",
	'mtp',	't',	2,	"integer",
	'hbp',	'b',	2,	"integer",
	'Dou',	'D',	2,	"character",
	'decontX',	'X',	2,	"character"
		),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
		cat("
		Usage example:
			Rscript securat_parameter.R -d hg19/ -p ANproject -o outdir/
		Options:
			--help		h	NULL		get this help
			--datadir	d	character	the cellRanger output datadir [forced]
			--metadata  f   character   the metadatafile 
			--outdir	o	character	output file dir [forced]
			--min	min	character	Filter cells based on minium gene numbers [default: 200]
			--max	max	character	 Filter cells based on max gene numbers [default: 10000]
			--mtp	mtp	character	Filter cells based on mito percent [default:20]
			--hbp	hbp	character	Filter cells based on HB percent [default:5]
			--minc	minc	character	Filter cells based on min umi 
			--maxc	maxc	character	Filter cells based on max umi
			--Dou	Dou	character	Filter cells based on double cells
			--decontX	X	character	Filter cells based on decontX
			\n")
			q(status=1)
			}


print(1)

if ( is.null(opt$min))	{ opt$min <- 200 }
if ( is.null(opt$max))	{ opt$max <- 1000000000 }
if ( is.null(opt$mtp))	{ opt$mtp <- 20 }
if ( is.null(opt$hbp))	{ opt$hbp <- 5}
if ( is.null(opt$minc))	{ opt$minc <- 200}
if ( is.null(opt$maxc))	{ opt$maxc <- 1000000000}
if ( is.null(opt$Dou))	{ opt$Dou <- "False"}
if ( is.null(opt$decontX))	{opt$decontX <- 1}

print(2)

A<<-opt$min
B<<-opt$max
CC<<-opt$minc
DD<<-opt$maxc
mt_percent<<-opt$mtp
HB_percent<<-opt$hbp
Doublet <<- opt$Dou
decontX <<- opt$decontX

print(3)
##读入RDS，并进行metadata的核查和替换
object <- readRDS(opt$datadir)
if ( is.null(opt$metadata)){
	warning("本次分析未提供新的metadata，将使用rds中原始的metadata信息")
}else{
	new_meta <- read.csv(opt$metadata,row.names = 1) #%>% column_to_rownames("Barcode")
	if (!identical(sort(rownames(new_meta)), sort(rownames(object@meta.data)))){
		missing_in_new <- setdiff(rownames(object@meta.data),rownames(new_meta))
		missing_in_old <- setdiff(rownames(new_meta),rownames(object@meta.data))
		cat("错误：细胞ID不匹配！\n")
		cat("新metadata数据中缺失的细胞:", paste(head(missing_in_new), collapse = ","), "...\n")
		cat("原数据中缺失的细胞:", paste(head(missing_in_old), collapse = ","), "...\n")
		stop("metadata数据核查失败，请检查细胞ID一致性")
	}else{
		cat("√细胞ID完全匹配\n")
	}
	new_meta <- new_meta[rownames(object@meta.data), ,drop = FALSE]
	object@meta.data <- new_meta
}



if (!("nCount_RNA" %in% names(object@meta.data))){stop("nCount_RNA列不在您的rds中，请使用安诺云工具Seurat质控来规范您的rds文件")}
if (!("nFeature_RNA" %in% names(object@meta.data))){stop("nFeature_RNA列不在您的rds中，请使用安诺云工具Seurat质控来规范您的rds文件")}
if (!("percent.mt" %in% names(object@meta.data))){warning("percent.mt列不在您的rds中，本次将不对线粒体比例进行过滤！！！")}
if (!("percent.hb" %in% names(object@meta.data))){warning("percent.hb列不在您的rds中，本次将不对红细胞比例进行过滤！！！")}
if (!("DF.classifications" %in% names(object@meta.data))){warning("DF.classifications列不在您的rds中，无法对双细胞进行过滤！！！")}
if (!("Contamination" %in% names(object@meta.data))){warning("Contamination列不在您的rds中，无法对rRNA污染进行过滤！！！")}


print('start filter')

object1 <- subset(object,subset = nFeature_RNA > A & nFeature_RNA < B &  nCount_RNA > CC &  nCount_RNA < DD )
if ("percent.mt" %in% names(object@meta.data)) {object1 <- subset(object1,subset = percent.mt < mt_percent)}
if ("percent.hb" %in% names(object@meta.data)) {object1 <- subset(object1,subset = percent.hb < HB_percent)}
if (Doublet =="True" & "DF.classifications" %in% names(object@meta.data)){object1 <- subset(object1, subset = DF.classifications == "Singlet")}
if ("Contamination" %in% names(object@meta.data)) {object1 <- subset(object1,subset = Contamination < decontX)}


file_stat <- paste(opt$outdir,'filter_stat_cell.csv',sep='/')

samples <- unique(object@meta.data$Sample)
del_cells_list <- list()

for (sample in samples) {
	sample_cells <- rownames(object@meta.data)[object@meta.data$Sample == sample]
	sample_meta <- subset(object, cells = sample_cells)
	#sample_meta <- object@meta.data[sample_cells, ]
	sample_cells_filtered <- intersect(sample_cells, rownames(object1@meta.data))
	sample_data <- data.frame(
	    sample = sample,
		total_cell = length(sample_cells),
		remaining_cell = length(sample_cells_filtered),
		low_nFeature = sum(sample_meta$nFeature_RNA <= A, na.rm = TRUE),
		high_nFeature = sum(sample_meta$nFeature_RNA >= B, na.rm = TRUE),
		low_nCount = sum(sample_meta$nCount_RNA <= CC, na.rm = TRUE),
		high_nCount = sum(sample_meta$nCount_RNA >= DD, na.rm = TRUE)
	)
	if ("percent.mt" %in% names(sample_meta@meta.data)) {sample_data$high_percent.mt <- as.character(nrow(subset(sample_meta@meta.data,percent.mt >= mt_percent)))}else{sample_data$high_percent.mt <-0}
	if ("percent.hb" %in% names(sample_meta@meta.data)) {sample_data$high_HB <- as.character(nrow(subset(sample_meta@meta.data,percent.hb >= HB_percent)))}else{sample_data$high_HB<-0}
	if (Doublet =="True" & "DF.classifications" %in% names(object@meta.data)) {sample_data$DF.classifications<- as.character(nrow(subset(sample_meta@meta.data,DF.classifications == "Doublet")))}else{sample_data$Doublet<-0}
	if ("Contamination" %in% names(sample_meta@meta.data)){sample_meta$Contamination <- as.character(nrow(subset(sample_meta@meta.data,Contamination >= decontX)))}else{sample_data$Contamination <- 0}
	sample_data$all_filtered_cell <- length(sample_cells)-length(sample_cells_filtered)
	#removed_cells<<-rbind(removed_cells,sample_data)
	del_cells_list[[sample]] <- sample_data
}
del_cells <- do.call(rbind, del_cells_list)
rownames(del_cells) <- NULL

write.csv(del_cells,file_stat,quote=F,row.names=F)


print('start plot')


qc_b <- FeatureStatPlot(object, stat.by = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb","Contamination"), group.by = 'Sample',ncol=1,stack = TRUE,add_box=TRUE) & NoLegend()&theme(axis.title.y = element_blank())
print(qc_b)
#ggsave(file.path(opt$outdir,"Vln_beforQC.pdf"),width = 12,height = 7)
ggsave(file.path(opt$outdir,"Vln_beforeQC.pdf"),width = 7+0.2*length(unique(object$Sample)),height =7 )
ggsave(file.path(opt$outdir,"Vln_beforeQC.png"),width = 7+0.2*length(unique(object$Sample)),height =7)

qc_a <- FeatureStatPlot(object1, stat.by = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb","Contamination"), group.by = 'Sample',ncol=1,stack = TRUE,add_box=TRUE) & NoLegend()&theme(axis.title.y = element_blank())
print(qc_a)
ggsave(file.path(opt$outdir,"Vln_afterQC.pdf"),width = 7+0.2*length(unique(object$Sample)),height = 7)
ggsave(file.path(opt$outdir,"Vln_afterQC.png"),width = 7+0.2*length(unique(object$Sample)),height = 7)

print('start saverds and metadata')

saveRDS(object1,paste0(opt$outdir,'/final_obj.rds'))

metadata_df <- object1@meta.data
#metadata_df <- rownames_to_column(metadata_df, var = "Barcode")
metadata_file<-paste0(opt$outdir,"/final_metadata.csv")
write.csv(metadata_df, file = metadata_file, row.names = TRUE, quote = FALSE)

print('finished all analysis!')
