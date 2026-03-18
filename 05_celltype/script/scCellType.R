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
suppressMessages(library(yaml))

bindir <- this.path::this.dir()
source(file.path(bindir,'func_scRNA_celltype_anno.R'))

parser = argparse::ArgumentParser(description = 'Unified script for cluster plots, ref markers, annotation and celltype plots.')
parser$add_argument('-c', '--config', dest = 'config', required = TRUE, help = 'main YAML config file')
parser$add_argument('-o', '--outdir', dest = 'outdir', required = TRUE, help = 'output directory')
parser$add_argument('--inputrds', dest = 'inputrds', default = '', help = 'override input rds path (optional)')
parser$add_argument('--annotated_rds', dest = 'annotated_rds', default = '', help = 'use an existing annotated rds for celltype plots/DE')
parser$add_argument('--do_cluster', dest = 'do_cluster', action = 'store_true', help = 'run cluster plots and cluster DE/enrich')
parser$add_argument('--do_refmarker', dest = 'do_refmarker', action = 'store_true', help = 'run reference marker plots')
parser$add_argument('--do_annotation', dest = 'do_annotation', action = 'store_true', help = 'run manual annotation')
parser$add_argument('--do_celltype', dest = 'do_celltype', action = 'store_true', help = 'run celltype umap plots and fraction plots')
parser$add_argument('--do_celltype_de', dest = 'do_celltype_de', action = 'store_true', help = 'run celltype DE and enrichment')

opt = parser$parse_args()

cfg <- yaml::read_yaml(opt$config)
ct_conf <- cfg[['05_celltype']]
if (is.null(ct_conf)) {
  stop('No 05_celltype section found in main config!')
}

db_conf <- cfg[['Database_Info']]

outdir <- opt$outdir
if(!dir.exists(outdir)){dir.create(outdir, recursive = TRUE)}

inputrds <- opt$inputrds
if (is.null(inputrds) || inputrds == '') {
  inputrds <- ct_conf$rdsfile
}
annotated_rds <- opt$annotated_rds

annofile <- ct_conf$annofile
clustercol <- ct_conf$cluster_col
celltype_col <- ct_conf$celltype_col
sample_col <- ct_conf$col_sample
group_col <- ct_conf$col_group
sample_order <- ct_conf$sample_order
group_order <- ct_conf$group_order
cluster_colors <- ct_conf$cluster_colors
celltype_levels <- ct_conf$celltype_levels
celltype_colors <- ct_conf$celltype_colors
refmarker.file <- ct_conf$refmarker_file
metacsv <- ct_conf$metacsv

if (!is.null(group_col)) {
  group_col <- strsplit(group_col, split = ',')[[1]]
}

species <- if (!is.null(db_conf) && !is.null(db_conf$species)) db_conf$species else "Human"

# 支持用户从 config 读取直接外置定义的 orgdb 和 org_kegg（适用于人/鼠外的杂合物种）
orgdb_custom <- if (!is.null(ct_conf) && !is.null(ct_conf$orgdb)) ct_conf$orgdb else ""
kegg_custom <- if (!is.null(ct_conf) && !is.null(ct_conf$organism_kegg)) ct_conf$organism_kegg else ""

if (orgdb_custom != "" && kegg_custom != "") {
  orgdb <- orgdb_custom
  organism_kegg <- kegg_custom
} else if (tolower(species) == "mouse") {
  orgdb <- "org.Mm.eg.db"
  organism_kegg <- "mmu"
} else {
  # 兜底 Human
  orgdb <- "org.Hs.eg.db"
  organism_kegg <- "hsa"
}

do_cluster <- isTRUE(opt$do_cluster) || isTRUE(ct_conf$do_cluster)
do_refmarker <- isTRUE(opt$do_refmarker) || isTRUE(ct_conf$do_refmarker)
do_annotation <- isTRUE(opt$do_annotation) || isTRUE(ct_conf$do_annotation)
do_celltype <- isTRUE(opt$do_celltype) || isTRUE(ct_conf$do_celltype)
do_celltype_de <- isTRUE(opt$do_celltype_de) || isTRUE(ct_conf$do_celltype_de)

if (!(do_cluster || do_refmarker || do_annotation || do_celltype || do_celltype_de)) {
  do_cluster <- TRUE
  do_refmarker <- TRUE
  do_annotation <- TRUE
  do_celltype <- TRUE
  do_celltype_de <- TRUE
}

if ((do_celltype || do_celltype_de) && !do_annotation && (is.null(annotated_rds) || annotated_rds == '')) {
  stop('do_celltype/do_celltype_de requires do_annotation or --annotated_rds.')
}

if (is.null(inputrds) || inputrds == '') {
  stop('inputrds is required (provide --inputrds or set 05_celltype.rdsfile).')
}

########## 0_preparations ################
if(!dir.exists(paste0(outdir,'/Rdata'))){dir.create(paste0(outdir,'/Rdata'),recursive = TRUE)}
if(!dir.exists(paste0(outdir,'/tmp'))){dir.create(paste0(outdir,'/tmp'),recursive = TRUE)}
if(do_cluster && !dir.exists(paste0(outdir,'/cluster_characterization'))){dir.create(paste0(outdir,'/cluster_characterization'))}
if(do_refmarker && !dir.exists(paste0(outdir,'/marker_expression'))){dir.create(paste0(outdir,'/marker_expression'))}
if(do_annotation && !dir.exists(paste0(outdir,'/annotation'))){dir.create(paste0(outdir,'/annotation'))}
if(do_celltype && !dir.exists(paste0(outdir,'/celltype_fraction'))){dir.create(paste0(outdir,'/celltype_fraction'))}
if(do_celltype && !dir.exists(paste0(outdir,'/celltype_characterization'))){dir.create(paste0(outdir,'/celltype_characterization'))}

if(!file.exists(inputrds)){
  stop(paste0(inputrds,' file not exists!'))
}

print(Sys.time())
data <- readRDS(inputrds)
DefaultAssay(data) <- "RNA"
if(!is.null(metacsv) && metacsv != '' && file.exists(metacsv)){
  meta <- read.csv(metacsv,stringsAsFactors = FALSE, row.names=1)
  data <- AddMetaData(data, meta)
}

metaraw <- data@meta.data

if(!is.null(sample_order) && sample_order != ''){
  sample_order <- strsplit(sample_order, split = ',')[[1]]
  data@meta.data[,sample_col] <- as.character(data@meta.data[,sample_col])
  data@meta.data[,sample_col] <- factor(data@meta.data[,sample_col], levels = sample_order)
}
if (sample_col != 'Sample'){
  data$Sample <- data@meta.data[,sample_col]
}

if(!is.null(group_order) && group_order != '' && !is.null(group_col)){
  gcol <- group_col[1]
  group_order <- strsplit(group_order, split = ',')[[1]]
  if (length(group_order[group_order %in% unique(data@meta.data[,gcol])]) != length(unique(data@meta.data[,gcol]))){
    stop(paste0('group_order not same as ',gcol,' in meta.data'))
  }
  data@meta.data[,gcol] <- as.character(data@meta.data[,gcol])
  data@meta.data[,gcol] <- factor(data@meta.data[,gcol], levels = group_order)
}

if (!is.null(group_col)){
  if (group_col[1] != 'Group'){
    data$Group <- data@meta.data[,group_col[1]]
  }
  plot_groups <- unique(c('Sample', group_col))
} else {
  plot_groups <- 'Sample'
}

columns <- colnames(data@meta.data)
columns <- columns[columns %in% plot_groups]
if(length(plot_groups) > length(columns)){
  stop('Not all plot_groups in the data!')
}

if(!is.null(clustercol) && clustercol != ''){
  if(clustercol %in% colnames(data@meta.data)){
    data$seurat_clusters <- data@meta.data[,clustercol]
  } else {
    stop(paste0(clustercol, ' Not in the metadata!'))
  }
}
Idents(data) <- data$seurat_clusters

if(!is.null(cluster_colors) && cluster_colors != ''){
  cluster_colors <- strsplit(cluster_colors,"\\s*,\\s*")[[1]]
}
ncluster <- length(unique(data$seurat_clusters))
if(length(cluster_colors) < ncluster){
  cluster_colors <- NULL
}

if(!is.null(celltype_colors) && celltype_colors != ''){
  celltype_colors <- strsplit(celltype_colors,"\\s*,\\s*")[[1]]
}

if(!is.null(celltype_levels) && celltype_levels != ''){
  celltype_levels <- strsplit(celltype_levels, "\\s*,\\s*")[[1]]
}

if(do_cluster){
  print(paste0(Sys.time(),' - cluster_characterization'))
  cluster_plots(data, outputDir=paste0(outdir,'/cluster_characterization'), groups=plot_groups, palcolor=cluster_colors)

  library(SeuratWrappers)
  diffgene <- cluster_de(data, rdsdir=file.path(outdir,'/Rdata'),
                         outdir=file.path(outdir,'cluster_characterization/4_Cluster.DE'))
  cluster_enrich(derds=file.path(outdir,'/Rdata/ClusterDEGs.rds'),
                 orgdb=orgdb, organism_kegg=organism_kegg,
                 enrichdir=paste0(outdir,'/cluster_characterization/5_Enrich'))
}

if(do_refmarker){
  if(is.null(refmarker.file) || refmarker.file == '' || !file.exists(refmarker.file)){
    print(paste0('refmarker file not exists: ', refmarker.file))
  } else {
    plot_refmarkers(data, refmarker.file, outdir=paste0(outdir,'/marker_expression'))
    files <- list.files(paste0(outdir,'/marker_expression'),pattern='2_FeaturePlot')
    example.file <- files[grep('png',files)][1]
    file.copy(paste0(outdir,'/marker_expression/',example.file),paste0(outdir,'/marker_expression/2_FeaturePlot.example.png'),)
  }
}

data_anno <- NULL
data.f <- NULL
celltype.keep <- NULL

if(!is.null(annotated_rds) && annotated_rds != ''){
  if(!file.exists(annotated_rds)){
    stop(paste0(annotated_rds,' file not exists!'))
  }
  data_anno <- readRDS(annotated_rds)
  DefaultAssay(data_anno) <- "RNA"
}

if(do_annotation){
  if(is.null(annofile) || annofile == '' || !file.exists(annofile)){
    stop('annofile is required for annotation.')
  }

  data_anno <- add_anno(data,outdir=paste0(outdir,'/annotation'),annofile,
                          celltype_levels,celltype_col=celltype_col)
  metaraw[,celltype_col] <- data_anno@meta.data[,celltype_col]

  ncts <- length(unique(data_anno@meta.data[,celltype_col]))
  if(length(celltype_colors) < ncts){
    celltype_colors <- palette_scp(unique(data_anno@meta.data[,celltype_col]),palette='Paired',NA_keep=TRUE)
  }

  celltype_umap_plots(data_anno,outdir=paste0(outdir,'/annotation/CellType_All'),celltype_col=celltype_col,groups=plot_groups,palcolor=celltype_colors)
  df <- read.table(annofile,sep='\t',header=TRUE)
  cols <- c('Cluster','CellType','Markers','reference','description')
  #判断是否有reference和description列
  if(!'reference' %in% colnames(df)){
    cols <- cols[-4]
  }
  if(!'description' %in% colnames(df)){
    cols <- cols[-5]
  } 
  df <- df[,cols]
  png(paste0(outdir,'/annotation/CellType_All/0_Cluster-CellType.anno.png'), width = 2000, height = 1600, res = 150)
  grid.table(df)
  dev.off()
  file.copy(annofile,paste0(outdir,'/annotation/CellType_All/0_Cluster-CellType.anno.xls'))

  anno.df <- read.csv(annofile,sep='\t')
  if('Markers' %in% colnames(anno.df)){
    anno_mks_plot1(data_anno,annofile,celltype_col=celltype_col,celltype.levels=celltype_levels,outdir=paste0(outdir,'/annotation/CellType_All'))
  }

  celltype.keep <- ct_conf$celltype_keep
  if (is.null(celltype.keep) || celltype.keep == '') {
    celltype.keep <- unique(anno.df$CellType)
  } else {
    celltype.keep <- strsplit(celltype.keep, "\\s*,\\s*")[[1]]
  }

  keep.cells <- rownames(data_anno@meta.data[data_anno@meta.data[,celltype_col] %in% celltype.keep,])
  data.f <- subset(data_anno,cells=keep.cells)
  data.f@meta.data[,celltype_col] <- droplevels(data.f@meta.data[,celltype_col])
  metaraw.f <- metaraw[keep.cells,]

  celltype_umap_plots(data.f,outdir=paste0(outdir,'/annotation/CellType_Keep'),celltype_col=celltype_col,groups=plot_groups,palcolor=celltype_colors)
  if('Markers' %in% colnames(anno.df)){
    anno_mks_plot2(data.f,annofile,celltype_col=celltype_col,celltype.keep=celltype.keep,outdir=paste0(outdir,'/annotation/CellType_Keep'))
  }

  saveRDS(data_anno,paste0(outdir,'/Rdata/Data-Annotation_CellType.rds'))
  saveRDS(data.f,paste0(outdir,'/Rdata/Data-Annotation_CellType.Keep.rds'))
  write.csv(metaraw.f,paste0(outdir,'/Rdata/Meta-Annotation_CellType.Keep.csv'),row.names = FALSE)
  write.csv(metaraw,paste0(outdir,'/Rdata/Meta-Annotation_CellType.csv'),row.names = FALSE)
}

if(do_celltype){
  if (is.null(data_anno)) {
    stop('do_celltype requires annotated data.')
  }
  if (is.null(data.f)) {
    if (is.null(celltype.keep) || length(celltype.keep) == 0) {
      celltype.keep <- unique(data_anno@meta.data[,celltype_col])
    }
    keep.cells <- rownames(data_anno@meta.data[data_anno@meta.data[,celltype_col] %in% celltype.keep,])
    data.f <- subset(data_anno,cells=keep.cells)
    data.f@meta.data[,celltype_col] <- droplevels(data.f@meta.data[,celltype_col])
  }
  if(length(celltype_colors) < length(unique(data_anno@meta.data[,celltype_col]))){
    celltype_colors <- palette_scp(unique(data_anno@meta.data[,celltype_col]),palette='Paired',NA_keep=TRUE)
  }
  celltype_umap_plots(data_anno,outdir=paste0(outdir,'/annotation/CellType_All'),celltype_col=celltype_col,groups=plot_groups,palcolor=celltype_colors)
  celltype_umap_plots(data.f,outdir=paste0(outdir,'/annotation/CellType_Keep'),celltype_col=celltype_col,groups=plot_groups,palcolor=celltype_colors)
  cell_fraction_plots(data_anno,out_prefix=paste0(outdir,'/celltype_fraction/1_CellType_All.in.'),
                      celltype_col=celltype_col,groups=plot_groups,do.boxplot=FALSE,palcolor=celltype_colors)
  cell_fraction_plots(data.f,out_prefix=paste0(outdir,'/celltype_fraction/2_CellType_Keep.in.'),
                      celltype_col=celltype_col,groups=plot_groups,do.boxplot=FALSE,palcolor=celltype_colors)
  roe_groups <- intersect(group_col,colnames(data.f@meta.data))
  if (!is.null(roe_groups)){
    for(group in roe_groups){
      Roe <- distribution_Roe(data.f@meta.data,celltype_column = celltype_col,celltype_level = celltype.keep,
                              condition_column = group,out_prefix=paste0(outdir,'/celltype_fraction/3_CellType_Keep.in.',group,'-'))
      ggsave(paste0(outdir,'/celltype_fraction/3_CellType_Keep.in.',group,'-Roe.pdf'),width=5,height=5)
      ggsave(paste0(outdir,'/celltype_fraction/3_CellType_Keep.in.',group,'-Roe.png'),width=5,height=5)
    }
  }

}

if(do_celltype_de){
  if (is.null(data.f)) {
    stop('do_celltype_de requires annotated data.')
  }
  library(SeuratWrappers)
  diffgene <- celltype_de(data.f,rdsdir=file.path(outdir,'/Rdata'),celltype_col=celltype_col,outdir=paste0(outdir,'/celltype_characterization/1_CellType_DE'))
  files <- list.files(paste0(outdir,'/celltype_characterization/1_CellType_DE'),pattern='4_FeatureExpr')
  example.file <- files[grep('png',files)][1]
  file.copy(paste0(outdir,'/celltype_characterization/1_CellType_DE/',example.file),
            paste0(outdir,'/celltype_characterization/1_CellType_DE/4_FeatureExprUMAP.example.png'),)

  celltype_enrich(derds=paste0(outdir,'/Rdata/CellTypeDEGs.rds'),
                  orgdb=orgdb,organism_kegg=organism_kegg,
                  enrichdir=paste0(outdir,'/celltype_characterization/2_CellType_Enrich'))
}

print(paste0(Sys.time(),' - scCellType done'))
