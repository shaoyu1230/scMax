qc_cluster <- function(seurat.obj, col_sample = 'Sample', col_group = c('Group','Patient')[1],nPCs=30,runtsne=FALSE,outdir='./'){
  if(!dir.exists(outdir)){
    dir.create(outdir,recursive = TRUE)
  }
  # reduction,clustering
  seurat.obj <- NormalizeData(object = seurat.obj, verbose = FALSE)
  print('Normalization finished!')
  seurat.obj <- FindVariableFeatures(object = seurat.obj,verbose = FALSE)
  seurat.obj <- ScaleData(object = seurat.obj,vars.to.regress = "percent.mt",verbose = FALSE)
  print('ScaleData finished!')
  seurat.obj <- RunPCA(object = seurat.obj,verbose = F)
  CellDimPlot(seurat.obj,reduction = "pca",group.by = col_sample) + coord_fixed(1:1)
  ggsave(paste0(outdir,"/1_PCA.pdf"),width = 6,height = 6)
  ggsave(paste0(outdir,"/1_PCA.png"),width = 6,height = 6)
  seurat.obj <- FindNeighbors(object = seurat.obj, reduction = "pca", dims = 1:nPCs, verbose = F)
  seurat.obj <- FindClusters(object = seurat.obj, resolution = 0.2, verbose = F)
  if(runtsne==TRUE){seurat.obj <- RunTSNE(object = seurat.obj, dims = 1:nPCs,check_duplicates = FALSE)}
  seurat.obj <- RunUMAP(object = seurat.obj, dims = 1:nPCs)
  # plots
  dimplots(seurat.obj,reduc='UMAP',outdir=outdir,col_sample=col_sample,col_group=col_group)
  dimplots(seurat.obj,reduc='TSNE',outdir=outdir,col_sample=col_sample,col_group=col_group)
  return(seurat.obj)
}


get_mtgenes <- function(seurat.obj,species=c('human','mouse','other')[3], genelist = NULL){
  # Check if the species is human or mouse
  if (species == "human") {
    mt_genes <- rownames(seurat.obj)[grep("^MT-", rownames(seurat.obj), ignore.case = TRUE)]
  } else if (species == "mouse") {
    mt_genes <- rownames(seurat.obj)[grep("^mt-", rownames(seurat.obj), ignore.case = TRUE)]
  } else if(!is.null(genelist)) {
    mt_genes <- intersect(genelist,rownames(seurat.obj))
    print("Using custom gene list for mitochondrial genes.\n")
    print(paste0("Found ", length(mt_genes)," mitochondrial genes in the custom list.\n"))
    print(mt_genes)
  }
  else {
    print("Species is not 'human', 'mouse', or provide a custom gene list. percent.mt will be set to 0.\n")
    mt_genes <- c()
  }
  return(mt_genes)
}

dimplots <- function(seurat.obj,reduc=c('UMAP','TSNE')[1],outdir='./',col_sample='Sample',col_group=c('Group','Patient')[1]){
  # cluster
  CellDimPlot(seurat.obj,group.by='seurat_clusters',pt.size=0.01,reduction=reduc)
  ggsave(paste0(outdir,'/2_',reduc,'-seurat_clusters.pdf'),width=7,height=5)
  ggsave(paste0(outdir,'/2_',reduc,'-seurat_clusters.png'),width=7,height=5)
  # sample 绘图
  CellDimPlot(seurat.obj,group.by=col_sample,pt.size=0.01,reduction=reduc)
  ggsave(paste0(outdir,'/3_',reduc,'-Sample.pdf'),width=7,height=5)
  ggsave(paste0(outdir,'/3_',reduc,'-Sample.png'),width=7,height=5)
  CellDimPlot(seurat.obj,group.by='seurat_clusters',split.by=col_sample,pt.size=0.01,ncol=3,reduction=reduc)
  nsample <- length(unique(seurat.obj[[col_sample]]))
  ht <- ceiling(nsample/3)*6
  ggsave(paste0(outdir,'/3_',reduc,'-Sample.split.pdf'),width=18,height=ht)
  ggsave(paste0(outdir,'/3_',reduc,'-Sample.split.png'),width=18,height=ht)
  # group 绘图
  if(is.null(col_group) || col_group==''){
    print("col_group is NULL or empty, will not plot group UMAP.")
  }else{
    for (group in col_group) {
      ngroup <- length(unique(seurat.obj[[group]]))
      ht <- ceiling(ngroup/3)*6
      CellDimPlot(seurat.obj,group.by=group,pt.size=0.01,reduction=reduc)
      ggsave(paste0(outdir,'/4_',reduc,'-',group,'.pdf'),width=6,height=5)
      ggsave(paste0(outdir,'/4_',reduc,'-',group,'.png'),width=6,height=5)
      CellDimPlot(seurat.obj,group.by='seurat_clusters',split.by=group,pt.size=0.01,reduction=reduc,ncol=3)
      ggsave(paste0(outdir,'/4_',reduc,'-',group,'.split.pdf'),width=18,height=ht)
      ggsave(paste0(outdir,'/4_',reduc,'-',group,'.split.png'),width=18,height=ht)
    }
  }
}

get_hbgenes <- function(seurat.obj,species=c('human','mouse','other')[3], genelist = NULL){
  # Check if the species is human or mouse
  if (species == "human") {
    hb_genes = rownames(seurat.obj)[grep("^HB[A|B|D|Q|Z]",rownames(seurat.obj),ignore.case = TRUE)]
  } else if (species == "mouse") {
    hb_genes = rownames(seurat.obj)[grep("^HB[A|B|D|Q|Z]",rownames(seurat.obj),ignore.case = TRUE)]
  } else if(!is.null(genelist)) {
    hb_genes <- intersect(genelist,rownames(seurat.obj))
    print("Using custom gene list for Hemoglobin genes.\n")
    print(paste0("Found ", length(hb_genes)," Hemoglobin genes in the custom list.\n"))
    print(hb_genes)
  }
  else {
    print("Species is not 'human', 'mouse', or provide a custom gene list. percent.mt will be set to 0.\n")
    hb_genes <- c()
  }
  return(hb_genes)
}


get_decontx <- function(seurat_data,prefix='./tmp/Sample1'){
  counts <- seurat_data@assays$RNA@counts
  res <- decontX(counts,seed = 2024)
  saveRDS(res,paste0(prefix,'.DecontX_res.rds'))
  seurat_data$Contamination <- res$contamination
  res <- seurat_data@meta.data[,c('orig.ident','Contamination')]
  write.csv(seurat_data@meta.data,file = paste0(prefix,'.DecontX_meta.csv'),row.names = T)
  return(res)
}


get_rate = function(num){
  if (num > 18000){
    rate = 0.14
  }else if(num > 15000){
    rate = 0.1
  }else if(num > 10000){
    rate = 0.076
  }else if (num > 9000){
    rate = 0.069
  }else if (num > 8000){
    rate = 0.061
  }else if (num > 7000){
    rate = 0.054
  }else if (num > 6000){
    rate = 0.046
  }else if (num > 5000){
    rate = 0.039
  }else if (num > 4000){
    rate = 0.031
  }else if (num > 3000){
    rate = 0.023
  }else{
    rate = 0.016
  }
  return(rate)
}
get_doublets = function(obj, rate=NULL, prefix = './tmp/Sample1'){
  if(any(c(is.null(rate), rate == "None"))){
    cellnum = ncol(obj)
    rate = get_rate(cellnum)
  }
  
  # parameters from DoubleFinder_lmSeuratV4.R
  pcs <- 1:20
  pN <- 0.25
  pK_search_pcs <- pcs
  cluster_resolution <- 0.5
  
  obj = NormalizeData(obj)
  obj = FindVariableFeatures(obj)
  obj = ScaleData(obj)
  obj = RunPCA(obj)
  obj = FindNeighbors(obj, dims = pcs)
  obj = RunUMAP(obj, dims = 1:10)
  obj = FindClusters(obj, resolution = cluster_resolution)
  
  message("\n[", Sys.time(), "]", " Start pK sweep!")
  sweep.res.list <- paramSweep(obj, PCs = pK_search_pcs, sct=FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  
  message("\n[", Sys.time(), "]", " Find pk number!")
  bcmvn <- find.pK(sweep.stats)
  best.pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  message("最佳 pK = ", best.pK)
  
  message("\n[", Sys.time(), "]", " Calculate homotypic proportion & expected doublets!")
  nCells <- ncol(obj)
  nExp <- round(rate * nCells)
  homotypic.prop <- modelHomotypic(obj$seurat_clusters)
  nExp.adj <- round(nExp * (1 - homotypic.prop))
  
  message("\n[", Sys.time(), "]", " Find doublet cell!")
  obj <- doubletFinder(obj, 
                       PCs = pcs, 
                       pN = pN, 
                       pK = best.pK, 
                       nExp = nExp.adj, 
                       reuse.pANN = FALSE,
                       sct = FALSE)
                       
  df.col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE)
  pANN.col <- grep("^pANN", colnames(obj@meta.data), value = TRUE)
  
  saveRDS(obj@meta.data, paste0(prefix, ".doubletFinder_meta.rds"))
  res = obj@meta.data[, c(pANN.col, df.col)]
  colnames(res) = c("doublets_score", "DF.classifications")
  return(res)
}

run_DF_for_sample <- function(seu, sample_name,
                              expectedDoubletRate = 0.075,
                              pcs = 1:20,
                              pN = 0.25,
                              pK_search_pcs = 1:20,
                              cluster_resolution = 0.5,
                              sample_col = "orig.ident") {

  message("\n==============================")
  message("开始处理样本：", sample_name)
  message("==============================")

  # 避免 subset 死板报错，使用安全的列提取方式
  cells_use <- colnames(seu)[seu@meta.data[[sample_col]] == sample_name]
  if(length(cells_use) == 0){
    warning(paste0("样本 ", sample_name, " 细胞数为0或标识未找到，跳过..."))
    return(seu)
  }
  sub <- subset(seu, cells = cells_use)
  message("细胞数：", ncol(sub))

  sub <- NormalizeData(sub)
  sub <- FindVariableFeatures(sub)
  sub <- ScaleData(sub)
  sub <- RunPCA(sub)
  sub <- FindNeighbors(sub, dims = pcs)
  sub <- RunUMAP(sub, dims = 1:10)
  sub <- FindClusters(sub, resolution = cluster_resolution)

  #### pK sweep
  message("  pK sweep ...")
  sweep.res <- paramSweep(sub, PCs = pK_search_pcs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  pK_table <- find.pK(sweep.stats)
  best.pK <- as.numeric(as.character(
    pK_table$pK[which.max(pK_table$BCmetric)]
  ))
  message("最佳 pK = ", best.pK)

  #### homotypic correction
  nCells <- ncol(sub)
  nExp <- round(expectedDoubletRate * nCells)
  homotypic.prop <- modelHomotypic(sub$seurat_clusters)
  nExp.adj <- round(nExp * (1 - homotypic.prop))

  #### DoubletFinder
  message("  运行 DoubletFinder ...")
  sub <- doubletFinder(
    sub,
    PCs = pcs,
    pN = pN,
    pK = best.pK,
    nExp = nExp.adj,
    reuse.pANN = FALSE,
    sct = FALSE
  )

  df.col <- grep("^DF.classifications", colnames(sub@meta.data), value = TRUE)
  pANN.col <- grep("^pANN", colnames(sub@meta.data), value = TRUE)

  cells <- colnames(sub)
  seu@meta.data[cells, "DF_result"] <- sub[[df.col]]
  seu@meta.data[cells, "PANN_score"] <- sub[[pANN.col]]

  return(seu)
}

##############################
#        循环样本分析
##############################

# for (s in samples) {
#   message(paste0("开始分析：", s))
#   seu <- run_DF_for_sample(seu, s)
# }



