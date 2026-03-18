# 共有基因富集分析
library(stringr)
GeneList2GOKEGG <- function(infile,orgdb='org.Rn.eg.db',organism_kegg='rno',
                            outdir,prefix='',runGO=TRUE,runKEGG=TRUE,title=''){
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
