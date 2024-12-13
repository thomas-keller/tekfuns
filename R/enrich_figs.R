



#' Filter a gene set enrichment result to eliminate genes with low absolute fold-change
#' @importFrom magrittr %>%
#' @param gse gene-set enrichment object from clusterprofiler
#' @param fcut minimum absolute foldchange of genes to be retained
#' @param showcat number of gene pathways to plot
#' @return cnetplot with filtered genes
#' @export
#' @examples todo
cnetfilt <- function(gse,fcut=2.0,showcat=3){
  gl=gse@geneList
  gl=gl[abs(gl)>=fcut]
  core.genes <- strsplit(as.data.frame(gse)[,"core_enrichment"] , "/")
  fcore<-sapply(lapply(core.genes,\(x) x[x %in% names(gl)] ),paste,collapse="/")
  fcore2=fcore[which(fcore!="")]

  g2=gse
  g2@result=g2@result[which(fcore!=""),]
  g2@result$core_enrichment=fcore2
  g2=clusterProfiler::setReadable(g2,org.Hs.eg.db,"ENTREZID")
  p=enrichplot::cnetplot(g2,color.params=list(foldChange=gl),showCategory=showcat)
  #p=cnetplot(g2)
  return(p)
}


#' automagically make a bunch of standard deseq2 and clusterprofiler plots
#'
#' @param dds DESeqDataset object
#' @param pcut gsea pvalue threshold
#' @param folder output result folder
#' @param fprefix output analysis prefix
#'
#' @return list of plots
#' @export
#'
#' @examples todo
rnaplots <- function(dds,pcut=0.05,fcut=2,folder=NULL,fprefix=NULL){
  if(!is.null(folder)){
    dir_create(folder)
  }
  rres=list()

  #filter out genes with low reads
  smallestGroupSize <- 3
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]

  #run deseq method, it will automatically perform several steps
  dds<-DESeq2::DESeq(dds)

  #extract the results
  res <- DESeq2::results(dds)
  #do some cleaning and filter out nonsignificant genes
  res05 <- na.omit(res)
  res05 <-as.data.frame(res05)
  res05 <- dplyr::filter(res05,padj<=0.05)
  #res05 <- dplyr::arrange(res05,desc(abs(log2FoldChange)))
  res05 <- dplyr::arrange(res05,padj)
  res05$ens=rownames(res05)
  ens=res05$ens
  res05$symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ens, column = c('SYMBOL'), keytype = 'ENSEMBL')
  res05$entrez <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ens,column = c('ENTREZID'),
                         keytype = 'ENSEMBL')
  res05$csymbol=ifelse(is.na(res05$symbol),res05$ens,res05$symbol)
  rres$res05=res05
  vsd=DESeq2::vst(dds,blind=TRUE)
  #pca
  p=DESeq2::plotPCA(vsd,intgroup="condition")
  rres$pca=p
  #dispersal estimate
  #its a base R plot, move to end
  #volcano
  p=EnhancedVolcano::EnhancedVolcano(res05,
                                   lab = res05$csymbol,
                                   x = 'log2FoldChange',y = 'pvalue',pCutoff = 10e-12, title = NULL,
                                   FCcutoff = 1.5,
                                   subtitle=NULL,
                                   cutoffLineType = 'twodash',
                                   cutoffLineWidth = 0.8,
                                   pointSize = 1.0,
                                   labSize = 3.0,
                                   colAlpha = 1,
                                   legendLabels=c('Not sig.','Log (base 2) FC','p-value','p-value & Log (base 2) FC'),
                                   legendPosition = 'right',
                                   legendLabSize = 10,
                                   legendIconSize = 5.0)
  rres$volc=p

  #heatmap
  #res05 arranged by padj now
  topg <- res05$ens[1:30]
  mat  <- SummarizedExperiment::assay(vsd)
  #mat  <- mat[row.names(mat) %in% topg, ]
  #do z scale
  mat <- t(scale(t(mat)))
  #resv=res05[row.names(res05) %in% row.names(mat),]
  #row.names(mat)=resv$symbol
  anno <- as.data.frame(colData(vsd)[, c("condition")])
  colnames(anno) <- "condition"
  #add rownames to anno so they match the column names from the expression matrix
  #this is needed so that the two dataframes can be matched up
  row.names(anno) <- colnames(mat)
  #row.names(mat)<-res05$csymbol[1:30]
  #0 point not always white, how to force:
  #thanks https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
  cols <- RColorBrewer::brewer.pal(10, "Spectral")
  plen<-50
  cols<- grDevices::colorRampPalette(rev(cols))(plen)
  mat=mat[1:50,]
  myBreaks <- c(seq(min(mat), 0, length.out=ceiling(plen/2) + 1),
                seq(max(mat)/plen, max(mat), length.out=floor(plen/2)))
  p=pheatmap::pheatmap(mat,annotation_col=anno,scale='row',show_rownames=FALSE,silent=TRUE,breaks=myBreaks,color=cols)
  rres$heatmap=p
  coldata=as.data.frame(colData(vsd))
  #gene dotplot
  tmat=DESeq2::counts(dds,normalized=TRUE)
  tmat=tmat[row.names(tmat) %in% topg[1:30],]
  sym=res05$csymbol[1:30]
  #row.names(tmat)=resv$symbol
  tmatd=as.data.frame(tmat)
  tmatd$symbol=sym
  expr_long=tidyr::pivot_longer(tmatd,-symbol,names_to='names',values_to='expression')
  expr_long=dplyr::inner_join(expr_long,coldata)
  expr_long$expression=expr_long$expression+1
  p=ggplot2::ggplot(expr_long,aes(x=symbol,y=expression,colour=condition))+geom_point()
  p=p+ggplot2::scale_x_discrete(guide=guide_axis(n.dodge=3))
  p=p+ggplot2::theme_bw()
  #convert y axis to log scale since values range over orders of magnitude
  p=p+ggplot2::scale_y_log10()
  #rename x and y labels
  p=p+ggplot2::xlab("gene symbol")+ggplot2::ylab("log10 counts")
  rres$topgdot=p
  #now to enrich plots
  fc=res05$log2FoldChange
  names(fc)=res05$entrez
  fc=sort(fc,decreasing=TRUE)

  m_df <- msigdbr::msigdbr(species = "Homo sapiens")
  m_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene)
  m_t2g$gs_name=stringr::str_remove(m_t2g$gs_name,"HALLMARK_")
  m_t2g$gs_name=stringr::str_remove(m_t2g$gs_name,"HALLMARK ")
  #replace the underscores with spaces to make wrapping easier
  m_t2g$gs_name=stringr::str_replace_all(m_t2g$gs_name,"_"," ")
  set.seed(1234)

  #GSEA function takes minimum two arguments; the fold-change values, and the pathway-gene database
  #if you take a look at the m_t2g dataframe we constructed from msigdbr, its has two columns, one with the pathway name (gs_name)
  #and one with the gene id (entrez_gene)

  gres=clusterProfiler::GSEA(fc,TERM2GENE=m_t2g,pvalueCutoff=pcut)
  #only use original enrich result for filtering genes with cnetfilt
  greso=gres
  #use setReadable to convert entrez ids to gene names
  gres=clusterProfiler::setReadable(gres,OrgDb=org.Hs.eg.db::org.Hs.eg.db,keyType="ENTREZID")

  p=enrichplot::dotplot(gres,x='NES')

  p=p+ggplot2::scale_y_discrete(labels=function(x) stringr::str_wrap(x,width=40))
  # next, lets make the x axis label more descriptive
  #you can use xlab and ylab for this if its a ggplot
  p=p+ggplot2::xlab("Normalized Enrichment Score")
  p=p+ggplot2::theme_bw()
  rres$endot=p

  p=cnetfilt(greso,fcut)
  rres$cnetfilt=p

  ego=enrichplot::pairwise_termsim(gres)
  p=enrichplot::emapplot(ego,edge.params=list(min=.1))
  rres$emap=p
  #save figures to folder/pdf
  #only if prefix and folder
  #cant be blank
  rres$gres=gres
  rres$fc=fc
  if(!is.null(fprefix) & !is.null(folder)){
    fname=glue::glue("./{folder}/{fprefix}_res05.csv")
    readr::write_csv(rres$res05,fname)
    fname=glue::glue("./{folder}/{fprefix}_gres.csv")
    readr::write_csv(as.data.frame(rres$gres),fname)
    fname=glue::glue("./{folder}/{fprefix}_enrichdot.pdf")
    ggplot2::ggsave(fname,plot=rres$endot,width=7,height=7)
    fname=glue::glue("./{folder}/{fprefix}_cnetfilt.pdf")
    ggplot2::ggsave(fname,plot=rres$cnetfilt,width=7,height=7)
    fname=glue::glue("./{folder}/{fprefix}_emap.pdf")
    ggplot2::ggsave(fname,plot=rres$emap,width=7,height=7)
    fname=glue::glue("./{folder}/{fprefix}_pca.pdf")
    ggplot2::ggsave(fname,plot=rres$pca,width=7,height=7)
    fname=glue::glue("./{folder}/{fprefix}_disp.pdf")
    ggplot2::ggsave(fname,plot=rres$disp,width=7,height=7)
    fname=glue::glue("./{folder}/{fprefix}_volc.pdf")
    ggplot2::ggsave(fname,plot=rres$volc,width=7,height=7)
    fname=glue::glue("./{folder}/{fprefix}_topgdot.pdf")
    ggplot2::ggsave(fname,plot=rres$topgdot,width=7,height=7)
    fname=glue::glue("./{folder}/{fprefix}_heat.pdf")
    pdf(fname,width=7,height=12)
    grid::grid.newpage()
    grid::grid.draw(rres$heatmap$gtable)
    dev.off()
    fname=glue::glue("./{folder}/{fprefix}_dispest.pdf")
    pdf(fname,width=7,height=7)
    DESeq2::plotDispEsts(dds)
    dev.off()
  }

  return(rres)

}


