



#' Filter a gene set enrichment result to eliminate genes with low absolute fold-change
#' @import org.Hs.eg.db
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
#' for transcript level analysis, need a tximeta object with bootstraps/ gibbs samples
#'
#' @param dds DESeqDataset object
#' @param sw tximeta object from salmon w/ inferential reps
#' @param regulons infer regulons with GENIE3 method? default FALSE
#' @param pcut gsea pvalue threshold
#' @param folder output result folder
#' @param fprefix output analysis prefix
#' @param nfcut gsea filter for network plots
#'
#' @return list of plots
#' @export
#'
#' @examples todo
rnaplots <- function(dds,sw=NULL,regulons=FALSE,pcut=0.05,nfcut=2,fcut=.5,folder=NULL,fprefix=NULL){
  if(!is.null(folder)){
    dir_create(folder)
  }
  rres=list()

  if(!is.null(sw)){
    #add fishpond transcript level analysis
    #input is tximeta object
    sw <- fishpond::scaleInfReps(sw) # scales counts
    sw <- fishpond::labelKeep(sw) # labels features to keep
    sw <- sw[SummarizedExperiment::mcols(sw)$keep,]
    set.seed(1)
    #still assuming metadata has a condition column of main interest
    sw <- fishpond::swish(sw, x="condition") # simplest Swish case
    cols <- c("log10mean","log2FC","pvalue","qvalue")
    most.sig <-SummarizedExperiment::mcols(sw) %>%
      as.data.frame() %>%
      dplyr::arrange(qvalue,desc(abs(log2FC)))
    sig <- SummarizedExperiment::mcols(sw)$qvalue < .05
    lo <- order(SummarizedExperiment::mcols(sw)$log2FC * sig)
    hi <- order(-SummarizedExperiment::mcols(sw)$log2FC * sig)
    #todo implement a ggplot for infrep plot
    #could move to a function
    sw <- tximeta::addIds(sw, "SYMBOL", gene=TRUE)
    #get top/bottom most genes by lfc
    gids=c(hi[1:5],lo[1:5])
    infReps <- SummarizedExperiment::assays(sw[gids,])[grep("infRep",
                                                           SummarizedExperiment::assayNames(sw))]
    infReps<- unlist(infReps)
    infReps=as.data.frame(infReps)
    infRepn=stringr::str_split_i(rownames(infReps),"\\.",1)
    infRept=stringr::str_split_i(rownames(infReps),"\\.",2)
    infReps$rep=infRepn
    infReps$tname=infRept
    infl=infReps %>% tidyr::pivot_longer(-all_of(c('rep','tname')),
                                         names_to='sample',
                                         values_to='expression')
    #infl$condition=stringr::str_split_i(infl$sample,"_",2)
    md=mcols(sw)[,c("tx_name","SYMBOL")]
    md=as.data.frame(md)
    cold=as.data.frame(colData(sw))
    infl=dplyr::inner_join(infl,md,by=dplyr::join_by(tname==tx_name))
    infl=dplyr::inner_join(infl,cold,by=dplyr::join_by(sample==names))
    p=ggplot2::ggplot(infl,aes(x=SYMBOL,
                               y=expression,
                               color=condition,
                               shape=sample))+scale_y_log10()+geom_jitter()
    p=p+ggpubr::theme_pubr()
    p=p+ggpubr::labs_pubr()
    rres$toptscript=p
    #facet wrap the top 5-10 genes/transcripts?
    # do gene  level, isoform test
    #add results to output

  } else{
    rres$toptscript=NULL
  }



  rname=resultsNames(dds)
  lrn=rname[length(rname)]
  #filter out genes with low reads
  smallestGroupSize <- 3
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]

  #run deseq method, it will automatically perform several steps
  dds<-DESeq2::DESeq(dds)

  #extract the results
  reso <- DESeq2::results(dds,alpha=.05)
  #will automatically use last coef in regression
  rname=resultsNames(dds)
  lrn=rname[length(rname)]
  res=DESeq2::lfcShrink(dds,res=reso,coef=lrn)
  res$unslfc=reso$log2FoldChange
  #do some cleaning and filter out nonsignificant genes
  res05 <- na.omit(res)
  res05 <-as.data.frame(res05)
  res05 <- dplyr::filter(res05,padj<=0.05) %>% dplyr::filter(abs(log2FoldChange)>=fcut)
  #res05 <- dplyr::arrange(res05,desc(abs(log2FoldChange)))
  res05 <- dplyr::arrange(res05,padj)
  res05$ens=rownames(res05)
  res05$lprank=res05$log2FoldChange*res05$pvalue
  ens=res05$ens
  res05$symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ens, column = c('SYMBOL'), keytype = 'ENSEMBL')
  res05$entrez <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ens,column = c('ENTREZID'),
                         keytype = 'ENSEMBL')
  res05$csymbol=ifelse(is.na(res05$symbol),res05$ens,res05$symbol)
  rres$res05=res05
  vsd=DESeq2::vst(dds,blind=FALSE)
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
                                   legendLabels=c('Not sig.','Log2 FC','p-value','p-value & Log2 FC'),
                                   legendPosition = 'bottom',
                                   legendLabSize = 10,
                                   legendIconSize = 5.0)
  rres$volc=p

  #heatmap
  #res05 arranged by padj now
  topg <- res05$ens[1:30]
  mat  <- SummarizedExperiment::assay(vsd)
  mat<-mat[row.names(mat) %in% res05$ens,]
  #mat  <- mat[row.names(mat) %in% topg, ]
  #do z scale
  #mat=mat[row.names(mat) %in% topg,]
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
  #mat=mat[1:30,]
  #row.names(mat)=res05$csymbol[1:30]
  myBreaks <- c(seq(min(mat), 0, length.out=ceiling(plen/2) + 1),
                seq(max(mat)/plen, max(mat), length.out=floor(plen/2)))

  #start with 4 clusters
  #cluster the whole expression
  numc=2
  cl=cutree(hclust(dist(mat)),numc)
  canno=data.frame(cluster=paste("Cluster",cl,sep=" "))
  rownames(canno)<-rownames(mat)
  fcl=list()
  enl=list()
  for(i in 1:numc){
    df=res05[cl==i,]
    df = df %>% arrange(pvalue)
    cname=paste('Cluster',i)
    fc=df$log2FoldChange
    names(fc)=df$entrez
    fc=sort(fc,decreasing=T)
    fc=fc[!is.na(names(fc))]
    ens=df$ens
    #print(ens)
    #print(fc)
    fcl[[cname]]=fc
    enl[[cname]]=ens


  }




  #p=pheatmap::pheatmap(mat,annotation_col=anno,annotation_row=canno,treeheight_row=0,treeheight_col=0,silent=TRUE,breaks=myBreaks,color=cols,show_rownames=F)
  #p=pheatmap::pheatmap(mat,annotation_col=anno,annotation_row=canno,treeheight_row=0,treeheight_col=0,breaks=myBreaks,color=cols,show_rownames=F)
  ca=ComplexHeatmap::HeatmapAnnotation(condition=anno$condition)
  ra=ComplexHeatmap::rowAnnotation(cluster=canno$cluster)
  p=ComplexHeatmap::Heatmap(mat,top_annotation=ca,
                            left_annotation=ra,
                            show_row_names=F,
                            heatmap_legend_param=list(title='Expression',
                                                      title_position = "leftcenter-rot"))
  ComplexHeatmap::draw(p)

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
  p=p+ggpubr::theme_pubr()
  #convert y axis to log scale since values range over orders of magnitude
  p=p+ggplot2::scale_y_log10()
  #rename x and y labels
  p=p+ggplot2::xlab("gene symbol")+ggplot2::ylab("log10 counts")
  rres$topgdot=p
  #now to enrich plots
  fc=res05$log2FoldChange
  names(fc)=res05$entrez
  fc=sort(fc,decreasing=TRUE)
  fc=fc[!is.na(names(fc))]
  fc=fc[!duplicated(names(fc))]




  set.seed(1234)

  m_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene)
  m_t2g$gs_name=stringr::str_remove(m_t2g$gs_name,"HALLMARK_")
  m_t2g$gs_name=stringr::str_remove(m_t2g$gs_name,"HALLMARK ")
  #replace the underscores with spaces to make wrapping easier
  m_t2g$gs_name=stringr::str_replace_all(m_t2g$gs_name,"_"," ")

  #try c2 as well for more fine grained (also potentially redundant) curated pathways
  m2_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2") %>%
    dplyr::select(gs_name, entrez_gene)
  m2_t2g$gs_name=stringr::str_replace_all(m2_t2g$gs_name,"_"," ")
  #GSEA function takes minimum two arguments; the fold-change values, and the pathway-gene database
  #if you take a look at the m_t2g dataframe we constructed from msigdbr, its has two columns, one with the pathway name (gs_name)
  #and one with the gene id (entrez_gene)

  gres=clusterProfiler::GSEA(fc,TERM2GENE=m_t2g,pvalueCutoff=pcut)

  gres2=clusterProfiler::GSEA(fc,TERM2GENE=m2_t2g,pvalueCutoff=pcut)
  #only use original enrich result for filtering genes with cnetfilt
  greso=gres
  #use setReadable to convert entrez ids to gene names


  if(nrow(gres)>0){
  gres=clusterProfiler::setReadable(gres,OrgDb=org.Hs.eg.db::org.Hs.eg.db,keyType="ENTREZID")
  p=enrichplot::dotplot(gres,x='NES')

  p=p+ggplot2::scale_y_discrete(labels=function(x) stringr::str_wrap(x,width=40))
  # next, lets make the x axis label more descriptive
  #you can use xlab and ylab for this if its a ggplot
  p=p+ggplot2::xlab("Normalized Enrichment Score")
  p=p+ggpubr::theme_pubr()
  rres$endoth=p

  p=cnetfilt(greso,nfcut)
  rres$cnetfilt=p

  ego=enrichplot::pairwise_termsim(gres)
  p=enrichplot::emapplot(ego,edge.params=list(min=.1))
  rres$emap=p
  }

  if(nrow(gres2@result)>0){
  gres2=clusterProfiler::setReadable(gres2,OrgDb=org.Hs.eg.db::org.Hs.eg.db,keyType="ENTREZID")
  p=enrichplot::dotplot(gres2,x='NES')
  p=p+scale_fill_continuous_diverging(mid=0)

  p=p+ggplot2::scale_y_discrete(labels=function(x) stringr::str_wrap(x,width=40))
  # next, lets make the x axis label more descriptive
  #you can use xlab and ylab for this if its a ggplot
  p=p+ggplot2::xlab("Normalized Enrichment Score")
  p=p+ggpubr::theme_pubr()
  rres$endotc2=p



  }
  #save figures to folder/pdf
  #only if prefix and folder
  #cant be blank
  rres$gres=gres
  rres$fc=fc

  rres$fcl=fcl
  rres$enl=enl
  cores=compareCluster(geneCluster=fcl,fun="gseGO",OrgDb=org.Hs.eg.db::org.Hs.eg.db,nPermSimple=10000,seed=42,pvalueCutoff=.1)

  #something overwrote the enrichplot dotplot, gives error unless prefix
  if(!is.null(cores)){
    cores=simplify(cores)
    p=enrichplot::dotplot(cores,color='NES')+scale_x_discrete(guide=guide_axis(n.dodge=2))
    p=p+scale_fill_continuous_diverging(mid=0)
    #kind of hacky way to rename the facets from activated/supresed to up/down reg
    #p$data$.sign=ifelse(p$data$.sign=='activated','Up-reg.','Down-reg.')
    #p$data$.sign=factor(p$data$.sign,c("Up-reg.","Down-reg."),ordered=T)
    p=p+xlab("Cluster")
    p=p+ggpubr::theme_pubr()
    rres$compbp=p
  } else{
    rres$compbp=NULL
  }
  cores=compareCluster(geneCluster=fcl,fun="GSEA",TERM2GENE=m_t2g,nPermSimple=10000,seed=42,pvalueCutoff=.1)
  #something overwrote the enrichplot dotplot, gives error unless prefix
  if(!is.null(cores)){
  p=enrichplot::dotplot(cores,color='NES')+scale_x_discrete(guide=guide_axis(n.dodge=2))
  p=p+scale_fill_continuous_diverging(mid=0)
  #kind of hacky way to rename the facets from activated/supresed to up/down reg
  #p$data$.sign=ifelse(p$data$.sign=='activated','Up-reg.','Down-reg.')
  #p$data$.sign=factor(p$data$.sign,c("Up-reg.","Down-reg."),ordered=T)
  p=p+xlab("Cluster")
  p=p+ggpubr::theme_pubr()
  rres$comph=p
  rres$comphst=as.data.frame(cores)

  } else{
    rres$comph=NULL
    rres$comphst=NULL
  }
  cores=compareCluster(geneCluster=fcl,fun="GSEA",TERM2GENE=m2_t2g,nPermSimple=10000,seed=42,pvalueCutoff=.1)
  #something overwrote the enrichplot dotplot, gives error unless prefix
  if(!is.null(cores)){
  p=enrichplot::dotplot(cores,color='NES')+scale_x_discrete(guide=guide_axis(n.dodge=2))
  p=p+colorspace::scale_fill_continuous_diverging(mid=0)
  #kind of hacky way to rename the facets from activated/supresed to up/down reg
  #p$data$.sign=ifelse(p$data$.sign=='activated','Up-reg.','Down-reg.')
  #p$data$.sign=factor(p$data$.sign,c("Up-reg.","Down-reg."),ordered=T)
  p=p+xlab("Cluster")
  p=p+ggpubr::theme_pubr()
  rres$compc2=p
  rres$compc2st

  } else{
    rres$compc2=NULL
    rres$compc2st=NULL
  }
  gpres=gprofiler2::gost(query=enl,organism='hsapiens',ordered_query=T,multi_query=TRUE,source=c("GO",'KEGG','REAC','CORUM'))
  rres$gpres=gpres
  hl=c()
  gptab=data.frame()
  gr=gpres$result %>% dplyr::filter(source=='GO:BP') %>% dplyr::pull(term_id)
  if(length(gr)>0){
    hl=c(hl,gr[1:2])
    r=gpres$result %>% dplyr::filter(source=='GO:BP') %>% dplyr::select(source,term_id,term_name,p_values)
    gptab=rbind(gptab,r[1:3,])
  }
  gr=gpres$result %>% dplyr::filter(source=='REAC') %>% dplyr::pull(term_id)
  if(length(gr)>0){
    hl=c(hl,gr[1:2])
    r=gpres$result %>% dplyr::filter(source=='REAC') %>% dplyr::select(source,term_id,term_name,p_values)
    gptab=rbind(gptab,r[1:3,])
  }
  gptab$term_name=stringr::str_wrap(gptab$term_name,40)
  rres$gptab=gptab
  pi=gprofiler2::gostplot(gpres, capped = TRUE, interactive = TRUE)
  rres$gp=pi
  #p=gostplot(gostres, capped = TRUE, interactive = FALSE)
  #rres$gp2=p
  p=gprofiler2::gostplot(gpres, capped = TRUE, interactive = FALSE)
  #pp <- gprofiler2::publish_gostplot(p, highlight_terms = hl,
  #                       width = NA, height = NA, filename = NULL )
  pp<-gprofiler2::publish_gostplot(p)
  rres$gp2=pp

  #get regulons w/ genie3 if regulons not false
  if(regulons){
    set.seed(54321)
    #could expand to handle candidate regulators
    wm<-GENIE3::GENIE3(mat,nCores=8)
    ll=GENIE3::getLinkList(wm)
    rres$ll=ll
  } else{
    rres$ll=NULL
  }

  if(!is.null(fprefix) & !is.null(folder)){
    fname=glue::glue("./{folder}/{fprefix}_res05.csv")
    readr::write_csv(rres$res05,fname)
    fname=glue::glue("./{folder}/{fprefix}_gresH.csv")
    readr::write_csv(as.data.frame(rres$gres),fname)
    fname=glue::glue("./{folder}/{fprefix}_compH.csv")
    readr::write_csv(as.data.frame(rres$comphst),fname)
    fname=glue::glue("./{folder}/{fprefix}_gresC2.csv")
    readr::write_csv(as.data.frame(rres$gres2),fname)
    fname=glue::glue("./{folder}/{fprefix}_compc2.csv")
    readr::write_csv(as.data.frame(rres$compc2st),fname)
    fname=glue::glue('./{folder/{fprefix}_reglist.csv')
    readr::write_csv(rres$ll,fname)
    fname=glue::glue("./{folder}/{fprefix}_toptscript.pdf")
    ggplot2::ggsave(fname,plot=rres$toptscript,width=7,height=10)
    fname=glue::glue("./{folder}/{fprefix}_gprof.pdf")
    ggplot2::ggsave(fname,plot=rres$g2,width=7,height=10)
    fname=glue::glue("./{folder}/{fprefix}_comph.pdf")
    ggplot2::ggsave(fname,plot=rres$comph,width=7,height=10)
    fname=glue::glue("./{folder}/{fprefix}_compc2.pdf")
    ggplot2::ggsave(fname,plot=rres$compc2,width=7,height=10)
    fname=glue::glue("./{folder}/{fprefix}_gprof.pdf")
    ggplot2::ggsave(fname,plot=rres$gp2,width=10,height=7)
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
    ComplexHeatmap::draw(rres$heatmap)
    dev.off()
    fname=glue::glue("./{folder}/{fprefix}_dispest.pdf")
    pdf(fname,width=7,height=7)
    DESeq2::plotDispEsts(dds)
    dev.off()

    fname=glue::glue("./{folder}/{fprefix}_toptscript.png")
    ggplot2::ggsave(fname,plot=rres$toptscript,width=7,height=7,units='in',dpi=600)
    fname=glue::glue("./{folder}/{fprefix}_H_enrichdot.png")
    ggplot2::ggsave(fname,plot=rres$endoth,width=7,height=7,units='in',dpi=600)
    fname=glue::glue("./{folder}/{fprefix}_C2_enrichdot.png")
    ggplot2::ggsave(fname,plot=rres$endotc2,width=7,height=7,units='in',dpi=600)
    fname=glue::glue("./{folder}/{fprefix}_cnetfilt.png")
    ggplot2::ggsave(fname,plot=rres$cnetfilt,width=7,height=7,units='in',dpi=600)
    fname=glue::glue("./{folder}/{fprefix}_emap.png")
    ggplot2::ggsave(fname,plot=rres$emap,width=7,height=7,units='in',dpi=600)
    fname=glue::glue("./{folder}/{fprefix}_pca.png")
    ggplot2::ggsave(fname,plot=rres$pca,width=7,height=7,units='in',dpi=600)
    fname=glue::glue("./{folder}/{fprefix}_disp.png")
    ggplot2::ggsave(fname,plot=rres$disp,width=7,height=7,units='in',dpi=600)
    fname=glue::glue("./{folder}/{fprefix}_volc.png")
    ggplot2::ggsave(fname,plot=rres$volc,width=7,height=7,units='in',dpi=600)
    fname=glue::glue("./{folder}/{fprefix}_topgdot.png")
    ggplot2::ggsave(fname,plot=rres$topgdot,width=7,height=7,units='in',dpi=600)
    fname=glue::glue("./{folder}/{fprefix}_heat.png")
    png(fname,width=7,height=12,units='in',res=600)
    ComplexHeatmap::draw(rres$heatmap)
    dev.off()
    fname=glue::glue("./{folder}/{fprefix}_dispest.pdf")
    png(fname,width=7,height=10,units='in',res=600)
    DESeq2::plotDispEsts(dds)
    dev.off()
  }

  return(rres)

}


