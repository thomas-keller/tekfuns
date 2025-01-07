#' convert limma results to dataframe with DESeq column names
#'
#' @param res_de limma toptable results
#' @param FDR false discovery rate cutoff
#'
#' @return dataframe with deseq-like result column names
#' @export
#'
#' @examples TODO
limma2df<-function(res_de,FDR=NULL){
  res<- cbind(rownames(res_de),res_de)
  names(res)[c(1,2,6)]=c('id','log2FoldChange','padj')
  res$id=as.character(res$id)
  res=res[order(res$padj),]
  if(!is.null(FDR)){
    res=res[!is.na(res$padj) & res$padj<=FDR,]
  }
  return(res)
}


#' prepare gsva results (after limma DE) for downstream genetonic
#'
#' @param obj limma toptable with gsva as input
#' @param res_de DE dataframe with ensembl ids
#' @param m_df annotations from msigdbr
#' @param gset genesets from gsva results (geneSets(gsva_res))
#' @param anno_df dataframe with ensembl to symbol mapping
#'
#' @return dataframe with standardized names for downstream genetonic
#' @export
#'
#' @examples TODO
shake_gsvaResult<-function(obj,res_de,m_df,gset,anno_df){
  m_dfu=m_df[!duplicated(m_df$gs_name),c("gs_name","gs_id","gs_description")]
  m_dfu=m_dfu[match(rownames(obj),m_dfu$gs_name),]
  anno_df2=anno_df[match(res_de$id,anno_df$gene_id),]

  gset2=gset[match(m_dfu$gs_name,names(gset))]
  gcomb=sapply(gset2,function(x) paste(anno_df2$gene_name[anno_df2$gene_name %in% x],collapse=','))
  gbg=sapply(gset2,function(x) length(x))
  gsde=sapply(gset2,function(x) sum(anno_df2$gene_name %in% x))


  mydf <- data.frame(
    gs_id = m_dfu$gs_id,
    gs_description = m_dfu$gs_name,
    gs_fulldesc=m_dfu$gs_description,
    gs_pvalue = obj$P.Value,
    gs_genes = gcomb,
    gs_de_count = gsde,
    gs_bg_count = gbg,
    gs_logFC = obj$logFC,
    gs_p.adjust = obj$adj.P.Val,
    stringsAsFactors = FALSE
  )

  rownames(mydf) <- mydf$gs_id

  return(mydf)

}

#' create a graph object based on msig annotations
#'
#'
#' Construct a gene-geneset-graph from the results of a functional enrichment
#' analysis using MSIGdb
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param res_de A `DESeqResults` object.
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be included
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be included in addition to
#' the top ones (via `n_gs`)
#' @param prettify Logical, controlling the aspect of the returned graph object.
#' If TRUE (default value), different shapes of the nodes are returned, based on
#' the node type
#' @param geneset_graph_color Character value, specifying which color should be
#' used for the fill of the shapes related to the gene sets.
#' @param genes_graph_colpal A vector of colors, also provided with their hex
#' string, to be used as a palette for coloring the gene nodes. If unspecified,
#' defaults to a color ramp palette interpolating from blue through yellow to red.
#'
#' @return An `igraph` object to be further manipulated or processed/plotted (e.g.
#' via [igraph::plot.igraph()] or
#' [visNetwork::visIgraph()][visNetwork::visNetwork-igraph])
#' @export
#'
#' @examples TODO
ggs_graphm <- function(res_enrich,
                       res_de,
                       annotation_obj = NULL,
                       gtl = NULL,
                       n_gs = 15,
                       gs_ids = NULL,
                       prettify = TRUE,
                       geneset_graph_color = "gold",
                       genes_graph_colpal = NULL) {
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }

  stopifnot(is.numeric(n_gs))
  stopifnot(is.logical(prettify))

  if (!is.null(genes_graph_colpal)) {
    if (!is(genes_graph_colpal, "character")) {
      stop(
        "Please check that you are correctly providing the color palette, ",
        "it should be encoded as a vector of colors specified as characters ",
        "(textual or hex codes)"
      )
    }

    if (!all(check_colors(genes_graph_colpal))) {
      stop(
        "You are providing your color palette in a format which ",
        "\ncan not be handled by `grDevices::col2rgb`. \n\n",
        "Try running `check_colors` on the palette object."
      )
    }
  }

  n_gs <- min(n_gs, nrow(res_enrich))

  enriched_gsids <- res_enrich[["gs_id"]]
  enriched_gsnames <- res_enrich[["gs_description"]]
  enriched_gsdescs <- res_enrich[['gs_fulldesc']]

  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(n_gs)], # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id] # the ones specified from the custom list
    )
  )

  enrich2list <- lapply(gs_to_use, function(gs) {
    go_genes <- res_enrich[gs, "gs_genes"]
    go_genes <- unlist(strsplit(go_genes, ","))
    return(go_genes)
  })
  names(enrich2list) <- res_enrich[gs_to_use, "gs_description"]

  list2df <- lapply(seq_along(enrich2list), function(gs) {
    data.frame(
      gsid = rep(names(enrich2list[gs]), length(enrich2list[[gs]])),
      gene = enrich2list[[gs]]
    )
  })
  list2df <- do.call("rbind", list2df)

  g <- igraph::graph_from_data_frame(list2df, directed = FALSE)

  nodeIDs_gs <- which(names(igraph::V(g)) %in% enriched_gsnames)
  nodeIDs_genes <- which(!(names(igraph::V(g)) %in% enriched_gsnames))

  igraph::V(g)$nodetype <- NA
  igraph::V(g)$nodetype[nodeIDs_gs] <- "GeneSet"
  igraph::V(g)$nodetype[nodeIDs_genes] <- "Feature"

  if (prettify) {
    # different shapes based on the node type

    # this does not work with visNetwork?
    # igraph::V(g)$value <- 15 # size? size2? or does this not work with the shapes I selected?
    # igraph::V(g)$value[nodeIDs_gs] <- 45

    # this one is handled correctly by visNetwork
    igraph::V(g)$shape <- c("box", "ellipse")[factor(igraph::V(g)$nodetype, levels = c("GeneSet", "Feature"))]


    # different colors for the gene nodes in function of their logFC
    fcs_genes <- res_de[annotation_obj$gene_id[match((igraph::V(g)$name[nodeIDs_genes]), annotation_obj$gene_name)], ]$log2FoldChange

    if (!is.null(genes_graph_colpal)) {
      mypal <- genes_graph_colpal
    } else {
      mypal <- rev(scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.4
      ))
    }

    igraph::V(g)$color[nodeIDs_genes] <- mosdef::map_to_color(fcs_genes, mypal, limits = c(-4, 4))
    igraph::V(g)$color[nodeIDs_gs] <- geneset_graph_color

    # title for tooltips
    igraph::V(g)$title <- NA
    igraph::V(g)$title[nodeIDs_gs] <- paste0(
      "<h4>",
      sprintf('<a href="https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/%s.html" target="_blank">%s</a>', enriched_gsnames[nodeIDs_gs], enriched_gsids[nodeIDs_gs]), "</h4><br>",
      igraph::V(g)$name[nodeIDs_gs], "<br><br>",
      sapply(enriched_gsdescs[nodeIDs_gs],function(x) paste0(strwrap(x,50),collapse='<br>'))
    )
    igraph::V(g)$title[nodeIDs_genes] <- paste0(
      "<h4>", igraph::V(g)$name[nodeIDs_genes], "</h4><br>",
      "logFC = ", format(round(fcs_genes, 2), nsmall = 2)
    )

  } else {
    igraph::V(g)$color[nodeIDs_genes] <- "#B3B3B3"
    igraph::V(g)$color[nodeIDs_gs] <- "#E5C494"
  }

  # re-sorting the vertices alphabetically
  rank_gs <- rank(igraph::V(g)$name[igraph::V(g)$nodetype == "GeneSet"])
  rank_feats <- rank(igraph::V(g)$name[igraph::V(g)$nodetype == "Feature"]) +
    length(rank_gs) # to keep the GeneSets first
  g <- igraph::permute(g, c(rank_gs, rank_feats))

  return(g)
}

#' Visually enhances a functional enrichment result table
#'
#' Creates a visual summary for the results of a functional enrichment analysis,
#' by displaying also the components of each gene set and their expression change
#' in the contrast of interest
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param res_de  A `DESeqResults` object.
#' @param annotation_obj A `data.frame` object with the feature annotation.
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be displayed.
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be displayed.
#' @param chars_limit Integer, number of characters to be displayed for each
#' geneset name.
#' @param plot_style Character value, one of "point" or "ridgeline". Defines the
#' style of the plot to summarize visually the table.
#' @param ridge_color Character value, one of "gs_id" or "gs_score", controls the
#' fill color of the ridge lines. If selecting "gs_score", the `z_score` column
#' must be present in the enrichment results table - see `get_aggrscores()` to do
#' that.
#' @param plot_title Character string, used as title for the plot. If left `NULL`,
#' it defaults to a general description of the plot and of the DE contrast.
#'
#' @return A `ggplot` object
#' @export
#'
#' @examples
#'
#' library("macrophage")
#' library("DESeq2")
#' library("org.Hs.eg.db")
#' library("AnnotationDbi")
#'
#' # dds object
#' data("gse", package = "macrophage")
#' dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage <- estimateSizeFactors(dds_macrophage)
#'
#' # annotation object
#' anno_df <- data.frame(
#'   gene_id = rownames(dds_macrophage),
#'   gene_name = mapIds(org.Hs.eg.db,
#'     keys = rownames(dds_macrophage),
#'     column = "SYMBOL",
#'     keytype = "ENSEMBL"
#'   ),
#'   stringsAsFactors = FALSE,
#'   row.names = rownames(dds_macrophage)
#' )
#'
#' # res object
#' data(res_de_macrophage, package = "GeneTonic")
#' res_de <- res_macrophage_IFNg_vs_naive
#'
#' # res_enrich object
#' data(res_enrich_macrophage, package = "GeneTonic")
#' res_enrich <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
#' res_enrich <- get_aggrscores(res_enrich, res_de, anno_df)
#' enhance_table(res_enrich,
#'   res_de,
#'   anno_df,
#'   n_gs = 10
#' )
#'
#' # using the ridge line as a style, also coloring by the Z score
#' res_enrich_withscores <- get_aggrscores(
#'   res_enrich,
#'   res_de,
#'   anno_df
#' )
#' enhance_table(res_enrich_withscores,
#'   res_de,
#'   anno_df,
#'   n_gs = 10,
#'   plot_style = "ridgeline",
#'   ridge_color = "gs_score"
#' )
enhance_tablem <- function(res_enrich,
                          res_de,
                          annotation_obj,
                          gtl = NULL,
                          n_gs = 50,
                          gs_ids = NULL,
                          chars_limit = 70,
                          plot_style = c("point", "ridgeline"),
                          ridge_color = c("gs_id", "gs_score"),
                          plot_title = NULL) {
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }

  plot_style <- match.arg(plot_style, c("point", "ridgeline"))
  ridge_color <- match.arg(ridge_color, c("gs_id", "gs_score"))

  n_gs <- min(n_gs, nrow(res_enrich))

  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(n_gs)], # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id] # the ones specified from the custom list
    )
  )

  gs_fulllist <- lapply(gs_to_use, function(gs) {
    genes_thisset <- res_enrich[gs, "gs_genes"]
    genes_thisset <- unlist(strsplit(genes_thisset, ","))

    genesid_thisset <- annotation_obj$gene_id[match(genes_thisset, annotation_obj$gene_name)]

    # removing the genes not finding a match in the annotation
    no_anno_match <- is.na(genesid_thisset)
    genes_thisset_anno <- genes_thisset[!no_anno_match]
    genesid_thisset_anno <- genesid_thisset[!no_anno_match]
    # ... and informing on which genes might be troublesome
    if (any(no_anno_match)) {
      message("Could not find a match in the annotation for some genes. ",
              "Please inspect your results in detail for geneset ",
              gs,
              " the gene(s) named: ",
              paste0(genes_thisset[no_anno_match], collapse = ", "))
    }

    res_thissubset <- res_de[genesid_thisset_anno, ]

    res_thissubset <- as.data.frame(res_thissubset)

    res_thissubset$gene_name <- genes_thisset_anno
    res_thissubset$gs_desc <- as.factor(res_enrich[gs, "gs_description"])
    res_thissubset$gs_id <- res_enrich[gs, "gs_id"]
    # return(as.data.frame(res_thissubset))
    return(res_thissubset)
  })
  gs_fulllist <- do.call(rbind, gs_fulllist)
  # message(dim(gs_fulllist)[1])

  #contrast/mcols doesnt have an analog for toptable
  #drop the title for now
  #probably easiest to just add as a manual title from selecting coef in limma results
  #this_contrast <- (sub(".*p-value: (.*)", "\\1", mcols(res_de, use.names = TRUE)["pvalue", "description"]))

  # to have first rows viewed on top
  gs_fulllist <- gs_fulllist[rev(seq_len(nrow(gs_fulllist))), ]
  gs_fulllist$gs_desc <- factor(gs_fulllist$gs_desc, levels = rev(levels(gs_fulllist$gs_desc)))
  max_lfc <- max(abs(range(gs_fulllist$log2FoldChange)))

  # common elements here
  gs_labels <- paste0(
    substr(as.character(unique(gs_fulllist$gs_desc)), 1, chars_limit),
    " | ", unique(gs_fulllist$gs_id)
  )

  if (plot_style == "point") {
    p <- ggplot(
      gs_fulllist, aes(
        x = .data$log2FoldChange,
        y = .data$gs_desc,
        fill = .data$gs_id,
        text = .data$gene_name
      )
    ) +
      scale_x_continuous(limits = c(-max_lfc, max_lfc)) +
      geom_point(alpha = 0.7, shape = 21, size = 2) +
      theme_minimal() +
      geom_vline(aes(xintercept = 0), col = "steelblue", alpha = 0.4) +
      theme(legend.position = "none") +
      scale_y_discrete(
        name = "",
        labels = gs_labels
      ) +
      labs(x = "log2 Fold Change")

  } else if (plot_style == "ridgeline") {

    if (ridge_color == "gs_score" & is.null(res_enrich$z_score)) {
      message("Fallback to plotting the ridgelines according to geneset id (Z score required)")
      ridge_color <- "gs_id"
    }

    if (ridge_color == "gs_score") {
      gs_fulllist$gs_zscore <- res_enrich$z_score[match(gs_fulllist$gs_id, res_enrich$gs_id)]
      p <- ggplot(
        gs_fulllist, aes(
          x = .data$log2FoldChange,
          y = .data$gs_desc,
          fill = .data$gs_zscore
        )
      ) +
        scale_x_continuous(limits = c(-max_lfc, max_lfc)) +
        scale_fill_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026") +
        ggridges::geom_density_ridges(
          aes(group = .data$gs_id),
          point_color = "#00000066",
          jittered_points = TRUE, scale = .95, rel_min_height = .01,
          point_shape = "|", point_size = 3, linewidth = 0.25,
          position = ggridges::position_points_jitter(height = 0)) +
        theme_minimal() +
        geom_vline(aes(xintercept = 0), col = "steelblue", alpha = 0.4) +
        scale_y_discrete(
          name = "",
          labels = gs_labels
        ) +
        labs(x = "log2 Fold Change")
    }
    else if (ridge_color == "gs_id") {
      p <- ggplot(
        gs_fulllist, aes(
          x = .data$log2FoldChange,
          y = .data$gs_desc,
          fill = .data$gs_id
        )
      ) +
        scale_x_continuous(limits = c(-max_lfc, max_lfc)) +
        ggridges::geom_density_ridges(
          aes(group = .data$gs_id),
          point_color = "#00000066",
          jittered_points = TRUE, scale = .95, rel_min_height = .01,
          point_shape = "|", point_size = 3, linewidth = 0.25,
          position = ggridges::position_points_jitter(height = 0)) +
        theme_minimal() +
        geom_vline(aes(xintercept = 0), col = "steelblue", alpha = 0.4) +
        theme(legend.position = "none") +
        scale_y_discrete(
          name = "",
          labels = gs_labels
        ) +
        labs(x = "log2 Fold Change")
    }
  }

  if (is.null(plot_title)) {
    #just drop the contrast in title for now since no analog in limma
    p <- p + ggtitle(paste0("Enrichment overview "))
    #p <- p + ggtitle(paste0("Enrichment overview - ", this_contrast))
  } else {
    p <- p + ggtitle(plot_title)
  }

  return(p)
}


