

convert_to_seurat <- function(adata) {

  # 1. let's get the adata object layers  
  counts <- t(Matrix::Matrix(py_to_r(adata$layers[["counts"]]), sparse = TRUE))
  
 # 2. let's get the metadata
  meta <- as.data.frame(py_to_r(adata$obs))
  cell_names <- rownames(meta)
  colnames(counts) <- cell_names
  
  
  # 3. now we create the seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = meta
  )
  
  # 4. let's add embeddings from obsm 
  for (key in adata$obsm_keys()) {
    mat <- py_to_r(adata$obsm[[key]])
    mat <- as.matrix(mat)
    
    # ensure rownames match cells
    rownames(mat) <- cell_names
    
    # Clean key name
    clean_name <- gsub("^X_|_", "", key)
    
    if (grepl("umap", tolower(clean_name))) {
      seurat_obj[[clean_name]] <- CreateDimReducObject(
        embeddings = mat,
        key = paste0(clean_name,"_"),
        assay = DefaultAssay(seurat_obj)
      )
    } else if (grepl("pca", tolower(clean_name))) {
      seurat_obj[[clean_name]] <- CreateDimReducObject(
        embeddings = mat,
        key = paste0(clean_name,"_"),
        assay = DefaultAssay(seurat_obj)
      )
    } else {
      seurat_obj@misc[[clean_name]] <- mat
    }
  }
  
  # --- 6. Variable features ---
  var_df <- as.data.frame(py_to_r(adata$var))
  if ("highly_variable" %in% colnames(var_df)) {
    hvf <- rownames(var_df)[which(var_df$highly_variable)]
    VariableFeatures(seurat_obj, assay = "RNA") <- hvf
  }
  
  return(seurat_obj)
}


plot_one_over_other <- function(umapdf, x, y, colvar, condition, group1, group2, integration) {
  p <- ggplot() +
    geom_point(data = umapdf %>% dplyr::filter(get(condition) == group1),
               aes(x = get(x), y = get(y), color = get(colvar)),
               alpha = 1, size = 0.5) +
    geom_point(data = umapdf %>% dplyr::filter(get(condition) == group2),
               aes(x = get(x), y = get(y), color = get(colvar)),
               alpha = 1, size = 0.5) +
    theme_classic() +
    labs(x = x, y = y, color = colvar) +
    ggtitle(integration) 
  return(p)
} 


plot_general <- function(umapdf, x, y, colvar) {
     p <- umapdf %>% 
    ggplot(aes(x = get(x), y = get(y), color = get(colvar))) +
    geom_point(alpha = 0.5, size = 0.5) +
    theme_classic() +
    labs(x = x, y = y, color = colvar) 
  
  return(p)
}

plot_with_patchwork <- function(colvar, condition, group1, group2) {

  plot_one_over_other(umap_etmr, "umapunintegrated_1", "umapunintegrated_2", colvar, condition, group1, group2, "unintegrated") +
    plot_one_over_other(umap_harm, "umapharmony_1", "umapharmony_2", colvar, condition, group1, group2, "harmony") +
    plot_layout(guides = "collect", axes = "collect") &
    plot_annotation(title = paste0(colvar, " distribution across integration methods in ETMR snRNA data")) 
}

plot_hist <- function(umapdf, x, colvar) {

  h <- umapdf %>%
    ggplot(aes(x = get(x),  fill = get(colvar))) +
    geom_bar(position = "fill") +
    theme_classic() +
    ggtitle(paste(colvar)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = x, fill = colvar, y = "Percent") +
    scale_y_continuous(labels = scales::percent)
  
  return(h)
}


hist_with_patchwork <- function(colvar) {
  plot_hist(umap_etmr, "unintergated_clusters", colvar, "unintegrated") +
    plot_hist(umap_harm, "leiden_harmony",  colvar,  "harmony") +
    plot_layout(guides = "collect", axes = "collect") &
    plot_annotation(title = paste0(colvar, " histogram across integration methods in ETMR snRNA data")) 
}



run_umap <- function (adata, pc=30, red="pca", kk=20, rr=0.8, integration = "unintegrated") {

      adata <- FindNeighbors(
        adata,
        dims = 1:pc,
        reduction = red,
        k.param = kk
      )
      
      adata <- FindClusters(
        adata,
        resolution = rr,
        algorithm = 4,
        random.seed = 2025,
        cluster.name = paste0(integration, "_clusters")
      )
      
      adata <- RunUMAP(
        adata,
        reduction = red,
        dims = 1:pc,
        n.neighbors = kk,
        seed.use = 2025,
        reduction.name = paste0("umap_",integration)
      )
  
    return(adata)    
    }

c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")


plot_hist_oneval <- function(umapdf, x, colvar, fill_c) {umapdf %>%
  group_by(.data[[x]]) %>%
  summarize(p_true = mean(.data[[colvar]], na.rm = TRUE)) %>%
  ggplot(aes(x = .data[[x]], y = p_true)) +
  geom_col(fill = fill_c) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  theme_classic() +
  labs(y = "Percentage")
}


get_pc_terms <- function(loadings, dims = 1:30, nfeatures = 20) {
  
  res_list <- lapply(dims, function(pc) {
    
    pc_vals <- loadings[, pc]
    
    pos <- sort(pc_vals, decreasing = TRUE)[1:nfeatures]
    neg <- sort(pc_vals, decreasing = FALSE)[1:nfeatures]
    
    tibble(
      PC = paste0("PC", pc),
      positive = paste(names(pos), collapse = ", "),
      negative = paste(names(neg), collapse = ", ")
    )
  })
  
  bind_rows(res_list)
}

plot_cells_3d_fixed <- function(cds,
                                dims = c(1,2,3),
                                reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
                                color_cells_by="cluster",
                                #group_cells_by=c("cluster", "partition"), #
                                genes=NULL,
                                show_trajectory_graph=TRUE,
                                trajectory_graph_color="black",
                                trajectory_graph_segment_size=5,
                                norm_method = c("log", "size_only"),
                                color_palette = NULL,
                                color_scale = "Viridis",
                                cell_size=25,
                                alpha = 1,
                                min_expr=0.1) {
  
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimension with",
                                      "reduction_method =", reduction_method,
                                      "before attempting to plot."))
  low_dim_coords <- reducedDims(cds)[[reduction_method]]
  if(!is.null(color_cells_by)) {
    assertthat::assert_that(color_cells_by %in% c("cluster", "partition",
                                                  "pseudotime") |
                              color_cells_by %in% names(colData(cds)),
                            msg = paste("color_cells_by must be a column in",
                                        "the colData table."))
  }
  
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers),
                          msg = paste("Either color_cells_by or markers must",
                                      "be NULL, cannot color by both!"))
  norm_method = match.arg(norm_method)
  
  if (show_trajectory_graph &&
      is.null(principal_graph(cds)[[reduction_method]])) {
    message("No trajectory to plot. Has learn_graph() been called yet?")
    show_trajectory_graph = FALSE
  }
  
  if (class(cds@colData[,color_cells_by])=="character") {
    cds@colData[,color_cells_by] = factor(cds@colData[,color_cells_by])
  }	
  
  gene_short_name <- NA
  sample_name <- NA
  
  x <- dims[[1]]
  y <- dims[[2]]
  z <- dims[[3]]
  
  S_matrix <- reducedDims(cds)[[reduction_method]]
  data_df <- data.frame(S_matrix[,c(dims)])
  
  colnames(data_df) <- c("data_dim_1", "data_dim_2", "data_dim_3")
  data_df$sample_name <- row.names(data_df)
  
  data_df <- as.data.frame(cbind(data_df, colData(cds)))
  
  if (color_cells_by == "cluster"){
    data_df$cell_color <- tryCatch({
      clusters(cds, reduction_method = reduction_method)[data_df$sample_name]},
      error = function(e) {NULL})
  } else if (color_cells_by == "partition") {
    data_df$cell_color <- tryCatch({
      partitions(cds,
                 reduction_method = reduction_method)[data_df$sample_name]},
      error = function(e) {NULL})
  } else if (color_cells_by == "pseudotime") {
    data_df$cell_color <- tryCatch({
      pseudotime(cds,
                 reduction_method = reduction_method)[data_df$sample_name]},
      error = function(e) {NULL})
  } else{
    data_df$cell_color <- colData(cds)[data_df$sample_name,color_cells_by]
  }
  
  ## Marker genes
  markers_exprs <- NULL
  if (!is.null(genes)) {
    if ((is.null(dim(genes)) == FALSE) && dim(genes) >= 2){
      markers <- unlist(genes[,1], use.names=FALSE)
    } else {
      markers <- genes
    }
    markers_rowData <-
      as.data.frame(subset(rowData(cds), gene_short_name %in% markers |
                             row.names(rowData(cds)) %in% markers))
    if (nrow(markers_rowData) >= 1) {
      cds_exprs <- SingleCellExperiment::counts(cds)[row.names(markers_rowData), ,drop=FALSE]
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))
      
      if ((is.null(dim(genes)) == FALSE) && dim(genes) >= 2){
        genes <- as.data.frame(genes)
        row.names(genes) <- genes[,1]
        genes <- genes[row.names(cds_exprs),]
        agg_mat <-
          as.matrix(my.aggregate.Matrix(cds_exprs,
                                        as.factor(genes[,2]),
                                        fun="sum"))
        agg_mat <- t(scale(t(log10(agg_mat + 1))))
        agg_mat[agg_mat < -2] <- -2
        agg_mat[agg_mat > 2] <- 2
        markers_exprs <- agg_mat
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        
        markers_exprs$feature_label <- markers_exprs$feature_id
        #markers_linear <- TRUE
      } else {
        cds_exprs@x <- round(10000*cds_exprs@x)/10000
        markers_exprs <- matrix(cds_exprs, nrow=nrow(markers_rowData))
        colnames(markers_exprs) <- colnames(SingleCellExperiment::counts(cds))
        row.names(markers_exprs) <- row.names(markers_rowData)
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        markers_exprs <- merge(markers_exprs, markers_rowData,
                               by.x = "feature_id", by.y="row.names")
        markers_exprs$feature_label <-
          as.character(markers_exprs$gene_short_name)
        markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <-
          markers_exprs$feature_id
        markers_exprs$feature_label <- factor(markers_exprs$feature_label,
                                              levels = markers)
      }
    }
  }
  
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name",
                     by.y="cell_id")
    data_df$expression <- with(data_df, ifelse(value >= min_expr, value, NA))
    sub1 <- data_df[!is.na(data_df$expression),]
    sub2 <- data_df[is.na(data_df$expression),]
    if(norm_method == "size_only"){
      p <- plotly::plot_ly(sub1) %>%
        plotly::add_trace(x = ~data_dim_1, y = ~data_dim_2, z = ~data_dim_3,
                          type = 'scatter3d', size=I(cell_size), alpha = I(alpha),
                          mode="markers", marker=list(
                            colorbar = list(title = "Expression", len=0.5),
                            color=~expression,
                            colors=color_scale,
                            line=list(width = 1,
                                      color = ~expression,
                                      colorscale=color_scale),
                            colorscale=color_scale)) %>%
        plotly::add_markers(x = sub2$data_dim_1, y = sub2$data_dim_2,
                            z = sub2$data_dim_3, color = I("lightgrey"),
                            size=I(cell_size),
                            marker=list(opacity = .4), showlegend=FALSE)
    } else {
      sub1$log10_expression <- log10(sub1$expression + min_expr)
      p <- plotly::plot_ly(sub1) %>%
        plotly::add_trace(x = ~data_dim_1, y = ~data_dim_2, z = ~data_dim_3,
                          type = 'scatter3d', size=I(cell_size), alpha = I(alpha),
                          mode="markers", marker=list(
                            colorbar = list(title = "Log10\nExpression", len=0.5),
                            color=~log10_expression,
                            colors=color_scale,
                            line=list(width = 1,
                                      color = ~log10_expression,
                                      colorscale=color_scale),
                            colorscale=color_scale)) %>%
        plotly::add_markers(x = sub2$data_dim_1, y = sub2$data_dim_2,
                            z = sub2$data_dim_3, color = I("lightgrey"),
                            size=I(cell_size),
                            marker=list(opacity = .4), showlegend=FALSE)
    }
  } else {
    if(color_cells_by %in% c("cluster", "partition")){
      if (is.null(data_df$cell_color)){
        p <- plotly::plot_ly(data_df, x = ~data_dim_1, y = ~data_dim_2,
                             z = ~data_dim_3, type = 'scatter3d',
                             size=I(cell_size), color=I("gray"),
                             mode="markers", alpha = I(alpha))
        message(paste("cluster_cells() has not been called yet, can't color",
                      "cells by cluster or partition"))
      } else{
        if(is.null(color_palette)) {
          N <- length(unique(data_df$cell_color))
          if(N > 8){
            color_palette <- as.character(paletteer::paletteer_d("rcartocolor::Pastel", N)) ###############################
          } else {
            color_palette <- as.character(paletteer::paletteer_d("rcartocolor::Pastel", N)) ###############################
          }
        }
        p <- plotly::plot_ly(data_df, x = ~data_dim_1, y = ~data_dim_2,
                             z = ~data_dim_3, type = 'scatter3d',
                             size=I(cell_size), color=~cell_color,
                             colors = color_palette,
                             mode="markers", alpha = I(alpha))
      }
    } else if(class(data_df$cell_color) == "numeric") {
      
      p <- plotly::plot_ly(data_df) %>%
        plotly::add_trace(x = ~data_dim_1, y = ~data_dim_2, z = ~data_dim_3,
                          type = 'scatter3d', size=I(cell_size), alpha = I(alpha),
                          mode="markers", marker=list(
                            colorbar = list(title = color_cells_by, len=0.5),
                            color=~cell_color,
                            colors=color_scale,
                            line=list(width = 1,
                                      color = ~cell_color,
                                      colorscale=color_scale),
                            colorscale=color_scale))
    } else {
      if(is.null(color_palette)) {
        N <- length(unique(data_df$cell_color))
        if(N > 8){
          color_palette <- as.character(paletteer::paletteer_d("rcartocolor::Pastel", N)) ###############################
        } else {
          color_palette <- as.character(paletteer::paletteer_d("rcartocolor::Pastel", N)) ###############################
        }
        names(color_palette) = levels(data_df$cell_color)
      }
      p <- plotly::plot_ly(data_df, x = ~data_dim_1, y = ~data_dim_2,
                           z = ~data_dim_3, type = 'scatter3d',
                           size=I(cell_size), color=~cell_color,
                           colors = color_palette,
                           mode="markers", alpha = I(alpha))
    }
  }
  p <- p %>%
    plotly::layout(scene = list(xaxis=list(title=paste("Component", x)),
                                yaxis=list(title=paste("Component", y)),
                                zaxis=list(title=paste("Component", z))))
  ## Graph info
  if (show_trajectory_graph) {
    
    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
      as.data.frame() %>%
      dplyr::select(prin_graph_dim_1 = x, prin_graph_dim_2 = y,
                    prin_graph_dim_3 = z) %>%
      dplyr::mutate(sample_name = rownames(.),
                    sample_state = rownames(.))
    
    dp_mst <- cds@principal_graph[[reduction_method]]
    
    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      dplyr::select(source = "from", target = "to") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select(source="sample_name",
                                       source_prin_graph_dim_1="prin_graph_dim_1",
                                       source_prin_graph_dim_2="prin_graph_dim_2",
                                       source_prin_graph_dim_3="prin_graph_dim_3"),
                       by = "source") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select(target="sample_name",
                                       target_prin_graph_dim_1="prin_graph_dim_1",
                                       target_prin_graph_dim_2="prin_graph_dim_2",
                                       target_prin_graph_dim_3="prin_graph_dim_3"),
                       by = "target")
    
    for (i in 1:nrow(edge_df)) {
      p <- p %>%
        plotly::add_trace(
          x = as.vector(t(edge_df[i, c("source_prin_graph_dim_1",
                                       "target_prin_graph_dim_1")])),
          y = as.vector(t(edge_df[i, c("source_prin_graph_dim_2",
                                       "target_prin_graph_dim_2")])),
          z = as.vector(t(edge_df[i, c("source_prin_graph_dim_3",
                                       "target_prin_graph_dim_3")])),
          color = trajectory_graph_color,
          line = list(color = I(trajectory_graph_color),
                      width = trajectory_graph_segment_size), mode = 'lines',
          type = 'scatter3d', showlegend = FALSE)
    }
  }
  p
}

"MetBrewer::Signac"