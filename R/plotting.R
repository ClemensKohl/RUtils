#' Plots a heatmap ordered by cell and gene clusters of the logcounts.
#' @export
plot_hm <- function(
  sce,
  cell_clusters,
  gene_clusters,
  show_column_title = character(0),
  show_row_title = character(0),
  row_title_rot = 90,
  show_row_names = FALSE
) {
    require(SingleCellExperiment)
    require(scran)
    require(ComplexHeatmap)
    require(ggthemes)
    require(grid)
    require(circlize)
    require(matrixStats)

    mat <- as.matrix(SingleCellExperiment::logcounts(sce[names(gene_clusters), names(cell_clusters)]))

    mat <- (mat - rowMeans(mat)) / matrixStats::rowSds(mat)
    mat[is.na(mat)] <- 0

    clust_order_cells <- c()
    un_cc <- sort(unique(cell_clusters))
    for (i in seq_along(un_cc)) {
        clust_order_cells <- c(clust_order_cells, which(cell_clusters == un_cc[i]))
    }

    clust_order_genes <- c()
    un_gc <- sort(unique(gene_clusters))
    for (i in seq_along(un_gc)) {
        clust_order_genes <- c(clust_order_genes, which(gene_clusters == un_gc[i]))
    }

    mat <- mat[clust_order_genes, clust_order_cells]

    ccs <- cell_clusters[colnames(mat)]
    gcs <- gene_clusters[rownames(mat)]

    col_fun <- circlize::colorRamp2(
        c(quantile(mat, 0.01), 0, quantile(mat, 0.99)),
        c("#2E5A87", "white", "#A90C38")
    )

    if (length(un_cc) <= 10) {
        pal <- ggthemes::tableau_color_pal("Tableau 10")
    } else {
        pal <- ggthemes::tableau_color_pal("Tableau 20")
    }

    annocol <- pal(length(un_cc))

    names(annocol) <- as.character(un_cc)

    column_ha <- ComplexHeatmap::HeatmapAnnotation(
        cluster = ccs,
        col = list(cluster = annocol),
        show_legend = c(FALSE),
        annotation_name_gp = grid::gpar(fontsize = 0)
    )

    row_ha <- ComplexHeatmap::rowAnnotation(
        cluster = gcs,
        col = list(cluster = annocol),
        show_legend = c(FALSE),
        annotation_name_gp = grid::gpar(fontsize = 0)
    )

    breaks <- c(min(-1, quantile(mat, 0.01)), 0, 2, 4, max(6, quantile(mat, 0.99)))

    hm <- ComplexHeatmap::Heatmap(mat,
        name = "z-score",
        heatmap_legend_param = list(
            col_fun = col_fun,
            title = "z-score",
            at = breaks,
            grid_width = unit(7, "mm"),
            legend_height = unit(5, "cm"),
            title_gp = grid::gpar(fontsize = 23),
            labels_gp = grid::gpar(fontsize = 18)
        ),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_row_names = show_row_names,
        row_title_side = "left",
        col = col_fun,
        top_annotation = column_ha,
        left_annotation = row_ha,
        use_raster = FALSE,
        row_split = gcs,
        column_split = ccs,
        column_title = show_column_title,
        column_title_gp = gpar(fontsize = 14),
        row_title = show_row_title,
        row_title_gp = gpar(fontsize = 14),
        row_title_rot = row_title_rot
    )

    return(hm)
}



#' @export
make_bic_hm <- function(bic,
                        max_bics = bic@Number,
                        algorithm = "Replace with algorithm name") {
    nbics <- bic@Number
    max_bics <- min(max_bics, bic@Number)

    clust_order_cells <- c()
    for (b in seq_len(nbics)) {
        idx <- which(bic@NumberxCol[b, ] == TRUE)
        idx <- idx[!idx %in% clust_order_cells]
        clust_order_cells <- c(clust_order_cells, idx)
    }

    idx <- seq_len(ncol(bic@NumberxCol))
    idx <- idx[!idx %in% clust_order_cells]
    clust_order_cells <- c(clust_order_cells, idx)

    clust_order_genes <- c()
    for (b in seq_len(nbics)) {
        idx <- which(bic@RowxNumber[, b] == TRUE)
        idx <- idx[!idx %in% clust_order_genes]
        clust_order_genes <- c(clust_order_genes, idx)
    }
    idx <- seq_len(nrow(bic@RowxNumber))
    idx <- idx[!idx %in% clust_order_genes]
    clust_order_genes <- c(clust_order_genes, idx)

    pal <- tableau_color_pal("Tableau 10")
    annocol <- pal(nbics)

    col_df <- bic@NumberxCol
    col_df[bic@NumberxCol] <- "y"
    col_df[!bic@NumberxCol] <- "n"
    col_df <- as.data.frame(t(col_df))
    col_df <- col_df[clust_order_cells, , drop = FALSE]
    col_df <- col_df[, seq_len(max_bics), drop = FALSE]

    col_df <- rev(col_df)

    anno_list <- list()
    for (b in seq_len(nbics)) {
        anno_list[[paste0("BC", b)]] <- c("y" = annocol[b], "n" = "white")
    }
    anno_list <- anno_list[seq_len(max_bics)]
    anno_list_cols <- anno_list

    column_ha <- HeatmapAnnotation(
        df = col_df,
        col = anno_list_cols,
        show_legend = c(FALSE),
        na_col = "white"
    )

    row_df <- bic@RowxNumber
    print(dim(row_df))
    row_df[bic@RowxNumber] <- "y"
    row_df[!bic@RowxNumber] <- "n"
    row_df <- as.data.frame(row_df)
    print(dim(row_df))
    row_df <- row_df[clust_order_genes, , drop = FALSE]
    row_df <- row_df[, seq_len(max_bics), drop = FALSE]

    row_df <- rev(row_df)

    anno_list_rows <- anno_list

    row_ha <- rowAnnotation(
        df = row_df,
        col = anno_list_rows,
        show_legend = c(FALSE),
        annotation_name_gp = gpar(fontsize = 0),
        na_col = "white"
    )

    col_fun <- colorRamp2(seq(quantile(cnts, 0.01), quantile(cnts, 0.99), length = 10), viridis(10))

    hm <- Heatmap(
        matrix = cnts,
        col = col_fun,
        name = " ",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_title = algorithm,
        top_annotation = column_ha,
        left_annotation = row_ha,
        use_raster = TRUE
    )

    return(hm)
}


# Get relative coordinates of points in a ggplot2 plot.
# FIXME: Prevent opening plot when calling it! ggplotGrob
#
# Solution adapted from Anwer by Allan Cameron at:
# https://stackoverflow.com/a/60857307/1376616
# TODO: Add documentation!
# TODO: Better adapt to your plotting with a changing panel size.
#' @export
get_x_y_values <- function(gg_plot) {
    img_dim      <- grDevices::dev.size("cm") * 10
    # Attempt to silence the empty plot window.
    png("NUL")
    gt <- ggplot2::ggplotGrob(gg_plot)
    dev.off()
    to_mm        <- function(x) grid::convertUnit(x, "mm", valueOnly = TRUE)
    n_panel      <- which(gt$layout$name == "panel")
    panel_pos    <- gt$layout[n_panel, ]
    panel_kids   <- gtable::gtable_filter(gt, "panel")$grobs[[1]]$children
    point_grobs  <- panel_kids[[grep("point", names(panel_kids))]]
    from_top     <- sum(to_mm(gt$heights[seq(panel_pos$t - 1)]))
    from_left    <- sum(to_mm(gt$widths[seq(panel_pos$l - 1)]))
    from_right   <- sum(to_mm(gt$widths[-seq(panel_pos$l)]))
    from_bottom  <- sum(to_mm(gt$heights[-seq(panel_pos$t)]))
    panel_height <- img_dim[2] - from_top - from_bottom
    panel_width  <- img_dim[1] - from_left - from_right
    xvals        <- as.numeric(point_grobs$x)
    yvals        <- as.numeric(point_grobs$y)
    yvals        <- yvals * panel_height + from_bottom
    xvals        <- xvals * panel_width + from_left
    data.frame(x = xvals/img_dim[1], y = yvals/img_dim[2])
}

# Plot sankey plot of dataframe (columns are steps/nodes)
# TODO: Add documentation
#
#' @export
plot_flow <- function(df, rm_redund = TRUE) {

    if (isTRUE(rm_redund)) {
        for (c in seq_len(ncol(df))) {
            if (c == 1) next
            if (all(df[, c] == df[, c - 1])) {
                sub_cls <- df[, -c]
            }
        }
    }

    sank <- ggsankey::make_long(df,
                                colnames(df))
    p <- ggplot2::ggplot(sank,
                        ggplot2::aes(x = x,
                                     next_x = next_x,
                                     node = node,
                                     next_node = next_node,
                                     fill = factor(node),
                                     label = node)) +
                      ggsankey::geom_sankey(node_color = 1, flow_alpha = 0.7)  +
                      ggsankey::geom_sankey_label(size = 3.5, color = 1, fill = "white") +
                      ggplot2::scale_fill_viridis_d(option = "A", alpha = 0.95) +
                      ggsankey::theme_sankey(base_size = 12) +
                      ggplot2::theme(legend.position = "none")
    return(p)
}

