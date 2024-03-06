#' Plots a heatmap ordered by cell and gene clusters of the logcounts.
plot_hm <- function(sce, cell_clusters, gene_clusters) {
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
        show_row_names = FALSE,
        row_title_side = "left",
        col = col_fun,
        top_annotation = column_ha,
        left_annotation = row_ha,
        use_raster = FALSE,
        row_split = gcs,
        column_split = ccs,
        column_title_gp = gpar(fontsize = 14),
        row_title_gp = gpar(fontsize = 14)
    )

    return(hm)
}



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