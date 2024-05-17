#' Download and preprocess Zeisel brain data.
#' @export
get_zeisel_brain <- function(subset_by_MHVG = FALSE) {
    require(SingleCellExperiment)
    require(scran)
    require(scater)
    require(scuttle)
    require(scRNAseq)

    set.seed(2358)

    sce <- scRNAseq::ZeiselBrainData()

    # Preprocessing
    mt_genes <- grepl("^mt-", rownames(sce), ignore.case = TRUE)

    qc_df <- scuttle::perCellQCMetrics(sce, subsets = list(Mito = mt_genes))

    reasons <- scuttle::perCellQCFilters(
        qc_df,
        sum.field = "sum",
        detected.field = "detected",
        sub.fields = c("subsets_Mito_percent")
    )

    colData(sce) <- cbind(colData(sce), qc_df)
    sce$discard <- reasons$discard

    sce <- sce[, !reasons$discard]

    cnts <- as.matrix(counts(sce))
    genes_detect <- rowSums(cnts > 0) > (ncol(cnts) * 0.01)
    sce <- sce[genes_detect, ]

    clust <- scran::quickCluster(sce)
    sce <- scran::computeSumFactors(sce, cluster = clust, min.mean = 0.1)
    sce <- scuttle::logNormCounts(sce)

    dec <- scran::modelGeneVar(sce)
    top_genes <- scran::getTopHVGs(dec, prop = 0.8)
    sce <- scran::fixedPCA(sce, subset.row = top_genes)
    sce <- scater::runUMAP(sce, dimred = "PCA")

    if (isTRUE(subset_by_MHVG)) {
        sce <- sce[top_genes, ]
    }

    return(sce)
}

#' Preprocess any SingleCellExperiment
#' @export
preprocess_data_set <- function(
    sce,
    org,
    mt_filter = TRUE,
    pct = 0.01,
    modPoisson = FALSE,
    rmbatch = FALSE,
    batchcol = NULL,
    truth = NULL,
    return_plots = FALSE) {

    plots <- list()

    if (is(sce, "Seurat")) {
        sce <- SingleCellExperiment(
            list(counts = GetAssayData(sce,
                slot = "count",
                assay = "RNA"
            )),
            colData = sce@meta.data
        )
    }


    # Remove genes not expressed in any cells
    sce <- sce[Matrix::rowSums(counts(sce)) > 0, ]

    template <- rownames(sce)[1]

    if (!isEmpty(grep("EN", template))) {
        names <- stringr::word(rownames(sce), 1, 1, "_")
        rownames(sce) <- names

        if (!isEmpty(grep("ENSG", template))) {
            org.db <- org.Hs.eg.db
        }

        if (!isEmpty(grep("ENSMU", template))) {
            org.db <- org.Mm.eg.db
            org <- "mm"
        }

        if (!isEmpty(grep("[.]", template))) {
            rownames(sce) <- stringr::word(rownames(sce), 1, 1, "[.]")
        }

        sce <- getSCEWithSymbols(sce, keytype = "ENSEMBL", org.db = org.db)
    } else {
        rowData(sce)$SYMBOL <- rownames(sce)
    }

    # Filter based on QC metrics

    if (org == "hs") {
        mt_genes <- grepl("^MT-", rowData(sce)$SYMBOL)
    } else if (org == "mm") {
        mt_genes <- grepl("^mt-", rowData(sce)$SYMBOL)
    }

    if (isTRUE(mt_filter)) {
        qc_df <- perCellQCMetrics(sce, subsets = list(Mito = mt_genes))

        reasons <- perCellQCFilters(qc_df,
            sum.field = "sum",
            detected.field = "detected",
            sub.fields = c("subsets_Mito_percent")
        )
    } else {
        qc_df <- perCellQCMetrics(sce, subsets = list(Mito = mt_genes))

        reasons <- perCellQCFilters(qc_df,
            sum.field = "sum",
            detected.field = "detected",
            sub.fields = NULL
        )
    }

    # QC plots

    colData(sce) <- cbind(colData(sce), qc_df)
    sce$discard <- reasons$discard

    if (!is.null(truth)) {
        plots[["p_counts"]] <- plotColData(sce, x = truth, y = "sum", colour_by = "discard") +
            scale_y_log10() +
            ggtitle("Total count") +
            theme(axis.text.x = element_text(angle = 90))

        plots[["p_detected"]] <- plotColData(sce, x = truth, y = "detected", colour_by = "discard") +
            scale_y_log10() +
            ggtitle("Detected features") +
            theme(axis.text.x = element_text(angle = 90))

        plots[["p_mito"]] <- plotColData(sce, x = truth, y = "subsets_Mito_percent", colour_by = "discard") +
            ggtitle("Mito percent") +
            theme(axis.text.x = element_text(angle = 90))
    }

    plots[["p_sumVSdetected"]] <- plotColData(sce, x = "sum", y = "detected", colour_by = "discard") +
        ggtitle("sum vs detected")

    plots[["p_hist_detected"]] <- ggplot(as.data.frame(colData(sce)), aes(x = detected)) +
        geom_histogram(fill = "#008080") +
        ggtitle("detected") +
        theme_bw()


    plots[["p_hist_sum"]] <- ggplot(as.data.frame(colData(sce)), aes(x = sum)) +
        geom_histogram(fill = "#008080") +
        ggtitle("sum") +
        theme_bw()




    # remove bad cells
    sce <- sce[, !reasons$discard]

    # filter out genes expressed in less than 1% cells
    bm <- counts(sce) > 0

    ncell <- round(ncol(sce) * pct)
    label <- paste0("pct", pct)

    idx <- which((rowSums(bm) > ncell))
    sce <- sce[idx, ]

    # Normalization
    clust.sce <- quickCluster(sce)
    sce <- computeSumFactors(sce, cluster = clust.sce, min.mean = 0.1)
    sce <- logNormCounts(sce)

    if (isTRUE(rmbatch)) {
        if (is.null(batchcol)) {
            stop("value of parameter batchlab is missing.")
        }
        cnt_bc <- sva::ComBat(counts(sce), batch = colData(sce)[, batchcol])
        if (min(cnt_bc) < 0) {
            cnt_bc <- cnt_bc - min(cnt_bc)
        }
        cnt_bc <- cnt_bc[rowSums(cnt_bc) > 0, colSums(cnt_bc) > 0]
        counts(sce) <- cnt_bc
        sce <- logNormCounts(sce)
    }


    # dim reduc.
    if (isTRUE(modPoisson)) {
        sce.dec <- modelGeneVarByPoisson(sce)
        sce.top <- getTopHVGs(sce.dec, prop = 0.2)
    } else {
        sce.dec <- modelGeneVar(sce)
        sce.top <- getTopHVGs(sce.dec, prop = 0.2)
    }


    sce <- denoisePCA(sce, technical = sce.dec, subset.row = sce.top)
    sce <- runUMAP(sce, dimred = "PCA")

    p_umap <- plotUMAP(sce, colour_by = truth)

    # doublet detection
    sce <- scDblFinder(sce, clusters = clust.sce)
    # table(sce$scDblFinder.class)

    plots[["p_umap_dblscore"]] <- plotUMAP(sce, colour_by = "scDblFinder.score")

    plots[["p_umap_dblclass"]] <- plotUMAP(sce, colour_by = "scDblFinder.class")

    if (isTRUE(return_plots)) {
        return(list("sce" = sce, "plots" = plots))
    }
    return(sce)
}


# change gene names to gene symbol
getSCEWithSymbols <- function(sce, keytype = "ENTREZID", org.db = org.Mm.eg.db) {
    if (keytype == "SYMBOL") {
        rowData(sce)$SYMBOL <- rownames(sce)
        return(sce)
    }
    require("AnnotationDbi")
    require("org.Hs.eg.db")
    require("org.Mm.eg.db")

    geneSymbols <- mapIds(org.db, keys = rownames(sce), column = "SYMBOL", keytype = keytype, multiVals = "first")

    rowData(sce)$SYMBOL <- geneSymbols
    return(sce)
}
