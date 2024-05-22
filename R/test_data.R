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
    # require(SingleCellExperiment)
    # require(scater)
    # require(scran)
    # require(scuttle)
    # require("AnnotationDbi")
    # require("org.Hs.eg.db")
    # require("org.Mm.eg.db")
    # require(scDblFinder)
    # require(Seurat)

    plots <- list()

    if (is(sce, "Seurat")) {
        sce <- SingleCellExperiment::SingleCellExperiment(
            list(counts = Seurat::GetAssayData(sce,
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
            org.db <- org.Hs.eg.db::org.Hs.eg.db
        }

        if (!isEmpty(grep("ENSMU", template))) {
            org.db <- org.Mm.eg.db::org.Mm.eg.db
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
        mt_genes <- base::grepl("^MT-", rowData(sce)$SYMBOL)
    } else if (org == "mm") {
        mt_genes <- base::grepl("^mt-", rowData(sce)$SYMBOL)
    }

    if (isTRUE(mt_filter)) {
        qc_df <- scuttle::perCellQCMetrics(sce, subsets = list(Mito = mt_genes))

        reasons <- scuttle::perCellQCFilters(qc_df,
            sum.field = "sum",
            detected.field = "detected",
            sub.fields = c("subsets_Mito_percent")
        )
    } else {
        qc_df <- scuttle::perCellQCMetrics(sce, subsets = list(Mito = mt_genes))

        reasons <- scuttle::perCellQCFilters(qc_df,
            sum.field = "sum",
            detected.field = "detected",
            sub.fields = NULL
        )
    }

    # QC plots

    SummarizedExperiment::colData(sce) <- cbind(SummarizedExperiment::colData(sce), qc_df)
    sce$discard <- reasons$discard

    if (!is.null(truth)) {
        plots[["p_counts"]] <- scater::plotColData(sce, x = truth, y = "sum", colour_by = "discard") +
            ggplot2::scale_y_log10() +
            ggplot2::ggtitle("Total count") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

        plots[["p_detected"]] <- scater::plotColData(sce, x = truth, y = "detected", colour_by = "discard") +
            ggplot2::scale_y_log10() +
            ggplot2::ggtitle("Detected features") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

        plots[["p_mito"]] <- scater::plotColData(sce, x = truth, y = "subsets_Mito_percent", colour_by = "discard") +
            ggplot2::ggtitle("Mito percent") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
    }

    plots[["p_sumVSdetected"]] <- scater::plotColData(sce, x = "sum", y = "detected", colour_by = "discard") +
        ggplot2::ggtitle("sum vs detected")

    plots[["p_hist_detected"]] <- ggplot2::ggplot(
        as.data.frame(SummarizedExperiment::colData(sce)), ggplot2::aes(x = detected)
    ) +
        ggplot2::geom_histogram(fill = "#008080") +
        ggplot2::ggtitle("detected") +
        ggplot2::theme_bw()


    plots[["p_hist_sum"]] <- ggplot2::ggplot(as.data.frame(SummarizedExperiment::colData(sce)), ggplot2::aes(x = sum)) +
        ggplot2::geom_histogram(fill = "#008080") +
        ggplot2::ggtitle("sum") +
        ggplot2::theme_bw()

    # remove bad cells
    sce <- sce[, !reasons$discard]

    # filter out genes expressed in less than 1% cells
    bm <- SingleCellExperiment::counts(sce) > 0

    ncell <- round(ncol(sce) * pct)
    label <- paste0("pct", pct)

    idx <- which((Matrix::rowSums(bm) > ncell))
    sce <- sce[idx, ]

    # Normalization
    clust.sce <- scran::quickCluster(sce)
    sce <- scran::computeSumFactors(sce, cluster = clust.sce, min.mean = 0.1)
    sce <- scuttle::logNormCounts(sce)

    if (isTRUE(rmbatch)) {
        if (is.null(batchcol)) {
            stop("value of parameter batchlab is missing.")
        }
        cnt_bc <- sva::ComBat(counts(sce), batch = SummarizedExperiment::colData(sce)[, batchcol])
        if (min(cnt_bc) < 0) {
            cnt_bc <- cnt_bc - min(cnt_bc)
        }
        cnt_bc <- cnt_bc[Matrix::rowSums(cnt_bc) > 0, Matrix::colSums(cnt_bc) > 0]
        SingleCellExperiment::counts(sce) <- cnt_bc
        sce <- scuttle::logNormCounts(sce)
    }


    # dim reduc.
    if (isTRUE(modPoisson)) {
        sce.dec <- scran::modelGeneVarByPoisson(sce)
        sce.top <- scran::getTopHVGs(sce.dec, prop = 0.2)
    } else {
        sce.dec <- scran::modelGeneVar(sce)
        sce.top <- scran::getTopHVGs(sce.dec, prop = 0.2, var.threshold = NULL)
    }


    sce <- scran::denoisePCA(sce, technical = sce.dec, subset.row = sce.top)
    sce <- scater::runUMAP(sce, dimred = "PCA")

    plots[["p_umap"]] <- scater::plotUMAP(sce, colour_by = truth)

    # doublet detection
    sce <- scDblFinder::scDblFinder(sce, clusters = clust.sce)
    # table(sce$scDblFinder.class)

    plots[["p_umap_dblscore"]] <- scater::plotUMAP(sce, colour_by = "scDblFinder.score")

    plots[["p_umap_dblclass"]] <- scater::plotUMAP(sce, colour_by = "scDblFinder.class")

    if (isTRUE(return_plots)) {
        return(list("sce" = sce, "plots" = plots))
    }
    return(sce)
}


#' change gene names to gene symbol
getSCEWithSymbols <- function(sce, keytype = "ENTREZID", org.db = org.Mm.eg.db) {
    if (keytype == "SYMBOL") {
        SummarizedExperiment::rowData(sce)$SYMBOL <- rownames(sce)
        return(sce)
    }
    require("AnnotationDbi")
    require("org.Hs.eg.db")
    require("org.Mm.eg.db")

    geneSymbols <- mapIds(org.db, keys = rownames(sce), column = "SYMBOL", keytype = keytype, multiVals = "first")

    SummarizedExperiment::rowData(sce)$SYMBOL <- geneSymbols
    return(sce)
}


#' Seurat vignette and conversion to SCE
#' @export
sce_pbmc3k <- function() {
    sce <- TENxPBMCData::TENxPBMCData(dataset = "pbmc3k")
    rownames(sce) <- make.unique(SummarizedExperiment::rowData(sce)$Symbol_TENx)
    colnames(sce) <- SummarizedExperiment::colData(sce)$Barcode

    pbmc <- SeuratObject::CreateSeuratObject(
        counts = as.matrix(SingleCellExperiment::counts(sce)),
        assay = "RNA",
        project = "pbmc3k",
        min.cells = 3,
        min.features = 200,
        meta.data = as.data.frame(SingleCellExperiment::colData(sce))
    )

    pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")

    # Filter data
    pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    no_zeros_rows <- Matrix::rowSums(pbmc, slot = "counts") > 0
    pbmc <- pbmc[no_zeros_rows, ]

    # Normalization
    pbmc <- Seurat::NormalizeData(pbmc,
        normalization.method = "LogNormalize",
        scale.factor = 10000,
        verbose = FALSE
    )

    pbmc <- Seurat::FindVariableFeatures(pbmc,
        # selection.method = "vst",
        nfeatures = 2000,
        verbose = FALSE
    )

    # Scaling
    all.genes <- rownames(pbmc)
    pbmc <- Seurat::ScaleData(pbmc,
        features = all.genes,
        verbose = FALSE
    )

    # Run PCA
    pbmc <- Seurat::RunPCA(pbmc,
        features = Seurat::VariableFeatures(object = pbmc),
        verbose = FALSE
    )

    # Cell clustering
    pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
    pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5, verbose = FALSE)

    pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10, verbose = FALSE)

    new.cluster.ids <- c(
        "0 - Naive CD4 T",
        "1 - CD14+ Mono",
        "2 - Memory CD4 T",
        "3 - B",
        "4 - CD8 T",
        "5 - FCGR3A+ Mono",
        "6 - NK",
        "7 - DC",
        "8 - Platelet"
    )

    names(new.cluster.ids) <- levels(pbmc)
    pbmc <- SeuratObject::RenameIdents(pbmc, new.cluster.ids)
    pbmc$cell_type <- SeuratObject::Idents(pbmc)

    return(Seurat::as.SingleCellExperiment(pbmc))
}
