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
