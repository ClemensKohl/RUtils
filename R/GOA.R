
#' Load the required gene set.
#' @description
#' Loads the speciefied gene set and subsets to the required organism.
#' @param set Name of the gene set. Currently only supports "CellMarker"
#' @param org Short name of the organism. "mm" for mouse, "hs" for human.
#' @returns
#' data frame with columns "cell_type" and "gene".
#' @exports
load_gene_set <- function(set = "CellMarker",
                          org) {

    stopifnot(org %in% c("mm", "hs"))

    if (set == "CellMarker") {

        gs <- CAbiNet::cellmarker_v2

        if (org == "mm") gs <- gs[gs$species == "Mouse", ]
        if (org == "hs") gs <- gs[gs$species == "Human", ]

        gs <- gs[, c("cell_name", "marker")]

    } else {
        stop("Other gene sets besides CellMarker are not implemented yet.")
    }

    colnames(gs) <- c("cell_type", "gene")

    return(gs)
}

#' Filter gene sets by size
#' @description
#' This function cuts removes gene sets below and above a certain size.
#' @param gene_sets A named list of gene sets (name is the gene set,
#' each element of the list contains a vector of genes)
#' @param min_size Min. number of genes in the gene set.
#' @param max_size Max number of genes in the gene set.
#' Set to Inf if you want to keep all genes.
#' @returns
#' Filtered list of gene sets.
#' @exports
filter_gene_sets <- function(gene_sets,
                             min_size = 10,
                             max_size = 500) {

    if (is.na(min_size) || is.null(min_size))
        min_size <- 1
    if (is.na(max_size) || is.null(max_size))
        max_size <- Inf #.Machine$integer.max

    ## index of gene_sets in used.
    ## logical
    gene_sets_length <- lengths(gene_sets)
    idx <- min_size <= gene_sets_length & gene_sets_length <= max_size

    return(gene_sets[idx])
}

#' Changes long format gene set data frame into a named list of gene sets.
#'
#' @description
#' Takes in a data frame with 2 columns (gene set name and gene name)
#' and accumulates all genes in the same gene set into a list element.
#'
#' @param gene_sets data frame with first column gene set names and the
#'  second column the gene symbol. Long format expected.
#'
#' @returns
#' A named list with gene sets as names and genes as vectors.
#' @exports
format_gene_sets <- function(gene_sets) {

    # Input gene_sets assumed to be a long format data frame
    # Bring data frame into shape that you can use for hypergeom test

    if (is(gene_sets, "tbl")) gene_sets <- as.data.frame(gene_sets)

    gene_list <- list()

    # unique gene set names
    gs <- sort(unique(drop(gene_sets[, 1])))

    for (i in seq_len(length(gs))) {

        sel <- which(gene_sets[, 1] == gs[i])
        genes <- drop(gene_sets[sel, 2])
        gene_list[[gs[i]]] <- genes
    }

    names(gene_list) <- gsub(" ", "_", names(gene_list))
    return(gene_list)
}

#' Adapted from `DOSE:::enricher_interal`.
#' Performs Gene Overrepresentation Analysis.
#'
#' @description
#' perform_goa takes a number of genes of interest (gois)
#' and a list of gene sets and performs gene overrepresentation
#' analysis.
#'
#' @param gois Genes of interest. Usually the co-clustered genes.
#' @param universe All genes in data set.
#' @param gene_sets Named list of gene sets and their genes.
#' @inheritParams filter_gene_sets
#'
#' @details
#' perform_goa performs a hypergeometric test on the gois and gene sets.
#' Gene sets are trimmed to genes that are present in universe
#'  (usually all the genes available in the data set) and trimmed to set
#'  min and max size.
#'
#' @returns
#' A data frame with the results of the go-analysis.
#'
#' @export
perform_goa <- function(gois,
                        gene_sets,
                        universe,
                        min_size,
                        max_size) {

    # subset gene sets to genes in universe
    gene_sets <- lapply(gene_sets, intersect, universe)

    # Subset gois to genes in gene sets
    all_gs <- unique(unlist(gene_sets))
    gois <- gois[gois %in% all_gs]

    # Filter out very small and very large gene sets
    # We do this after subsetting the gois.
    gene_sets <- filter_gene_sets(gene_sets = gene_sets,
                                  min_size = min_size, # 10
                                  max_size = max_size) # 500

    # number of clustered genes in each gene set
    gois_in_set <- sapply(gene_sets, intersect, gois)

    # Remove gene sets with 0 gois in them
    gois_in_set <- filter_gene_sets(gene_sets = gois_in_set,
                                    min_size = 1,
                                    max_size = Inf)

    # subset gene sets to thos with gois in them
    # make sure the two sets are the same order.
    gs_names <- sort(unique(names(gois_in_set)))
    gene_sets <- gene_sets[gs_names]
    gois_in_set <- gois_in_set[gs_names]

    # Build parameter data frame
    ngois <- length(gois)                   # clustered genes
    group1 <- lengths(gene_sets)            # the size of gene sets
    group2 <- length(all_gs)                # total genes in gene sets
    overlap <- lengths(gois_in_set)         # nr gois in gene_set

    phyper_df <- data.frame(
        gois_in_set = overlap - 1,          # white balls drawn / gois in gs
        genes_in_set = group1,              # total white balls / genes in gs
        genes_universe = group2 - group1,   # total black balls / all genes in gene set - gs
        ngois = ngois                       # balls drawn / number gois
    )
    rownames(phyper_df) <- names(gene_sets)

    # Hypergeometric test
    pvalues <- apply(phyper_df, 1, function(n) {
        stats::phyper(n[1], n[2], n[3], n[4], lower.tail = FALSE)
    })

    # adjust p-values
    p_adj <- stats::p.adjust(pvalues, method = "BH")


    ## gene ratio and background ratio
    gene_ratio <- apply(data.frame(a = overlap, b = ngois), 1, function(x) {
        paste(x[1], "/", x[2], sep = "", collapse = "")
    })

    bg_ratio <- apply(data.frame(a = group1, b = group2), 1, function(x) {
        paste(x[1], "/", x[2], sep = "", collapse = "")
    })

    # return results either as data frame or list
    enrich_res <- data.frame(gene_set = names(gene_sets),
                             pval = pvalues,
                             padj = p_adj,
                             GeneRatio = gene_ratio,
                             BgRatio = bg_ratio,
                             ngois_in_set = overlap,
                             ngenes_in_set = group1,
                             ngois = ngois,
                             ngenes_in_sets = group2)

    rownames(enrich_res) <- NULL

    ord <- order(enrich_res$padj)

    return(enrich_res[ord, ])
}

#' Perform gene set overrepresentation analysis
#' for each bicluster and annotate cells based on
#' the best match.
#'
#' @description
#' per_cluster_goa loads the required gene set, formats it and performs
#' gene overrepresentation analysis for each bicluster in the cabic
#' object.
#'
#' @param gene_clusters Named vector of gene clusters
#'  as obtained from `caclust`.
#' @inheritParams perform_goa
#' @inheritParams load_gene_set
#'
#' @return
#' A list contain the goa results for each cluster.
#'
#' @export
per_cluster_goa <- function(gene_clusters,
                            universe,
                            org,
                            set = "CellMarker",
                            min_size = 10,
                            max_size = 500) {

    # Load gene sets
    gs <- load_gene_set(set = set, org = org)
    gene_sets <- format_gene_sets(gs)

    gc_names <- sort(unique(gene_clusters))

    gc_list <- lapply(X = gc_names, FUN = function(x) {
        names(gene_clusters[which(gene_clusters == gc_names[x])])
    })

    names(gc_list) <- gc_names

    # Perform goa for each cluster
    goa_res <- list()

    for (c in seq_len(length(gc_list))){

        gene <- gc_list[[c]]
        clst_name <- names(gc_list)[c]

        goa <- perform_goa(gois = gene,
                           gene_sets = gene_sets,
                           universe = universe,
                           min_size = min_size,
                           max_size = max_size)

        goa_res[[clst_name]] <- goa
    }

    return(goa_res)
}

#' turns string of a ratio in a number.
#' @description
#' Adapted from DOSE::parse_ratio.
#'
#' @param ratio a string of the form "1/2"
#'
#' @returns
#' The numeric value the string represented.
#' @exports
parse_ratio <- function(ratio) {

    ratio <- sub("^\\s*", "", as.character(ratio))
    ratio <- sub("\\s*$", "", ratio)
    numerator <- as.numeric(sub("/\\d+$", "", ratio))
    denominator <- as.numeric(sub("^\\d+/", "", ratio))
    return(numerator / denominator)
}

#' Plot gene overrepresentation analysis results.
#'
#' @description
#' Plots the gene overrepresentation analysis results for each cluster.
#'
#' @param goa_res List of goa results for each cluster.
#' @param nres Number of top gene sets to plot.
#'
#' @returns
#' A ggplot2 object. Dotplot of the top `nres` gene sets per cluster.
#' @export
plot_goa_res <- function(goa_res,
                         nres = 3) {

    goa_res <- lapply(goa_res, function(x, nc = nres) {
        subs <- min(nrow(x), nres)
        x[seq_len(subs), ]
    })

    goa_res <- dplyr::bind_rows(goa_res, .id = "cluster")
    goa_res$GeneRatio <- parse_ratio(goa_res$GeneRatio)

    goa_res$gene_set <- factor(goa_res$gene_set,
                               levels = unique(goa_res$gene_set))

    p <- ggplot2::ggplot(goa_res, ggplot2::aes(
        x = cluster,
        y = gene_set,
        size = GeneRatio,
        color = padj
    )) +
    ggplot2::geom_point() +
    viridis::scale_color_viridis(direction = -1) +
    ggplot2::labs(
        title = paste0("Top ", nres, " gene sets per cluster"),
        x = "Cluster",
        y = "Gene Set"
    ) +
    ggplot2::theme_bw()

    return(p)
}
