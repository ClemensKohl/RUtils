#TODO: Fix requirements!

# TODO: Add documentation
# TODO: Make more general!

#' @exports
de_comp <- function(biclust, ref_sce, truth_col, topdegs = 200){


      colLabels(ref_sce) <- factor(colData(ref_sce)[,truth_col])
      degs <- scran::scoreMarkers(ref_sce)

      gcs <- biclust@RowxNumber
      gc_clusters <- colnames(gcs)

      intermat <- matrix(NA,
                         nrow = length(degs),
                         ncol = ncol(gcs))

      rownames(intermat) <- names(degs)
      colnames(intermat) <- colnames(gcs)


      for (ct in seq_len(length(degs))){

        degenes <- degs[[ct]]
        degenes <- degenes[order(degenes$min.logFC.cohen, decreasing = TRUE),]
        degenes <- degenes[seq_len(topdegs),]
        degenes <- as.data.frame(degenes) %>%
                      tibble::rownames_to_column("gene")

        inter <- vector(mode="numeric", length = ncol(gcs))

        for (n in seq_len(ncol(gcs))){
            marker_genes <- rownames(gcs)[which(gcs[,n] == TRUE)]
            overlap_genes <- length(base::intersect(degenes$gene, marker_genes))
            inter[n] <- overlap_genes/topdegs
        }

        intermat[ct,] <- inter

      }

      return(intermat)

}

# TODO: Add documentation
# TODO: Make more general!

#' @exports
plot_de_intersection <- function(intermat){

  intermat <- as.data.frame(intermat)
  cnames <- colnames(intermat)
  p <- intermat %>%
          tibble::rownames_to_column("cell_type") %>%
          tidyr::pivot_longer(
              cnames,
              names_to = "cluster",
              values_to = "perc_intersection") %>%
          ggplot2::ggplot(aes(x=cluster, y = cell_type, fill = perc_intersection)) +
              ggplot2::geom_tile() +
              viridis::scale_fill_viridis(limits = c(0, 1))+
              ggplot2::theme_bw() +
              ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1))

  return(p)
}
