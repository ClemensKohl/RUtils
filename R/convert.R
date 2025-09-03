#' Convert degrees to radians.
#' @export
deg2rad <- function(x) {
  r <- (x * pi) / 180
  return(r)
}

#' Convert radians to degrees.
rad2deg <- function(x) {
  d <- (x * 180) / pi
  return(d)
}

#' Converts two vectors of gene and cell clusters to biclustlib results.
#'
#' @param cell_clusters Named vector of cell clusters.
#' @param gene_clusters Named vector of gene clusters.
#'
#' @return
#' An object of type "Biclust".
#'
#' @export
bic_to_biclust <- function(cell_clusters, gene_clusters, params = list()) {
  require("biclust")

  ctypes <- sort(unique(cell_clusters))
  gtypes <- sort(unique(gene_clusters))
  bitypes <- union(ctypes, gtypes)

  number <- length(bitypes)

  if (number == 0) {
    number_x_col <- matrix(0)
    row_x_number <- matrix(0)
  } else {
    number_x_col <- do.call(
      rbind,
      lapply(bitypes, function(x) {
        cell_clusters == x
      })
    )
    row_x_number <- do.call(
      cbind,
      lapply(bitypes, function(x) {
        gene_clusters == x
      })
    )
  }

  rownames(row_x_number) <- names(gene_clusters)
  colnames(row_x_number) <- paste0("BC_", bitypes)

  rownames(number_x_col) <- paste0("BC_", bitypes)
  colnames(number_x_col) <- names(cell_clusters)

  bic <- new(
    "Biclust",
    "Parameters" = params,
    "RowxNumber" = row_x_number,
    "NumberxCol" = number_x_col,
    "Number" = number,
    "info" = list("Generated from cell and gene clusters.")
  )

  return(bic)
}
