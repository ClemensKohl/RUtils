#' Removes everything within a sphere
#'  of the defined quantile of the vector norm.
#'
#' @param x matrix of row vectors
#' @param qcutoff quantile.
#'
#' @returns
#' Matrix of the vectors longer than the defined cutoff.
#' @export
sphere_cutoff <- function(x, qcutoff = 0.8) {
  xn <- row_norm(x)
  q <- quantile(xn, qcutoff)
  x <- x[xn > q, ]
  return(x)
}

#' Returns the indices of all points with a norm outside of a sphere
#'  with radius of the defined quantile of the vector norm.
#'
#' @param x matrix of row vectors
#' @param qcutoff quantile.
#'
#' @returns
#' Indices of points lying outside of sphere
#' @export
ca_sphere_idx <- function(x, qcutoff = 0.8) {
  xn <- row_norm(x)
  q <- quantile(xn, qcutoff)
  idx <- which(xn > q)

  return(idx)
}


#' prin_coords: principal coordinates (row-vectors)
#' D: singular vectors.
#' @export
prin_to_std <- function(prin_coords, D) {
  if (is.null(dim(prin_coords))) {
    stopifnot(length(prin_coords) == length(D))
    std_coords <- prin_coords / D
  } else {
    stopifnot(ncol(prin_coords) == length(D))
    std_coords <- sweep(
      prin_coords,
      2,
      D,
      "/"
    )
  }

  return(std_coords)
}
