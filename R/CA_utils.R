#' Removes everything within a sphere
#'  of the defined quantile of the vector norm.
#'
#' @param x matrix of row vectors
#' @param qcutoff quantile.
#'
#' @returns
#' Matrix of the vectors longer than the defined cutoff.
sphere_cutoff <- function(x, qcutoff = 0.8) {
    xn <- row_norm(x)
    q <- quantile(xn, qcutoff)
    x <- x[xn > q, ]
    return(x)
}


#' prin_coords: principal coordinates (row-vectors)
#' D: singular vectors.
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
