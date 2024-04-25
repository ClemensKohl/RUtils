#' Total Least Squares Regression for a set of points.
#' @param points A matrix of points. Rows are observations.
#' @exports
tls_regression <- function(points) {
    require(irlba)

    if (nrow(points) == 1) {
        reg_line <- points / row_norm(points)
    } else {
        suppressWarnings({
            reg_line <- irlba::irlba(points, nv = 1, right_only = TRUE)$v
        })
    }

    return(reg_line)
}

#' Get slope of a set of lines.
#' @exports
slope <- function(lines, dims = 1:2) {
    if (is.null(nrow(lines))) {
        d <- lines[dims[2]] / lines[dims[1]]
    } else {
        d <- lines[, dims[2]] / lines[, dims[1]]
    }
    return(d)
}
