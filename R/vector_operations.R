#' Calculate the norm of a row-vector.
#' @exports
row_norm <- function(x) {
    if (is.matrix(x)) {
        norm <- sqrt(rowSums(x^2))
    } else if (is.null(dim(x))) {
        norm <- sqrt(sum(x^2))
    } else {
        stop("Uknown object.")
    }

    return(norm)
}

#' Get distance of points to a line.
#' @exports
dist_to_line <- function(X, lines) {
    x_norm <- row_norm(X)

    if (is.matrix(lines)) {
        lines <- t(lines)
    }
    proj <- (X %*% lines)
    dist <- x_norm^2 - proj^2
    dist[dist < 0] <- 0
    dist <- sqrt(dist)

    return(dist)
}

#' Get euclidean distance between two vectors.
#' @exports
euclidean_dist <- function(a, b) {
    dist <- sqrt(sum((a - b)^2))
    return(dist)
}

#' Get cosine between row-vectors A and B.
#' @exports
get_cosine <- function(A, B) {
    norma <- row_norm(A)
    normb <- row_norm(B)

    if (is.null(dim(A)) && is.null(dim(B))) {
        costheta <- (A %*% B) / (norma * normb)
        costheta <- min(max(costheta, -1), 1)
    } else {
        if (!is.null(dim(B))) B <- t(B)

        costheta <- (A %*% B) / (outer(norma, normb, FUN = "*"))
        costheta <- pmin(pmax(costheta, -1), 1)
    }

    return(costheta)
}

#' Get angle between row-vectors A and B.
#' @exports
get_angle <- function(A, B) {
    angle <- acos(get_cosine(A, B))
    return(angle)
}


# vectorized version of the TS-SS
#' @exports
ts_ss <- function(A, B) {
    require("pdist")

    norma <- row_norm(A)
    normb <- row_norm(B)

    # Angle
    theta <- (A %*% t(B)) / (outer(norma, normb, FUN = "*"))
    theta <- pmin(pmax(theta, -1), 1)
    theta <- acos(theta) + deg2rad(10)

    # triangle similarity
    ts <- 0.5 * outer(norma, normb, FUN = "*") * sin(theta)

    # euclidean distance
    ed <- as.matrix(pdist::pdist(A, B))
    # Magnitude difference
    md <- abs(outer(norma, normb, FUN = "-"))

    # sector similarity
    ss <- (ed + md)^2 * (theta / 2)

    # TS-SS
    tsss <- ts * ss

    return(tsss)
}
