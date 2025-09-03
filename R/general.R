#' @export
is.empty <- function(x) {
  return(isTRUE(length(x) == 0 & !is.null(x)))
}

#' @export
nan_idxs <- function(x, out = "indx") {
  nans <- which(is.nan(x))
  cidx <- floor(nans / nrow(x)) + 1
  ridx <- nans %% nrow(x)

  if (out == "indx") {
    idxs <- list("ridx" = ridx, "cidx" = cidx)
  } else if (out == "list") {
    idxs <- vector(mode = "list", length = length(nans))
    for (i in seq_len(length(nans))) {
      idxs[[i]] <- c(ridx[i], cidx[i])
    }
  }
  return(idxs)
}

#' @export
rand_idx <- function(X, k) {
  idxs <- sample(seq_len(nrow(X)), size = k)
  return(idxs)
}

#' Convert a factor to a numeric.
#' @details
#' Assumes that the factors are numbers such as "1".
#' @param f A vector of factors.
#' @returns
#' The factors converted to numbers.
#' @export
f2n <- function(f) {
  f <- as.numeric(as.character(f))
  stopifnot(is.numeric(f))

  return(f)
}
