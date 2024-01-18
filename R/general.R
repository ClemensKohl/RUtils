is.empty <- function(x) {
    return(isTRUE(length(x) == 0 & !is.null(x)))
}

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
