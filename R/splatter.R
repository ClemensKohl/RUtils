
#' @exports
get_gene_clusters <- function(splatter_sim) {

  rd <- as.data.frame(rowData(splatter_sim))

  de_fac <- dplyr::select(rd, starts_with("DEFacGroup"))
  de_fac <- as.matrix(de_fac)

  rMs <- rowMaxs(de_fac)
  idxs <- purrr::map(seq_len(nrow(de_fac)), function(x) which(de_fac[x,] == rMs[x]))

  de_genes <- de_fac > 1 | de_fac < 1 # doesnt consider that gene can be DE in 2 groups
  de_bool <- matrix(FALSE,
                    nrow = nrow(de_fac),
                    ncol = ncol(de_fac),
                    dimnames = list(rownames(de_fac),
                                    colnames(de_fac)))
  for (x in seq_len(nrow(de_fac))) {

    if(any(de_genes[x,])) {
    vals <- de_fac[x,][de_genes[x,]]
    max_val <- which(abs(vals) == max(abs(vals)))

    de_bool[x, which(de_fac[x,] == vals[max_val])] <- TRUE
    }else {
      next
    }
  }

  grp <- as.numeric(de_bool %*% seq_len(ncol(de_bool)))
  grp_fact <- as.factor(grp)
  names(grp_fact) <- rownames(de_bool)

  return(grp_fact)
}
