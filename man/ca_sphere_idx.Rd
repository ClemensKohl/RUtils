% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CA_utils.R
\name{ca_sphere_idx}
\alias{ca_sphere_idx}
\title{Returns the indices of all points with a norm outside of a sphere
with radius of the defined quantile of the vector norm.}
\usage{
ca_sphere_idx(x, qcutoff = 0.8)
}
\arguments{
\item{x}{matrix of row vectors}

\item{qcutoff}{quantile.}
}
\value{
Indices of points lying outside of sphere
}
\description{
Returns the indices of all points with a norm outside of a sphere
with radius of the defined quantile of the vector norm.
}
