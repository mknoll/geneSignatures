% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fract.R
\name{plotFract}
\alias{plotFract}
\title{Plot calculated fractions of singature
overlaps}
\usage{
plotFract(obj, type = "cont")
}
\arguments{
\item{obj}{singGene object}

\item{type}{plot continous fractions (cont)
or classified sample matrix (cat)}
}
\value{
list with elements GRP (classified) and
N (continuous fraction)
}
\description{
Plonts continuous or categorical
classifications of samples based on signatures
and returns the respective matrices
}
