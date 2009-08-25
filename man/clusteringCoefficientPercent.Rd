\name{clusteringCoefficientPercent}
\alias{clusteringCoefficientPercent}
\title{
  Calculate the clustering coefficient measured as the percentage of non-zero correlation coefficients
}
\description{
  Given a correlation matrix the PCIT algorithm is applied. Linear indices are
  returned for those correlations found to be non-significant.
}
\usage{
  clusteringCoefficientPercent(adj)
}
\arguments{
  \item{adj}{ - An adjacency matrix. Calculating the clustering coefficient
  percentage only makes sense if some connections are zero  i.e. no connection. }
}
\value{
  A numerical between 0 and 100.
}
\author{
  Nathan S. Watson-Haigh
}
