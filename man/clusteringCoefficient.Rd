\name{clusteringCoefficient}
\alias{clusteringCoefficient}
\title{
  Calculate the clustering coefficient for each node/gene of an adjacency matrix
}
\description{
  The clustering coefficient is defined as the proportion of existing
  connections from the total possible (Watts and Strogatz, 1998).
}
\usage{
  clusteringCoefficient(adj)
}
\arguments{
  \item{adj}{ - An adjacency matrix. Calculating the clustering coefficient
  only makes sense if some connections are zero  i.e. no connection.}
}
\value{
  A vector of clustering coefficients for each node/gene of the adjacency matrix.
}
\references{
  D.J. Watts and S.H. Strogatz. (1998) Collective dynamics of 'small-world' networks. Nature.  393(6684). 440-442.
}
\author{
  Nathan S. Watson-Haigh
}
