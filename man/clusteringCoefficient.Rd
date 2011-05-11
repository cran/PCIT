\name{clusteringCoefficient}
\alias{clusteringCoefficient}
\title{
	Calculate the clustering coefficient
}
\description{
	Calculate the clustering coefficient for an adjacency matrix. By default, the local clustering coefficient
	is calculated.
}
\usage{
	clusteringCoefficient(adj, FUN='localClusteringCoefficient', \dots)
}
\arguments{
	\item{adj}{ - An adjacency matrix. Calculating the clustering coefficient
	only makes sense if some connections are zero  i.e. no connection.}
	\item{FUN}{ - The function for calculating the clustering coefficient.}
	\item{\dots}{ - Arguments to pass to FUN}
}
\value{
	The clustering coefficient(s) for the adjacency matrix.
}
\examples{
	data(PCIT)
	m <- m[1:200,1:200]        # just use a small subset of the data
	result <- pcit(m)
	m[idx(result)] <- 0
	
	clusteringCoefficient(m)
}
\seealso{
	\code{\link{localClusteringCoefficient}}
	\code{\link{clusteringCoefficientPercent}}
}
\author{
	Nathan S. Watson-Haigh
}
