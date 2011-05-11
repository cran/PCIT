\name{clusteringCoefficientPercent}
\alias{clusteringCoefficientPercent}
\title{
	Calculate the clustering coefficient as a percentage
}
\description{
	Given an adjacency matrix, calculate the clustering coefficient as a percentage of non-zero adjacencies.
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
\examples{
	data(PCIT)
	m <- m[1:200,1:200]        # just use a small subset of the data
	result <- pcit(m)
	m[idx(result)] <- 0
	
	clusteringCoefficientPercent(m)
}
\seealso{
	\code{\link{clusteringCoefficient}}
}
\author{
	Nathan S. Watson-Haigh
}
