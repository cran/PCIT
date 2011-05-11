\name{localClusteringCoefficient}
\alias{localClusteringCoefficient}
\title{
	Calculate the local clustering coefficient
}
\description{
	Calculate the local clustering coefficient for each node in an adjacency matrix.
	The clustering coefficient is defined as the proportion of existing
	connections from the total possible (Watts and Strogatz, 1998).
}
\usage{
	localClusteringCoefficient(adj)
}
\arguments{
	\item{adj}{ - An adjacency matrix. Calculating the clustering coefficient
	only makes sense if some connections are zero  i.e. no connection.}
}
\value{
	A vector of local clustering coefficients for each node/gene of the adjacency matrix.
}
\references{
	D.J. Watts and S.H. Strogatz. (1998) Collective dynamics of 'small-world' networks. Nature.  393(6684). 440-442.
}
\examples{
	data(PCIT)
	m <- m[1:200,1:200]        # just use a small subset of the data
	result <- pcit(m)
	m[idx(result)] <- 0
	
	localClusteringCoefficient(m)
}
\seealso{
	\code{\link{clusteringCoefficient}}
}
\author{
	Nathan S. Watson-Haigh
}
