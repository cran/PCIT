\name{idx}
\alias{idx}
\title{
	Get indicies for significant edges
}
\description{
	Get the indicies for the significant edges in a network.
}
\usage{
	idx(result)
}
\arguments{
	\item{result}{ - A result object returned from pcit() }
}
\value{
	Linear indices are returned for those correlations found to be significant.
}
\examples{
	data(PCIT)
	m <- m[1:200,1:200]        # just use a small subset of the data
	result <- pcit(m)
	
	idx <- idx(result)
}
\author{
	Nathan S. Watson-Haigh
}
