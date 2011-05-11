\name{data}
\alias{data}
\alias{m}
\non_function{}
\keyword{
	datasets
}
\title{
	Demo Data
}
\description{
	Example/mock correlation data set for use as input for running the demo.
}
\usage{
	data(PCIT)
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
