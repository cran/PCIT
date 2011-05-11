\name{idxInvert}
\alias{idxInvert}
\title{
	Invert linear indices from a matrix
}
\description{
	Given a matrix from which which a subset of linear indices were obtained,
	invert those indices.
}
\usage{
	idxInvert(m, idx)
}
\arguments{
	\item{m}{ - A matrix from which idx is a subset of linear indices, OR the number of rows/columns from such a matirx }
	\item{idx}{ - A vector containing a subset of linear indices from the matrix m }
}
\author{
	Nathan S. Watson-Haigh
}
\examples{
	m <- matrix(1, 5, 5)
	diag(m) <- 0
	m
	idx <- which(m==0)
	idx
	idxInvert(m, idx)
	idxInvert(5, idx)
}
