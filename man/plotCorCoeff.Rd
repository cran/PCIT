\name{plotCorCoeff}
\alias{plotCorCoeff}
\title{
	Plot superimposed histograms of correlation coefficients
}
\description{
	Given a complete correlation matrix and a list of linear indices, superimpose
	the distributions of correlation coefficients, defined by their indices, on top of the
	distribution of all correlations.
	
	Only data from the upper triangle is used to plot the distributions since a correlation matrix should be symmetrical. This
	means that any specified indices which fall within the diagonal or the lower
	triangle are effectively ignored.
}
\usage{
	plotCorCoeff(m, idx, col, breaks="Scott", \dots)
}
\arguments{
	\item{m}{ - A correlation matrix representing the raw data. }
	\item{idx}{ - A list of indices for a subset of correlations to be superimposed on the plot. }
	\item{col}{ - A vector of colours of equal length to idx. These are used to colour each of the distributions specified by the indices in idx. }
	\item{breaks}{ - Defaults to "Scott", see the \code{\link{hist}} documentation for info and other options. }
  \item{\dots}{ - Additional parameters to be passed to rect(). }
}
\examples{
	data(PCIT)
	m <- m[1:200,1:200]        # just use a small subset of the data
	result <- pcit(m)
	
	op <- par(mfrow=c(2,1))
	plotCorCoeff(m, list("PCIT Significant" = idx(result)), col=c("black"))
	plotCorCoeff(m, list("PCIT Significant" = idx(result),
		"abs(adj) >= 0.7" = which(abs(m) >= 0.7)), col=c("black", "red"))
	par(op)
}
\author{
	Nathan S. Watson-Haigh
}
\seealso{
	\code{\link{hist}}
}

