\name{pcitMemoryRequirement}
\alias{pcitMemoryRequirement}
\title{
	Calculate the memory requirement for running PCIT
}
\description{
	This function attempts to determine the amount of computer memory that would
	be required to run PCIT of a given correlation matrix.
}
\usage{
	pcitMemoryRequirement(m, units=c("MB", "bytes", "KB", "GB", "TB"),
		nCopies=3)
}
\arguments{
	\item{m}{ - A correlation matrix on which PCIT may be run, OR the number of rows/columns from such a matirx }
	\item{units}{ - The units of RAM memory to be use for the returned value }
	\item{nCopies}{ - The maximum number of copies of the correlation matrix which PCIT holds at any one time }
}
\author{
	Nathan S. Watson-Haigh
}
\examples{
	m <- matrix(1, 20, 20)
	diag(m) <- 0
	m
	pcitMemoryRequirement(m, "KB")
	pcitMemoryRequirement(10000, "GB")
}
\seealso{
	\code{\link{maxMatrixSize}}
}
