\name{maxMatrixSize}
\alias{maxMatrixSize}
\title{
	Calculate the maximum correlation matrix size which PCIT can handle
}
\description{
	This function attempts to determine the maximum sized correlation matrix
	which can be handled by PCIT given an amount of computer memory for the
	serial implementation.
}
\usage{
	maxMatrixSize(ram, units=c("MB", "bytes", "KB", "GB", "TB"), nCopies=3)
}
\arguments{
	\item{ram}{ - The amount of RAM memory available }
	\item{units}{ - The units in which RAM was specified }
	\item{nCopies}{ - The maximum number of copies of the correlation matrix which PCIT holds at any one time }
}
\author{
	Nathan S. Watson-Haigh
}
\examples{
	maxMatrixSize(1, "GB")
	maxMatrixSize(512, "MB")
}
\seealso{
	\code{\link{pcitMemoryRequirement}}
}
