\name{pcit}
\alias{pcit}
\title{
	Apply the PCIT algorithm
}
\description{
	Given a correlation matrix the PCIT algorithm (Reverter & Chan 2008) is applied to identify
	significant correlations. If a parallel environment running Rmpi is
	detected, a parallel implementation will be run unless force.serial=TRUE
}
\usage{
	pcit(m, force.serial=FALSE, force.parallel=FALSE, nslaves=NULL,
		verbose=getOption("verbose"),
		tol.type=c("mean", "min", "max", "median"),
		pass.type=c("file", "memory", "db"))
}
\arguments{
	\item{m}{ - A correlation matrix. }
	\item{force.serial}{ - A boolean to indicate if the serial implementation of PCIT should be forced. }
	\item{force.parallel}{ - A boolean to indicate if the parallel implementation of PCIT should be forced. }
	\item{nslaves}{ - The number of slaves to spawn. By default, as many slaves as possible are spawned. UNTESTED OPTION. }
	\item{verbose}{ - A boolean to indicate if verbose output should be used. }
	\item{tol.type}{ - The type of tolerance measure to be used in PCIT. Current options are "mean", "min" and "max". }
	\item{pass.type}{ - The type of approach used to pass the correlation matrix from the master CPU to the slave CPUs. Current options are "file", "memory" and "db". }
}
\value{
	Linear indices are returned for those correlations found to be significant.
}
\references{
	Reverter, A. & Chan, E.K., 2008. Combining partial correlation and an information theory approach to the reversed-engineering of gene co-expression networks. Bioinformatics, btn482.
}
\examples{
	data(PCIT)
	m <- m[1:200,1:200]        # just use a small subset of the data
	
	result <- pcit(m)
}
\author{
	Nathan S. Watson-Haigh
}
