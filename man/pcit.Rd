\name{pcit}
\alias{pcit}
\title{
  Apply the PCIT algorithm to a correlation matrix
}
\description{
  Given a correlation matrix the PCIT algorithm is applied to identify
  non-significant correlations. If a parallel environment running Rmpi is
  detected, a parallel implementation will be run unless force.serial=TRUE
}
\usage{
  result <- pcit(m, force.serial=FALSE, force.parallel=FALSE, nslaves=NULL,
    pass.file=TRUE, verbose=getOption("verbose"), tol.type=c("mean", "min", "max", "median"))
}
\arguments{
  \item{m}{ - A correlation matrix. }
  \item{force.serial}{ - A boolean to indicate if the serial implementation of PCIT should be forced. }
  \item{force.parallel}{ - A boolean to indicate if the parallel implementation of PCIT should be forced. }
  \item{nslaves}{ - The number of slaves to spawn. By default, as many slaves as possible are spawned. UNTESTED OPTION. }
  \item{pass.file}{ - A boolean to indicate if the correlation matrix should be passed between the master and slave CPUs as a file. }
  \item{verbose}{ - A boolean to indicate if verbose output should be used. }
  \item{tol.type}{ - The type of tolerance measure to be used in PCIT. Current options are "mean", "min" and "max". }
}
\value{
  Linear indices are returned for those correlations found to be non-significant.
}
\references{
  Reverter, A. & Chan, E.K., 2008. Combining partial correlation and an information theory approach to the reversed-engineering of gene co-expression networks. Bioinformatics, btn482.
}
\author{
  Nathan S. Watson-Haigh
}
