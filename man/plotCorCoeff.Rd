\name{plotCorCoeff}
\alias{plotCorCoeff}
\title{
  Superimpose a histogram of significant correlation coefficients over all
  correlation coefficents
}
\description{
  Given a complete correlation matrix and a list of non-significant correlation
  indices, plot the distribution of significant correlations on top of the
  distribution of all correlations.
}
\usage{
  plotCorCoeff(m, nonSignif, leg.txt = c("All", "PCIT Significant"), \dots)
}
\arguments{
  \item{m}{ - A correlation matrix. }
  \item{nonSignif}{ - Indices of all non-significant correlations. }
  \item{leg.txt}{ - A vector of labels to be used in the figure legend. }
  \item{\dots}{ - Additional parameters to be passed to rect(). }
}
\author{
  Nathan S. Watson-Haigh
}
