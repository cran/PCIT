\name{getEdgeList}
\alias{getEdgeList}
\title{
  Converts an adjacency matrix into edge list representation
}
\description{
  Given an adjacency matrix, converts it into edge list representation. This
  edge list can be written to a file for easy import into other software such
  as cytoscape. 
  
  Only the upper triangle is returned as it is assumed the matrix is symmetric.
  The edge list is returned as a data frame with 3 columns: 'From', 'To' and
  'Weight'.
}
\usage{
  getEdgeList(m, rm.zero=TRUE)
}
\arguments{
  \item{m}{ - An adjacency matrix. }
  \item{rm.zero}{ - A boolean to indicate whether zero weight edges should be excluded from the edge list. }
}
\examples{
	data(PCIT)
	m <- m[1:200,1:200]        # just use a small subset of the data
	
	el <- getEdgeList(m)
	
	# modify the edge list to include some useful attributes for cytoscape
	el$sign[el$Weight<0] <- '-'
	el$sign[el$Weight>0] <- '+'
	el$Weight <- abs(el$Weight)
	# write the edge list stuff to a file suitable for import into cytoscape
	write.table(el, file="el.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
}
\author{
  Nathan S. Watson-Haigh
}
\seealso{
  \code{\link{hist}}
}

