\name{defineTasks}
\alias{defineTasks}
\title{
  Defines a list of tasks to be undertaken by slave CPUs
}
\description{
  This function defines a list of PCIT tasks to be undertaken by slave CPUs in
  a parallel environment. By default, all possible 3-way combinations of genes
  are split such that each slave CPU has a single task to perform and these
  tasks each take approximately the same amount of time to complete. This
  limits the amount of bandwidth used for passing messages between the master
  and slave CPUs and ensures all slaves complete their task at about the same
  time.
}
\usage{
  defineTasks(nGenes, nSlaves, tasksPerSlave = 1, plot = FALSE)
}
\arguments{
  \item{nGenes}{ - The total number of genes. }
  \item{nSlaves}{ - The number of slave CPUs. }
  \item{tasksPerSlave}{ - The number of tasks to make for each slave CPU. }
  \item{plot}{ - A boolean as to whether to generate a plot showing the total
                 "work" to be carried out in each task. }
}
\value{
  A list of tasks to be performed by the slave CPUs. Tasks are defined by a
  range of x values which represent the row/col position of a gene on the
  adjacency matrix. A range of x values of 1:5 for a task means that the slave
  CPU will perform PCIT on only the three-way combinations of genes that have
  genes 1 to 5 in the 1st position of all possible combinations. For small
  values of x, the amount of "work" is more than that required for larger
  values of x. This is because there are more 3-way combinations of genes with
  gene 1 in the first position, than there are with gene 1000 in the 1st
  position.
}
\author{
  Nathan S. Watson-Haigh
}
