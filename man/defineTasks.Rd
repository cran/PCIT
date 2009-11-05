\name{defineTasks}
\alias{defineTasks}
\title{
  Define a list of tasks for slave CPUs
}
\description{
  This function defines a list of PCIT tasks to be undertaken by slave CPUs in
  a parallel environment. The set of calculations performed for each trio of
	genes are entirely independent. Workload is distributed evenly among the
	slaves in order to maximise speed by assigning sets of gene trios to a task.
	The task is then carried out by a slave CPU. This assignment of gene trios
	into sets is based on which gene is in the first position in a given trio of
	genes. For example, if we have 5 genes (A, B, C, D and E) we have a set
	containing 6 gene trios with A in the first position (ABC, ABD, ABE, ACD, ACE
	and ADE), a set containing 3 gene trios with B in the first position (BCD,
	BCE and BDE) and a set containing 1 gene trio with C in the first position
	(CDE) giving a total of 10 gene trios (given by C(n,3), where n is the total
	number of genes) assigned to 1 of 3 sets. The computational workload for each
	trio is the same, thus the workload for a set of trios is directly
	proportional to the number of trios it contains, C(n-m,k-1), Where n is the
	total number of genes and the gene constrained in the first position of a
	gene trio is in the mth row/column in the correlation matrix and takes values
	in the interval [1,n-2].
	
	By default, the number of tasks to be created is equal to the number for
	slave CPUs. Therefore, multiple gene trio sets are assigned to each task but
	in such a way as to balance the work load. Since gene trio sets created with
	small values of m have many more gene trio members than those created with
	larger values of m. A plot can be produced which shows the amount of work
	taken for a gene trio set constrained with the mth gene in the first
	position. The plot also shows how the cumlative work for all n-2 gene trio
	sets is divided into nSlave*tasksPerSlave equally sized tasks, so they can be
	completed in approximately the same amount of time.

}
\usage{
  defineTasks(n, nSlaves, tasksPerSlave = 1, plot = FALSE)
}
\arguments{
  \item{n}{ - The total number of genes. }
  \item{nSlaves}{ - The number of slave CPUs. }
  \item{tasksPerSlave}{ - The number of tasks to create for each slave CPU. }
  \item{plot}{ - A boolean as to whether to generate a plot showing the total
                 "work" to be carried out by each task. }
}
\author{
  Nathan S. Watson-Haigh
}
\examples{
	defineTasks(n=100, nSlaves=5, plot=TRUE)
}
