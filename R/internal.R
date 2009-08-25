.setup_mpi <- function(nslaves=NULL) {
  # Default to running serially.  Only run in parallel if Rmpi
  # is installed AND there's multiple slaves AND pcit is called
  # from the master.
  run_parallel <- FALSE
  
  # Load the MPI Environment if not already there.
  if (!is.loaded('mpi_initialize')) {
    require('Rmpi', quietly=TRUE)
  }
  
  # Now, if Rmpi is loaded, make sure we have a bunch of slaves
  if (is.loaded('mpi_initialize')) {
    if (mpi.comm.size() > 2) {
      # It's already loaded and set up
      return(TRUE)
    }
    
    # Spawn as many slaves as possible
    if (mpi.comm.size() < 2) {
      if(is.null(nslaves)) {
        answer <- try(mpi.spawn.Rslaves(), silent=TRUE)
      } else {
        answer <- try(mpi.spawn.Rslaves(nslaves=nslaves), silent=TRUE)
      }
      if (class(answer) == "try-error") {
        cat(geterrmessage(), "\nWARN: loaded Rmpi library but couldn't spawn any slaves. Falling back to serial version of functions.\n")
        return(FALSE)
      }
    }
    
    # If there's multiple slaves, AND this process is the master,
    # run in parallel.
    if (mpi.comm.size() > 2) {
      if (mpi.comm.rank() == 0) {
        run_parallel <- TRUE
      }
    } else {
      cat("WARN: I could only find", mpi.comm.size()-1, "slaves to work with. Parallel programming usually requires >= 2 slaves, but I'll continue anyway!\n")
      run_parallel <- TRUE
    }
  }
  
  return(run_parallel)
}

.pcit <- function(m, x=1:(nrow(m)-2), tol.type=c("mean", "min", "max", "median"), verbose=getOption("verbose")) {
  # m is a matrix of direct correlations between genes using:
  # m <- similarity(d, type="none")
  
  if (! isSymmetric(m) ) {
    stop("The pcit() function requires a square, symmetrical matrix as input")
  }
  
  # use "mean" as the default for calculating the tolerances
  tol.type <- match.arg(tol.type)
  switch(tol.type,
    "mean" = { tol.type <- 1 },
    "min" = { tol.type <- 2 },
    "max" = { tol.type <- 3 },
    "median" = {
      tol.type <- 4
      stop("Tolerance calculation using median is not implemented in the PCIT Fortran code.")
    }
  )
  # set NA's to zero
  index <- is.na(m)
  if (sum(index)>0) {
    m[index] <- 0
    warning(paste(sum(index), " values were found to be N/A and set to zero.", sep=""))
  }
  
  m_partials <- m
  
  result <- .Fortran( "pcit", correlations=as.single(m), partial_correlations=as.single(m_partials), nGenes=as.integer(nrow(m)), xVals=as.integer(x), nXVals=as.integer(length(x)), tolType=as.integer(tol.type) )
  
  result$partial_correlations <- matrix(result$partial_correlations,nrow(m),ncol(m))
  result$nonsignif <- which(result$partial_correlations == 0, arr.ind=TRUE)
  
  # remove items from the results list which we don't want to return to the calling function
  result[c("correlations", "partial_correlations", "nGenes", "xVals", "nXVals", "tolType")] <- NULL
  
  # return only the idecies for those that PCIT found to be significant
  return(result)
}

.hist <- function(m, nonSignif, col=c("gray75", "black"), leg.txt=NULL, breaks="Scott", ...) {
  # plot just the data from upper tiangle of the 1st data matrix
  d1 <- m[upper.tri(m)]
  all.h <- hist(d1, plot=FALSE, breaks=breaks,)#, main="Density Plot of Raw Correlation Coefficients", xlab="Correlation Coefficient")
  
  # create a copy of the matrix, and zero out those positions not significant according to PCIT
  m2 <- m
  m2[nonSignif] <- 0
  # get just non-zero data from the upper triangle
  d2 <- m2[upper.tri(m2)]
  d2 <- d2[which(d2 != 0)]
  yh <- hist(d2, plot=FALSE, breaks=all.h$breaks)
  
  plot.new()
  plot.window(xlim=range(all.h$breaks, yh$breaks),
    ylim=range(0, all.h$counts, yh$counts))
  rect(all.h$breaks[-length(all.h$breaks)], 0,
    all.h$breaks[-1], all.h$counts, col=col[1], ...)
  rect(yh$breaks[-length(yh$breaks)], 0,
    yh$breaks[-1], yh$counts, col=col[2], ...)
  axis(1)
  axis(2)
  title(main="Density Distribution of Correlation Coefficients", xlab="Correlation Coefficient", ylab="Frequency")
  
  if(! is.null(leg.txt)) {
    legend(-1,max(all.h$counts), fill=col, legend=leg.txt)
  }
}

.clusterCoef <- function(x) {
  coef <- (length(which(x>0))-1)/(length(x)-1)
  return(coef)
}

.sub2ind <- function(x, y, nrow, ncol=NULL) {
  ## Returns a linear index for the (x,y) coordinates passed in.
  if (is.matrix(x) || is.data.frame(x)) {
    stopifnot(ncol(x) == 2)
    if (!missing(y)) {
      if (missing(nrow)) {
        nrow <- y
      } else {
        ncol <- nrow
        nrow <- y
      }
    }
    y <- x[,2]
    x <- x[,1]
  }
  
  if (is.matrix(nrow)) {
    d <- dim(nrow)
    nrow <- d[1]
    ncol <- d[2]
  } else if (is.null(ncol)) {
    stop("Dimensions of matrix under-specified")
  }
  
  # Sanity check to ensure we got each var doing what it should be doing
  if (length(x) != length(y) || length(nrow) != 1 || length(ncol) != 1) {
    stop("I'm confused")
  }
  
  ((x - 1) + ((y - 1) * nrow)) + 1
  
}

.Last <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0){
      #print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    #print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}
