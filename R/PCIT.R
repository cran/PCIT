pcit <- function(m, force.serial=FALSE, force.parallel=FALSE, nslaves=NULL, pass.file=TRUE, verbose=getOption("verbose"), tol.type=c("mean", "min", "max", "median") ) {
  # Create data structure to store the results
  results <- list()
  results$nonsignif <- NULL
  
  if(force.serial) {
    if(force.parallel) {
      cat("WARN: You specified both force.serial=TRUE and force.parallel=TRUE. I'll assume you want me to choose which implementation to use!\n")
      
    } else {
      run_parallel <- FALSE
    }
  } else {
    # If serial isn't forced, we'll try and run pcit in parallel, run_parallel, will tell us if there is a suitable parallel environment or not
    run_parallel <- .setup_mpi(nslaves=nslaves)
  }
  
  if(!run_parallel) {
    # no parallel environment, or serial was forced
    if(force.parallel) {
      stop("Unable to force the parallel implementation of pcit()")
    }
    if(verbose) { cat("Running serial implementation of pcit().\n") }
    results <- .pcit(m, verbose=verbose, tol.type=tol.type)
    results$nonsignif <- .sub2ind(results$nonsignif, nrow=nrow(m), ncol=ncol(m))
    
  } else {
    # lets run the parallel implementation
    if(verbose) { cat("Running parallel implementation of pcit().\n") }
    
    # We're gonna pass the data to the slaves in a file if requested
    if(pass.file) {
      RData <- tempfile()
      save(m, file=RData)
    } else {
      cat("WARN: setting pass.file=FALSE may provide a small amount of speedup by passing data around in memory.\n  However, there is a limit to the size of the data set that can be passed in memory.\n  If you receive an 'Error: serialization is too large to store in a raw vector', you have reached this limit and should use pass.file=TRUE\n")
    }
    
    nSlaves <- mpi.comm.size()-1
    
    # Function the slaves will call to perform a validation on the
    # fold equal to their slave number.
    # Assumes the following object(s) will be passed to it by the master: m
    slavefunction <- function() {
      # load the pcit library in each slave
      require('PCIT', quietly=TRUE)
      
      # Note the use of the tag for sent messages: 
      #     1=ready_for_task, 2=done_task, 3=exiting 
      # Note the use of the tag for received messages: 
      #     1=task, 2=done_tasks 
      junk <- 0
      done <- 0
      while (done != 1) {
        # Signal being ready to receive a new task 
        mpi.send.Robj(junk,0,1)
        
        # Receive a object list from master
        objList <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
        task_info <- mpi.get.sourcetag()
        tag <- task_info[2]
        
        if (tag == 1) {
          # do a task, taking data/objects from the objList
          # Find thos entries that are zero'd out in the m_sub_pcit matrix and convert them back to indices for the full matrix
          # Lets not check to see see which zero's were changed only by this slave as we'd need more time and memory to do so - i.e. a nGenes x nGenes Boolean matrix for each check
          if (pass.file) {
            # load the whole data matrix from file and run PCIT() on it
            load(RData)
            results <- .pcit(m, x=objList$x, verbose=verbose, tol.type=tol.type)
            results$nonsignifOffset <- 1
            
          } else {
            # get the portion of the data matrix and do PCIT() on it
            results <- .pcit(objList$m, x=objList$x, verbose=verbose, tol.type=tol.type)
            # convert array indicies returned by pcit into indices equivilent to the full matrix
            results$nonsignif <- results$nonsignif + objList$nonsignifOffset - 1
            results$nonsignifOffset <- objList$nonsignifOffset
            
          }
          
          results$nonsignif <- .sub2ind(results$nonsignif, nrow=nrow(m), ncol=ncol(m))
          
          # Send a results message back to the master
          mpi.send.Robj(results,0,2)
          
        } else if (tag == 2) {
          # maybe send date() info back to the master so it can compile debugging info on the distribution of times each slave took?
          done <- 1
          if (verbose) {
            cat("Slave", mpi.comm.rank(), "- FINISHED:", date(), "\n")
          }
        }
      # We'll just ignore any unknown messages
      }
      
      mpi.send.Robj(junk,0,3)
    }
    
    # We're in the parent.  
    # send data object(s) to all the slaves
    mpi.bcast.Robj2slave(tol.type)
    
    if (pass.file) {
      mpi.bcast.Robj2slave(RData)
    }
    # Send function(s) to all the slaves
    mpi.bcast.Robj2slave(slavefunction)
    # Call the function in all the slaves to get them ready to
    # undertake tasks
    mpi.bcast.cmd(slavefunction())
    
    # define a list of tasks
    tasks <- defineTasks(nGenes=nrow(m), nSlaves=nSlaves)
    
    junk <- 0
    exited <- 0
    child_errors <- FALSE
    while (exited < nSlaves) { 
      # Receive a message from a slave 
      message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
      message_info <- mpi.get.sourcetag() 
      slave_id <- message_info[1] 
      tag <- message_info[2] 
      
      if (tag == 1) {
        # slave is ready for a task. Give it the next task, or tell it tasks 
        # are done if there are none.
        # If errors have been found in children, don't give any new tasks
        if (child_errors) {
          cat("WARN: Telling slave ",slave_id," to exit\n");
          mpi.send.Robj(junk,slave_id,2)
        }
        if (length(tasks) > 0) {
          if (pass.file) {
            # send the job with data loaded from file
            mpi.send.Robj(tasks[[1]], slave_id, 1)
            
          } else {
            # send only a portion of the data matrix to be as memory efficient as possible
            xMin <- min(tasks[[1]]$x)
            # get the minimaml size matrix to pass to the slave
            # off-set the values of x to correspond to the correct values for m.sub
            x.new <- tasks[[1]]$x - xMin + 1
            
            # send this info as a list of objects to the slave to deal with
            objList <- list(m=m[xMin:nrow(m), xMin:nrow(m)], indOffset=xMin, x=x.new)
            
            # Send a task, and then remove it from the task list
            mpi.send.Robj(objList, slave_id, 1)
          }
          
          tasks[[1]] <- NULL 
          
        } else {
          mpi.send.Robj(junk, slave_id, 2)
        }
        
      } else if (tag == 2) { 
        # The message contains results
        # compile each subset into a master result
        
        # for linear ind which are positions of non-significnat connections - basically, a zero in one/both of results and message
        results$nonsignif <- unique(c(results$nonsignif, message$nonsignif))
        
        # for linear ind which are positions of significant connections - basically, a 1 in both results and message
        # BUT - only from the same overlapping sub matrix - need to append those indices from outside the overlapping region of the larger matrix
        # the following DOES NOT work!
        #results$nonsignif <- c(results$nonsignif, message$nonsignif)[duplicated(c(results$nonsignif, message$nonsignif))]
        
        #results$nonsignif <- append(results$nonsignif, message$nonsignif)
        
      } else if (tag == 3) {
        # A slave has closed down.
        exited <- exited + 1
      } else if (tag == 4) {
        exited <- exited + 1
        child_errors <- TRUE
        cat("ERROR: Slave ", slave_id, ": ", message,"\n")
      }
    }
    
    # clean up
    if (pass.file) {
      unlink(RData)
    }
    if (child_errors) {
      stop("Errors in R slaves.  See last.warnings or check slave logs")
    }
  }
  # Back to doing things in common with both serial and parallel implementations
  
  return(results)
}

# based on equal "work" per task
defineTasks <- function(nGenes, nSlaves, tasksPerSlave=1, plot=FALSE) {
  
  # calculate how much "work" each value of x will need
  work <- vector()
  for(x in 1:(nGenes-2)) {
    work[x] <- ((nGenes-1-x)*(nGenes-x))/2    # this is just the sum of 1:nGenes-2 but more efficient for large values of nGenes
  }
  
  slice <- function(v, n){
    subtot <- floor(sum(v)/n)
    cumtot <- cumsum(v)
    p <- rep(0,n-1)
    for(i in 1:(n-1)) {
      p[i] <- max(which(cumtot < (subtot*i)))
    }
    p <- append(p, 0, after=0)
    p <- append(p, length(work))
    
    return(p)
  }
  
  
  work_slices <- slice(work, nSlaves*tasksPerSlave)
  
  tasks <- vector('list')
  colours <- vector()
  for(i in 1:(length(work_slices)-1)) {
    start <- work_slices[i]+1
    end <- work_slices[i+1]
    
    tasks[[length(tasks)+1]] <- list(x=(start:end))
    colours <- append(colours,rep(i,end-start+1))
  }
  
  if (plot) {
    par(mfrow=c(2:1))
    plot(work, col=colours, ylab="Work", xlab="x", main=paste("Work per x - Colours Showing Split of x into",nSlaves*tasksPerSlave,"tasks"))
    plot(cumsum(work)/max(cumsum(work))*100, col=colours, pch=20, ylab="Proportion of Total Work", xlab="x", main=paste("Cumlative Work - Colours Showing Split of x into",nSlaves*tasksPerSlave,"tasks"))
    abline(h=c(cumsum(work)[c(work_slices[2:(length(work_slices)-1)])]/max(cumsum(work))*100))
  }
  
  return(tasks)
}

clusteringCoefficient <- function(adj) {
  m_coef <- apply(adj, 1, FUN=.clusterCoef)
  return(m_coef)
}

clusteringCoefficientPercent <- function(adj) {
  d <- adj[upper.tri(adj)]
  
  d1 <- d[which(d != 0)]
  return(length(d1) / length(d) *100)
}

plotCorCoeff <- function(m, nonSignif, leg.txt=c("All", "PCIT Significant"), ...) {
  .hist(m, nonSignif, leg.txt=leg.txt, ...)
}

# below are my scribbles for creating an S4 complient class/methods etc
#setClass("pcit",
#  representation(
#    dend = "NULL"
#  ),
#  contains="matrix",
#  prototype=prototype(
#    new("matrix")
#  )
#)

#setMethod("show", "pcit",
#  function(object) {
#    cat("pcit object\n")
#    cat("  Type       :", class(object), "\n")
#    cat("  Size       :", paste(dim(object), collapse="x"), "\n")
#    cat("  Dendrogram :", if(is.null(object@dend)) { "FALSE" } else { "TRUE" }, "\n")
#    cat("\n")
#  }
#)

#setMethod("image", "pcit",
#  function(x, labRow=FALSE, labCol=FALSE, scale="none", ...) {
#    heatmap(
#      x@.Data,
#      Rowv = x@dend,
#      Colv = x@dend,
#      labRow = labRow,
#      labCol = labCol,
#      scale = scale,
#      ...
#    )
#  }
#)

#setValidity("pcit",
#  function(object) {
#    retval <- NULL
#    if(! isSymmetric(object)) {
#      retval <- c(retval, "Data matrix must be symmetrical")
#    }
#    if(is.null(retval)) return(TRUE)
#    else return(retval)
#  }
#)
