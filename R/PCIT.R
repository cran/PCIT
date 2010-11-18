pcit <- function(m, force.serial=FALSE, force.parallel=FALSE, nslaves=NULL, verbose=getOption("verbose"), tol.type=c("mean", "min", "max", "median"), pass.type=c("file", "memory", "db") ) {
	pass.type <- match.arg(pass.type)
	tol.type <- match.arg(tol.type)
	
	# Create data structure to store the results
	results <- list()
	results$idx <- NULL
	
	if(force.serial && force.parallel) {
		cat("  WARN: You specified both force.serial=TRUE and force.parallel=TRUE. I'll assume you want me to choose which implementation to use!\n")
	} else if(force.parallel) {
		# attempt to setup parallel environment
		run_parallel <- .setup_mpi()
		if(!run_parallel) { stop("  Unable to force the parallel implementation of pcit()") }
	} else if(force.serial) {
		run_parallel <- FALSE
	} else {
		# neither force.serial or force.parallel set to true
		# attempt to autodetect parallel environment
		run_parallel <- .setup_mpi()
	}
	
	
	if (run_parallel) {
		if(verbose) { cat("  Running parallel implementation of pcit().\n") }
	} else {
		if(verbose) { cat("  Running serial implementation of pcit().\n") }
		results <- .pcit(m, verbose=verbose, tol.type=tol.type)
		results$idx <- .sub2ind(results$idx, nrow=nrow(m), ncol=ncol(m))
		return(results)
	}
	
	# we only get here if we have a detected a parallel environment
	# We're gonna pass the data to the slaves using one of the possible options
	switch(pass.type,
			file = {
				RData <- tempfile()
				save(m, file=RData)
			},
			memory = {
				cat("WARN: setting pass.file=FALSE may provide a small amount of speedup by passing data around in memory.\n  However, there is a limit to the size of the data set that can be passed in memory.\n  If you receive an 'Error: serialization is too large to store in a raw vector', you have reached this limit and should use pass.type=\"file\" (the default)\n")
			},
			db = {
				stop("Using a DB to pass data is not yet implemented.\n")
			}
	)
	
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
				switch(pass.type,
						file = {
							# load the whole data matrix from file and run PCIT() on it
							load(RData)
							results <- .pcit(m, x=objList$x, verbose=verbose, tol.type=tol.type)
							results$nonsignifOffset <- 1
						},
						memory = {
							results <- .pcit(m, x=objList$x, verbose=verbose, tol.type=tol.type)
							results$nonsignifOffset <- 1
						},
						db = {
							stop("Using a DB to pass data is not yet implemented.\n")
						}
				)
				
				results$idx <- .sub2ind(results$idx, nrow=nrow(m), ncol=ncol(m))
				
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
	
	switch(pass.type,
			file = {
				mpi.bcast.Robj2slave(RData)
			},
			memory = {
				# This means each slave has it's own copy of the master correlation matrix
				mpi.bcast.Robj2slave(m)
			},
			db = {
				stop("Using a DB to pass data is not yet implemented.\n")
			}
	)
	
	# Send function(s) to all the slaves
	mpi.bcast.Robj2slave(slavefunction)
	# Call the function in all the slaves to get them ready to
	# undertake tasks
	mpi.bcast.cmd(slavefunction())
	
	# define a list of tasks
	tasks <- defineTasks(n=nrow(m), nSlaves=nSlaves)
	
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
				switch(pass.type,
						file = {
							# send the job with data loaded from file
							mpi.send.Robj(tasks[[1]], slave_id, 1)
						},
						memory = {
							mpi.send.Robj(tasks[[1]], slave_id, 1)
						},
						db = {
							stop("Using a DB to pass data is not yet implemented.\n")
						}
				)
				
				tasks[[1]] <- NULL 
				
			} else {
				mpi.send.Robj(junk, slave_id, 2)
			}
			
		} else if (tag == 2) { 
			# The message contains results
			# compile each subset into a master result
			
			# for linear ind which are positions of non-significnat connections - basically, a zero in one/both of results and message
			if(verbose) {
				cat("Master: calculating unique indicies ... ")
			}
			
			if (is.null(results$idx)) {
				# if the master doesn't have any indices yet, just assign the result from the slave
				# The slave will only have zero'd out connections within the sub matrix it was dealing with
				results$idx <- message$idx
			} else {
				# combining data from slave to that which the master already knows about
				results$idx <- intersect(results$idx, message$idx)
			}
			
			if(verbose) {
				cat("DONE\n")
			}
			
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
	if(verbose) {
		cat("Master: completed pcit().\n")
	}
	# clean up
	switch(pass.type,
			file = {
				unlink(RData)
			},
			memory = {
				
			},
			db = {
				stop("Using a DB to pass data is not yet implemented.\n")
			}
	)
	if (child_errors) {
		stop("Errors in R slaves.  See last.warnings or check slave logs")
	}
	
	.freeSlaves()
	
	return(results)
}

# based on equal "work" per task
defineTasks <- function(n, nSlaves, tasksPerSlave=1, plot=FALSE) {
	
	# calculate how much "work" each value of x will need
	# this is equivilant to nCm(n-m+1,3) - nCm(n-m,3) where n = number of genes and m is in the set [1:n-2]
	m <- 1:(n-2)
	work <- ((n-m)*(n-m-1))/2
	
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
	
	
	slice_boundaries <- slice(work, nSlaves*tasksPerSlave)
	
	tasks <- vector('list')
	colours <- vector()
	for(i in 1:(length(slice_boundaries)-1)) {
		start <- slice_boundaries[i]+1
		end <- slice_boundaries[i+1]
		
		tasks[[length(tasks)+1]] <- list(x=(start:end))
		colours <- append(colours,rep(i,end-start+1))
	}
	
	if (plot) {
		par(mfrow=c(2:1))
		# Colours indicate the sets of gene trios assigned to the same task. By default, there is 1 task per slave CPU
		plot(work, col=colours, pch=20, ylab="Work", xlab="Gene Trio Set (m)", main="Work per Set of Gene Trios")
		
		plot(cumsum(work)/max(cumsum(work))*100, col=colours, pch=20, ylab="Cumlative Work (%)", xlab="Gene Trio Set (m)", main="Cumlative Work for Sets of Gene Trios")
		abline(h=c(cumsum(work)[c(slice_boundaries[2:(length(slice_boundaries)-1)])]/max(cumsum(work))*100))
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

plotCorCoeff <- function(m, idx, col=c("black"), breaks="Scott", ...) {
	col.default <- "grey75"
	
	if(length(col) != length(names(idx))) {
		stop("The number of colours (", length(col) ,") does not match the number of indices list elements (", length(names(idx)) ,").")
	}
	
	# get data from just one triangle and only that which is not NA
	dat <- m[upper.tri(m) & !is.na(m)]
	all.hist <- hist(dat, plot=FALSE, breaks=breaks)
	
	# get the break information from all.hist
	breaks <- all.hist$breaks
	
	# plot the distribution for all m values
	plot.new()
	plot.window(xlim=range(breaks),
			ylim=range(0, all.hist$counts))
	rect(breaks[-length(breaks)], 0,
			breaks[-1], all.hist$counts, col=col.default, ...)
	
	# for each vector of indices in the idx list, superimpose a distribution with the same breaks as all.hist
	for (i in 1:length(idx)) {
		dat <- m[intersect(idx[[i]], which(upper.tri(m)))]
		
		i.hist <- hist(dat, plot=FALSE, breaks=breaks)
		rect(i.hist$breaks[-length(i.hist$breaks)], 0,
				i.hist$breaks[-1], i.hist$counts, col=col[i], ...)
		
	}
	
	axis(1)
	axis(2)
	title(main="Density Distribution of Correlation Coefficients", xlab="Correlation Coefficient", ylab="Frequency")
	
	# Still need to plot the legend
	if( length(names(idx)) >= 1 ) {
		legend(-1, max(all.hist$counts), fill=c(col.default, col), legend=c("All", names(idx)))
	}
}

idx <- function(result) {
	return(result$idx)
}

idxInvert <- function(m, idx) {
	if (class(m)=="numeric" | class(m)=="integer") {
		nNodes <- m
	} else {
		nNodes <- try(nrow(m), silent=TRUE)
		if (class(nNodes) == "try-error" | is.null(nNodes)) {
			cat("ERROR: argument 'm' must be a numeric OR an object on which nrow() can be performed.\n\n", geterrmessage())
			return(FALSE)
		}
	}
	
	return(setdiff(1:(nNodes^2), idx))
}

pcitMemoryRequirement <- function(m, units=c("MB", "bytes", "KB", "GB", "TB"), nCopies=3) {
	double.bytes <- 8
	units <- match.arg(units)
	
	if (class(m)=="numeric" | class(m)=="integer") {
		nNodes <- m
	} else {
		nNodes <- try(nrow(m), silent=TRUE)
		if (class(nNodes) == "try-error" | is.null(nNodes)) {
			cat("ERROR: argument 'm' must be a numeric OR an object on which nrow() can be performed.\n\n", geterrmessage())
			return(FALSE)
		}
	}
	
	switch(units,
			"bytes" = { denominator <- 1024^0 },
			"KB" = { denominator <- 1024^1 },
			"MB" = { denominator <- 1024^2 },
			"GB" = { denominator <- 1024^3 },
			"TB" = { denominator <- 1024^4 }
	)
	
	ram <- nNodes^2*double.bytes*nCopies / denominator
	
	return(list("RAM" = ram, "units" = units))
}

maxMatrixSize <- function(ram, units=c("MB", "bytes", "KB", "GB", "TB"), nCopies=3) {
	units <- match.arg(units)
	switch(units,
			"bytes" = { units <- 0 },
			"KB" = {units <- 1 },
			"MB" = { units <- 2 },
			"GB" = { units <- 3 },
			"TB" = { units <- 4 }
	)
	
	return(floor(sqrt(ram*(1024^units)/(8*nCopies))))
	
}

getEdgeList <- function(m, rm.zero=TRUE) {
	# returns just the upper/lower? triangle
	
	# USAGE:
	# > edgeList <- getEdgeList(adjMatrix)
	# > write.table(edgeList, file='edgeList.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	
	# from=character(length=((nrow(adj)*nrow(adj))-10)/2), to=character(length=((nrow(adj)*nrow(adj))-10)/2), weight=as.single(vector(mode="numeric", length=((nrow(adj)*nrow(adj))-10)/2))
	from_idx <- vector(mode='integer', length=((nrow(m)^2)-nrow(m))/2)
	to_idx <- vector(mode='integer', length=((nrow(m)^2)-nrow(m))/2)
	weight <- vector(mode='numeric', length=((nrow(m)^2)-nrow(m))/2)
	
	result <- .Fortran('getedgelist', mat=as.single(m), nGenes=as.integer(nrow(m)), from_idx=as.integer(from_idx), to_idx=as.integer(to_idx), weight=as.single(weight))
	
	# TODO convert result$from_idx and result$to_idx to character stings according to rownames(adj) and colnames(adj)
	edgeList <- data.frame(From=rownames(m)[result$from_idx], To=colnames(m)[result$to_idx], Weight=as.vector(result$weight))
	#colnames(adjList) <- c('From', 'To', 'Weight')
	
	return(edgeList[edgeList$Weight != 0,])
}
