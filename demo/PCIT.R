data(PCIT)

m <- m[1:200,1:200]        # just use a small subset of the data

# lets see the size and a subset of the correlation matrix
dim(m)
m[1:10,1:10]

# apply the PCIT algorithm, forcing the serial implementation - Should take about 45s on a 2.3Ghz CPU
system.time(result_serial <- pcit(m, force.serial=TRUE))

# apply the PCIT algorithm, using the parallel implementation if we can, otherwise fall back to the serial implementation
system.time(result <- pcit(m))

# pcit() doesn't return the ind's in the same order when done in parallel and serial
# check that we got the same answer using both functions
all.equal(m[idx(result_serial)], m[idx(result)])

# create a copy of the correlation matrix and set all cells to zero which have been identified by PCIT as being non-significant
m.new <- m
ind <- idx(result_serial)
m.new[ind] <- 0

cc <- clusteringCoefficient(m.new)
ccp <- clusteringCoefficientPercent(m.new)
ccp

op <- par(mfrow=c(3,2))
plot(density(m[upper.tri(m)]), main="Density Plot of Raw Correlation Coefficients", xlab="Correlation Coefficient")
hist(cc, main="Connectivity Distribution", xlab="Proportion of Connections", ylab="Number of Genes")
hist(cc*length(cc), main="Connectivity Distribution", xlab="Number of Connections", ylab="Number of Genes")
plotCorCoeff(m, list("PCIT Significant" = ind), col=c("black"))
plotCorCoeff(m, list("1st 20000"=c(1:20000), "PCIT Significant" = ind), col=c("black", "red"))
plotCorCoeff(m, list("PCIT Significant" = ind, "1st 20000"=c(1:20000)), col=c("black", "red"))
par(op)
