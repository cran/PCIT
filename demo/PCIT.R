data(PCIT)

# lets see the size and a subset of the correlation matrix
dim(m3)
m3[1:20,1:20]

# apply the PCIT algorithm, forcing the serial implementation - Should take about 45s on a 2.3Ghz CPU
system.time(m3_pcit_serial <- pcit(m3, force.serial=TRUE))

# apply the PCIT algorithm, using the parallel implementation if we can, otherwise fall back to the serial implementation
system.time(m3_pcit <- pcit(m3))

# pcit() doesn't return the ind's in the same order when done in parallel and serial
# check that we got the same answer using both functions
all.equal(m3[m3_pcit_serial$ind], m3[m3_pcit$ind])

# create a copy of the correlation matrix and set all cells to zero which have been identified by PCIT as being non-significant
m3.new <- m3
ind <- m3_pcit_serial$ind
m3.new[ind] <- 0

cc <- clusteringCoefficient(m3.new)
ccp <- clusteringCoefficientPercent(m3.new)
ccp

op <- par(mfrow=c(2,2))
plot(density(m3[upper.tri(m3)]), main="Density Plot of Raw Correlation Coefficients", xlab="Correlation Coefficient")
hist(cc, main="Connectivity Distribution", xlab="Proportion of Connections", ylab="Number of Genes")
hist(cc*length(cc), main="Connectivity Distribution", xlab="Number of Connections", ylab="Number of Genes")
plotCorCoeff(m3, ind)
par(op)
