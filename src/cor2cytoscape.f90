! Written by Aaron Mcdonough
! modified by Nathan S Watson-Haigh
! only outputs one triangle from the matrix
SUBROUTINE getedgelist(mat, nGenes, from_idx, to_idx, weight)
	INTEGER :: nGenes, removeZero, outputRow=0, from_idx(((nGenes*nGenes)-nGenes)/2), to_idx(((nGenes*nGenes)-nGenes)/2)
	REAL(KIND=4) :: mat(nGenes,nGenes), cutoff, weight(((nGenes*nGenes)-nGenes)/2)
	
	DO i=1,nGenes,+1
		DO j=i+1,nGenes,+1
			outputRow = outputRow+1
			from_idx(outputRow) = i
			to_idx(outputRow) = j
			weight(outputRow) = mat(i,j)
		END DO
	END DO
END SUBROUTINE
