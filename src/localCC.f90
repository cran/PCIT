! Written by Nathan S Watson-Haigh
! Based on the following equation: C_i = E_i/[k_i(k_i-1)/2]
!	Where: C_i is the clustering coef of node i
!	       E_i is the sum of the abs adjacencies that exists between direct neighbours of node i
!	       k_i is the number of direct neighbours of node i
!	k_i(k_i - 1)/2 is the maximum possible sum of abs adjacencies that could exist between the direct neighbours of node i.

SUBROUTINE localCC(mat, nGenes, cc, E, k)
	INTEGER :: nGenes, k(nGenes), idx(nGenes)
	REAL(KIND=4) :: mat(nGenes,nGenes), E(nGenes), cc(nGenes)
	LOGICAL :: neighbours(nGenes)

	DO i=1,nGenes,+1
		idx(i) = i
	END DO

	DO i=1,nGenes,+1
		neighbours = (mat(i,1:nGenes) /= 0 .OR. mat(1:nGenes,i) /= 0)	! determine which nodes are neighbours of node i
		k(i) = count(neighbours)										! calculate how many neighbours node i has
		IF (k(i)<2) THEN												! node i has fewer than 2 neighbours, skip this calculation as the E_i will be zero
			CYCLE
		END IF

		! select a subnetwork containing only nodes that are neighbours of node i
		! then sum the adjacencies to get E_i
		E(i) = sum( mat( pack(idx, neighbours), pack(idx, neighbours) ) )

	END DO
	cc = E / ( k * (k-1) )

END SUBROUTINE
