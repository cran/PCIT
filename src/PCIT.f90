! Function for calculating the local tolerance given the 3 pairwise correlations between X, Y and Z
REAL(KIND=4) FUNCTION tolerance(xy, xz, yz, tolType)
	IMPLICIT  NONE
	
	REAL(KIND=4) :: xy, xz, yz
	INTEGER :: tolType
	
	REAL(KIND=4) :: xy_g_z, xz_g_y, yz_g_x
	xy_g_z = (xy - xz*yz) / sqrt((1-xz**2) * (1-yz**2))
	xz_g_y = (xz - xy*yz) / sqrt((1-xy**2) * (1-yz**2))
	yz_g_x = (yz - xy*xz) / sqrt((1-xy**2) * (1-xz**2))
	
	IF ( tolType == 1 ) THEN
		! mean
		tolerance = (abs(xy_g_z/xy) + abs(xz_g_y/xz) + abs(yz_g_x/yz)) / 3
	ELSE IF ( tolType == 2 ) THEN
		! min
		tolerance = MIN(abs(xy_g_z/xy), abs(xz_g_y/xz), abs(yz_g_x/yz))
	ELSE IF ( tolType == 3 ) THEN
		!max
		tolerance = MAX(abs(xy_g_z/xy), abs(xz_g_y/xz), abs(yz_g_x/yz))
	END IF
	
	RETURN
END FUNCTION tolerance

SUBROUTINE pcit(cor, pcor, nGenes, xVals, nXVals, tolType)
	
	INTEGER :: nGenes, tolType
	REAL(KIND=4) :: cor(nGenes,nGenes), pcor(nGenes,nGenes)
	REAL(KIND=4) :: xy, xz, yz
	REAL(KIND=4) :: tol
	
	INTEGER :: x,y,z  ! loop indices
	INTEGER :: nXVals, xVals(nXVals)
	
	DO i=1, nXVals, +1  ! loop over all supplied xVals
		x = xVals(i)
		DO y = x+1, nGenes-1, +1
			xy = cor(x,y)
			
			DO z = y+1, nGenes, +1
				xz = cor(x,z)
				yz = cor(y,z)
				
				tol = tolerance(xy, xz, yz, tolType)
				
				IF ( abs(xy) < abs(xz*tol) .AND. abs(xy) < abs(yz*tol) )  THEN
					pcor(x,y) = 0    ! upper triangle
					pcor(y,x) = 0    ! lower triangle
				END IF
				IF ( abs(xz) < abs(xy*tol) .AND. abs(xz) < abs(yz*tol) )  THEN
					pcor(x,z) = 0    ! upper triangle
					pcor(z,x) = 0    ! lower triangle
				END IF
				IF ( abs(yz) < abs(xy*tol) .AND. abs(yz) < abs(xz*tol) )  THEN
					pcor(y,z) = 0    ! upper triangle
					pcor(z,y) = 0    ! lower triangle
				END IF
			END DO
		END DO
	END DO
END SUBROUTINE

