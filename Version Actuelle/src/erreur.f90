MODULE mod_erreur

USE numerics
USE mod_correction
USE mod_operateur

!---Mettre un paramètre d'erreur moyenne dans le fichier données

CONTAINS 


FUNCTION Init_D(W)

	IMPLICIT NONE 
	
	REAL(rp),DIMENSION(2,Nx) :: W
	REAL(rp), DIMENSION(2,Nx) :: W_new
	INTEGER :: i
	
	INTEGER, DIMENSION(2) :: Init_D

	
	Init_D(1) = -1
	W_new(:,:) = W(:,:)
	err_tot = 0._rp
	
	CALL Contrainte(W_new(:,:))
	
	Do i = 1,Nx
		err_tot = abs(W(2,i) - W_new(2,i))/Nx + err_tot
		IF(abs(W(2,i) - W_new(2,i)) > err_moy .and. Init_D(1) == -1) THEN
			Init_D(1) = i
			Init_D(2) = i+3  !-- Le domaine doit avoir une longueur d'au minimum 4 points
		END IF 
		IF(abs(W(2,i) - W_new(2,i)) > err_moy .and. Init_D(1) /= -1 .and. Init_D(2) < i) THEN 
			Init_D(2) = i
		END IF 
	END DO 
	
END FUNCTION Init_D



FUNCTION nouvelle_erreur(W_loc)
	IMPLICIT NONE 
	
	REAL(rp),DIMENSION(2,Nx) :: W_loc
	REAL(rp), DIMENSION(2,Nx) :: W_new
	
	REAL(rp) :: nouvelle_erreur
	INTEGER :: i

	
	W_new(1:2,1:Nx) = W_loc(1:2,1:Nx)
	nouvelle_erreur = 0._rp
	
	CALL Contrainte(W_new(1:2,1:Nx))
	
	DO i = 1,Nx
		nouvelle_erreur = abs(W_loc(2,i) - W_new(2,i))/Nx + nouvelle_erreur
	END DO 	


END FUNCTION nouvelle_erreur


END MODULE mod_erreur
