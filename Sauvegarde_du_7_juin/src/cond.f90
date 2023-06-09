MODULE mod_cond 

USE numerics
USE mod_vit_top


IMPLICIT NONE 

CONTAINS 


!===========================================================================
!
!	FONCTION  DES CONDITIONS INITIAL DU PROBLÈME
!
! @Entrées : X(Nx), Le vecteur de discretisation du domaine 
!
! @Sorties : W(4,Nx), Le Vecteur des solutions approchées
!			 H, La hauteur d'eau + fond dans le domaine (h_0 dans l'article)
!
!===========================================================================

SUBROUTINE Cond_init(W,X,H)

	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: W
	REAL(rp), DIMENSION(:), INTENT(IN) :: X
	REAL(rp), INTENT(INOUT) :: H
	REAL(rp) :: top !penser à mettre la topologie dans un tableau 
	INTEGER :: i
	
	
	CALL h_zer(W,X,H)
	CALL Contrainte(W(:,1:Nx))

END SUBROUTINE Cond_init		

	
SUBROUTINE h_zer(W,X,H)
	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: W
	REAL(rp), DIMENSION(:), INTENT(IN) :: X
	REAL(rp), INTENT(INOUT) :: H
	INTEGER :: i
	REAL(rp) :: k,eps
	
CALL h_theo(X,0._rp,W(1,:))
END SUBROUTINE h_zer 




SUBROUTINE hu_zer(W,X,H)
	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: W
	REAL(rp), DIMENSION(:), INTENT(IN) :: X
	REAL(rp), INTENT(INOUT) :: H
	REAL(rp), DIMENSION(Nx) :: u_zer
	INTEGER :: i
	
CALL u_theo(x,0._rp,u_zer)

W(2,:) = u_zer(:)

END SUBROUTINE hu_zer 




SUBROUTINE Cond_bord(W)
	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: W
	INTEGER :: i,Ne
	
	Ne = size(W(1,:)) 

	!IF( c2 > 0._rp) THEN 
	     W(1,1) = lB !W(1,2)
	     W(2,1) = -alpha*(lb-W(1,2))/(2*dx)!0._rp ! W(2,2)
	     
	     W(1,Ne) = rB !W(1,Ne-1)
	     W(2,Ne) = -alpha*(W(1,Ne-1)-rB)/(2*dx) !W(2,Ne-1)
	   	
END SUBROUTINE Cond_bord


!=========================================================================
!
!	FONCTION DE CONTRAINTES POUR LA DISPERSION
!
!=========================================================================

SUBROUTINE Contrainte(W_loc)
	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: W_loc
	INTEGER :: Ne, i
	
	Ne = size(W(2,:)) !-- Permet une utilisation locale de la contrainte 
	
	DO i = 2,Ne-1
		W_loc(2,i) = -alpha*(W_loc(1,i+1) - W_loc(1,i-1))/(2*dx) 
	END DO 

	W_loc(2,1) = 0._rp
	W_loc(2,Ne) = 0._rp
END SUBROUTINE Contrainte


!=========================================================================
!
!	FONCTION DE CALCUL DES SOLUTIONS THEORIQUES
!
!=========================================================================


SUBROUTINE h_theo(x,t,sol)
	REAL(rp), DIMENSION(:) :: x
	REAL(rp)	       :: t
	REAL(rp), DIMENSION(Nx) :: sol
	INTEGER :: i	
	REAL(rp) :: k, c
	REAL(rp), PARAMETER :: pi = acos(-1._rp)
	
	REAL(rp) :: deb_sin, fin_sin
	REAL(rp) :: deb_pente, fin_pente
	
deb_sin = 0.25_rp!(Xmax - Xmin)/4._rp
fin_sin = 0.75_rp	!3._rp*(Xmax - Xmin)/4._rp

deb_pente = 1._rp
fin_pente = 2._rp

sol(:)=0._rp

DO i = 1,Nx
	if (Cas_test == 1) then 
		sol(i) = sin(pi/2._rp + pi*(X(i)-Xmin)/(Xmax - Xmin))
		lB = 1._rp
		rB = -1._rp
	else if(Cas_test == 2) then 
		if(X(i) > deb_sin .and. X(i) < fin_sin) sol(i) = sin(pi*2._rp*(X(i) - deb_sin)/(fin_sin - deb_sin))
		lB = 0._rp
		rB = 0._rp
	else
		lB = 1._rp
		rB = 0._rp
		if(X(i) < deb_pente) then 
			sol(i) = lB
		else if(X(i) > fin_pente) then 
			sol(i) = rB
		else
			sol(i) = fin_pente - X(i)
		end if 
	end if 
END DO 
	

END SUBROUTINE h_theo

SUBROUTINE u_theo(x,t,sol)
	REAL(rp), DIMENSION(:), INTENT(IN) :: x
	REAL(rp)	      , INTENT(IN) :: t
	REAL(rp), DIMENSION(:), INTENT(OUT):: sol
	INTEGER :: i	
	REAL(rp) :: k, c
	REAL(rp), DIMENSION(Nx):: h_sol
	
DO i = 1,Nx
	
	sol(i) = 0._rp
	
END DO 

END SUBROUTINE u_theo


END MODULE mod_cond 
