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
	h_temp(:) = W(1,:)
	

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
	     W(1,1) = 0._rp !W(1,2)
	     W(2,1) = 0._rp !W(2,2)
	     
	     W(1,Ne) = 0._rp !W(1,Ne-1)
	     W(2,Ne) = 0._rp !W(2,Ne-1)
	   	
END SUBROUTINE Cond_bord


!=========================================================================
!
!	FONCTION DE CONTRAINTES POUR LA DISPERSION
!
!=========================================================================

SUBROUTINE Contrainte(W)
	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: W
	INTEGER :: Ne, i
	
	Ne = size(W(2,:)) !-- Permet une utilisation locale de la contrainte 
	
	DO i = 2,Ne-1
		W(2,i) = -alpha*(W(1,i+1) - W(1,i-1))/(2*dx) 
	END DO 

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

	real(rp),parameter::&
		x_L=0.75_rp&
		,x_R=1.25_rp&
		,sig=1.e-1_rp

sol(:)=0._rp
DO i = 1,Nx
	if (X(i) > x_L .and. X(i) < x_R) sol(i) = 1._rp
!	sol(i) = max(exp(-((x(i)-x_L)/sig)**2)+exp(-((x(i)-x_R)/sig)**2),sol(i))
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
