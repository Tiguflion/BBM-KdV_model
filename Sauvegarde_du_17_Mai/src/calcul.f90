MODULE mod_calcul

USE numerics
USE mod_cond 
USE mod_flux
USE mod_vit_top

IMPLICIT NONE 


CONTAINS 



!==================================================
!	FONCTION QUI CALCULE LE PAS DE TEMPS
!
! @Entrée : W(Nx,2), Le vecteur des valeurs de la solution approchée
!
! @Sortie : dt, Le pas de temps

!==================================================
 
SUBROUTINE calcul_pas_de_temps(W,dt) 
      REAL(rp), DIMENSION(:,:), INTENT(in)  :: W
      REAL(rp),                 INTENT(out) :: dt

      INTEGER :: i
      REAL(rp)    :: hl, hr
      REAL(rp)    :: ul, cl, ur, cr
      REAL(rp)    :: aux

      dt = min(dx/abs(c1),dx/abs(c2))
      dt = min(dt,Tfin - date)
	
END SUBROUTINE calcul_pas_de_temps


SUBROUTINE Calcul_Transport()
INTEGER :: i
REAL(rp) :: SolL,SolR

CALL flux_Lax_Friedrichs(W(:,:),Flux(:,:))

	DO i = 2,Nx-1
		!WRITE(*,*) dt/dx*(Flux(:,i) - Flux(:,i-1))
		W(:,i) = W(:,i) - dt/dx*(Flux(:,i) - Flux(:,i-1))
	END DO 
!READ(*,*)
END SUBROUTINE Calcul_transport 




END MODULE mod_calcul 
