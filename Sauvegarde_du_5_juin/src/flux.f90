MODULE mod_flux 

USE numerics
USE mod_cond 
USE mod_vit_top
IMPLICIT NONE 

CONTAINS 


!======================================================================
!	FONCTION QUI CALCULE LE FLUX DE GODUNOV WELL BALANCED DANS TOUT 
!	LE DOMAINE POUR UN CAS SANS FRICTION
!
! @Entrée : W(Nx,2), Le vecteur des valeurs de la solution approchée
!			X(Nx), Le vecteur de discretisation du domaine 
!			g, La constante de gravité (g = 9.81 m.s^{-2})
!			Nx, Le nombre de volumes d'approximations
!
! @Sortie : Flux, Le flux de Godunov Well-Balanced 
!======================================================================

SUBROUTINE flux_Lax_Friedrichs(W,Flux)
	REAL(rp), DIMENSION(:,:), INTENT(IN) :: W
	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: Flux 
	REAL(rp), DIMENSION(2) :: v
	
	INTEGER :: i 
	
	v(1) = c1
	v(2) = c2
	Flux(1:2,1:Nx) = 0._rp
	
	DO i = 1,Nx-1
		Flux(:,i) = v(:)/2._rp*(W(:,i+1)+W(:,i)) + dx/(2*dt)*(W(:,i) - W(:,i+1))
	END DO 

END SUBROUTINE flux_Lax_Friedrichs


SUBROUTINE flux_upwind_scal(W,a,Flux)
	REAL(rp), 		INTENT(IN) :: a
	REAL(rp), DIMENSION(:), INTENT(IN) :: W
	REAL(rp), DIMENSION(:), INTENT(INOUT) :: Flux
	
	INTEGER :: i
	
	
	REAL(rp) :: am, ap
	
	Flux(1:Nx) = 0._rp
	
	am = min(a,0._rp)
	ap = max(a,0._rp)
	
	DO i = 1,Nx-1
		Flux(i) = am*W(i+1) + ap*W(i)
	END DO 
	
END SUBROUTINE flux_upwind_scal



SUBROUTINE flux_upwind_burg(W,Flux)

	REAL(rp), DIMENSION(:), INTENT(IN) :: W
	REAL(rp), DIMENSION(:), INTENT(INOUT) :: Flux
	REAL(rp) :: a
	REAL(rp) :: ap, am
	INTEGER :: i
	

	Flux(1:Nx) = 0._rp
	
	DO i = 1,Nx-1
		a = (W(i) + W(i+1))/2._rp
		am = min(0._rp,a)
		ap = max(0._rp,a)
		Flux(i) = W(i)*ap + W(i+1)*am
	END DO 
		

	
END SUBROUTINE flux_upwind_burg
!=========================================================================
!	FONCTION QUI CALCULE LE SECOND MEMBRE DE L'EQUATION DE SHALLOW WATER
!	POUR UN SCHEMAS DE GODUNOV WELL BALANCED SANS FRICTION
!
!@Entrées : xl, le point d'approximation gauche de l'intervalle
!			xr, le point d'approximation droit de l'intervalle
!			W_l, les valeurs d'approximation sur l'intervalle précédent
!			W_r, les valeurs d'approximation sur l'intervalle suivant 
!
!@Sorties : sol, le vecteur du second membre de l'approximation
!
!=========================================================================

SUBROUTINE second_membre_godunov(xl,xr,W_l,W_r,sol)
	REAL(rp), INTENT(IN) :: xl,xr
	REAL(rp), DIMENSION(4), INTENT(IN) :: W_l,W_r
	REAL(rp), DIMENSION(4), INTENT(OUT) :: sol
	REAL(rp) :: sol1,sol2
	CALL topo(xl,sol1)
	CALL topo(xr,sol2)
	sol(1) = 0.0_rp
	sol(2) = .5_rp*(W_l(1)+W_r(1))*(sol2-sol1)/dx
	sol(3) = 0.0_rp
	sol(4) = 0.0_rp
END SUBROUTINE second_membre_godunov


END MODULE mod_flux
