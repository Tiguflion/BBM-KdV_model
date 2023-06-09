PROGRAM main

USE numerics
USE mod_gnuplot
USE mod_cond 
USE mod_flux
USE mod_calcul
USE mod_dat 
USE mod_correction
USE mod_erreur

IMPLICIT NONE

INTEGER :: i

REAL(rp) :: temp,t1,t2
REAL(rp) :: t_sortie = 0._rp
REAL(rp) :: dt_sortie 
REAL(rp), DIMENSION(:,:), ALLOCATABLE :: W_n



CALL cpu_time(t1)

CALL init_var(niter,date,dx) !-- Initialise tableaux et variable + conditions initiales

IF(video) THEN
	CALL init_animate(dt_sortie) !-- Initialise l'animation
END IF

ALLOCATE(W_n(2,Nx))

E_0 = sum(W(1,:)**2 + W(2,:)**2)


OPEN(UNIT = 34, FILE = 'energie.dat',STATUS= 'UNKNOWN')

niter = 0

DO WHILE (date < Tfin ) !-- La seconde condition vérifie s'il y a du NaN

	IF(video) THEN 
		CALL ecrit_frame(t_sortie,dt_sortie,date) !-- Ecriture de l'image pour l'animation
	END IF 
	
	CALL calcul_pas_de_temps(W,dt)
	
	date = date + dt

	
	!CALL Calcul_transport()
	CALL Calcul_Burgers()
	
	D(:) = Init_D(W)
	
	W_n(:,:) = W(:,:)

	
	IF(dispersion .and. D(1) /= -1) THEN
		DO WHILE(err_tot > err_moy ) 
			W_n(:,D(1):D(2)) = Correction(W_n(:,D(1):D(2))) 
			CALL Contrainte(W_n(1:2,D(1):D(2))) !-- Bug de segmentation lié à cette contrainte
			D(1) = max(1,D(1)-1)
			D(2) = min(Nx,D(2)+1)
			err_tot = nouvelle_erreur(W_n(:,:))
		END DO
	END IF
 
	W(1:2,1:Nx) = W_n(1:2,1:Nx)
	
	CALL Cond_bord(W) !Calcul de u(:,1), u(:,Nx) et w(:,1) (si c_2>0)

	CALL ecrit_energie()
	

write(*,*) 'Temps :',date, ' Temps final :', Tfin 

END DO 


CLOSE(34)
!======================================================================
!
!		SORTIES DANS LE TERMINAL ET AFFICHAGE DES GRAPHIQUES 
!
!
!======================================================================


CALL sol_dat()
!
CALL CPU_TIME(t2)
!
IF(video) THEN 
	CALL Animation()
	!CALL suppr_fichier_data_frame()
END IF 
WRITE(*,*) 'Temps de la résolution : ', t2 - t1

END PROGRAM main 
