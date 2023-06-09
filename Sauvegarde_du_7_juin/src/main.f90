PROGRAM main

USE numerics
USE mod_gnuplot
USE mod_cond 
USE mod_flux
USE mod_calcul
USE mod_dat 
USE mod_correction
USE mod_advection

IMPLICIT NONE

INTEGER :: i

REAL(rp) :: temp,t1,t2
REAL(rp) :: t_sortie = 0._rp
REAL(rp) :: dt_sortie 
REAL(rp) :: Eu_n = 0._rp
REAL(rp) :: Ew_n = 0._rp

CALL cpu_time(t1)

CALL init_var(niter,date,dx) !-- Initialise tableaux et variable + conditions initiales

!IF(video) THEN
	CALL init_animate(dt_sortie) !-- Initialise l'animation
!END IF


E_0 = sum(W(1,:)**2 + W(2,:)**2)

OPEN(UNIT = 34, FILE = 'energie.dat',STATUS= 'UNKNOWN')

niter = 0

DO WHILE (date < Tfin .and. W(1,Nx/2) ==  W(1,Nx/2) ) !-- La seconde condition vérifie s'il y a du NaN

	IF(video) THEN 
		CALL ecrit_frame(t_sortie,dt_sortie,date) !-- Ecriture de l'image pour l'animation
	END IF 
	
	CALL calcul_pas_de_temps(W,dt)
	
	date = date + dt

	
	!CALL Calcul_transport()
	CALL Calcul_Burgers()
	
	IF(dispersion) THEN
		CALL Correction(W(:,2:Nx-1)) 
		CALL Contrainte(W(:,2:Nx-1))
	END IF
	
	CALL Cond_bord(W) !Calcul de u(:,1), u(:,Nx) et w(:,1) (si c_2>0)

	Eu_n = 0._rp
	Ew_n = 0._rp
	Eu_n = sum(W(1,:)**2)
	Ew_n = sum(W(2,:)**2)
	WRITE(34,*) date, Eu_n/E_0,Ew_n/E_0


write(*,*) 'Temps :',date, ' Temps final :', Tfin 

END DO 

!======================================================================
!
!		SORTIES DANS LE TERMINAL ET AFFICHAGE DES GRAPHIQUES 
!
!
!======================================================================


CALL sol_dat()

CALL CPU_TIME(t2)

IF(video) THEN 
	CALL Animation()
	!CALL suppr_fichier_data_frame()
END IF 
WRITE(*,*) 'Temps de la résolution : ', t2 - t1

END PROGRAM main 
