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



CALL cpu_time(t1)

CALL init_var(niter,date,dx) !-- Initialise tableaux et variable + conditions initiales

CALL init_animate(dt_sortie) !-- Initialise l'animation


OPEN(UNIT = 34, FILE = 'niv_eau.dat',STATUS= 'UNKNOWN')

niter = 0

DO WHILE (date < Tfin)

	CALL ecrit_frame(t_sortie,dt_sortie,date) !-- Ecriture de l'image pour l'animation

	CALL calcul_pas_de_temps(W,dt)

	date = date + dt

	
	!CALL Calcul_transport()
	CALL Calcul_Burgers()
	
	IF(dispersion) THEN
		 CALL Correction(W(:,2:Nx-1)) 
		 CALL Contrainte(W(:,2:Nx-1)) ! calcul uniquement w_n+1 en fonction de u, peut se faire localement
	END IF
	
	CALL Cond_bord(W) !Calcul de u(:,1), u(:,Nx) et w(:,1) (si c_2>0)


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

CALL Animation()


WRITE(*,*) 'Temps de la résolution : ', t2 - t1

END PROGRAM main 
