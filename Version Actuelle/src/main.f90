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
REAL(rp) :: t_sortie
REAL(rp) :: dt_sortie
INTEGER :: nb_sortie = 10
INTEGER :: compt2 = 0

CALL system ("rm -f out.mp4")

CALL cpu_time(t1)
CALL init_var(niter,date,dt,dx)

CALL init_tab()

CALL Cond_init(W(:,:),X(:),H)

dt_sortie = Tfin/nb_sortie
t_sortie = 0._rp

CALL Init_fichier_data_frame()
CALL init_frame_data_names()



OPEN(UNIT = 34, FILE = 'niv_eau.dat',STATUS= 'UNKNOWN')

compt2 = 0
niter = 0

DO WHILE (date < Tfin)
	
!	h_temp(:) = W(1,:)

	CALL calcul_pas_de_temps(W,dt)

	date = date + dt
	niter = niter+1
	
	CALL Calcul_transport()

		
!	if (dispersion) CALL Correction(W(:,18:45)) ! calcul uniquement u_n+1

	CALL Cond_bord(W) !Calcul de u(:,1), u(:,Nx) et w(:,1) (si c_2>0)

!	if (dispersion) CALL Contrainte() ! calcul uniquement w_n+1 en fonction de u


if(t_sortie < date) then 
	CALL init_frame_data_names()
	CALL h_theo(x,date,h_temp)
	CALL Ecrits_data()
	CALL script_sauv_png("Hauteur d'eau à T = ")
	
	t_sortie = t_sortie + dt_sortie
end if 

write(*,*) 'Temps :',date, ' Temps final :', Tfin 

END DO 


CALL h_theo(x,Tfin,h_temp)
!======================================================================
!
!		SORTIES DANS LE TERMINAL ET AFFICHAGE DES GRAPHIQUES 
!
!
!======================================================================


CALL sol_dat()
CALL CPU_TIME(t2)

!
!CALL Animation()
!
!CALL suppr_fichier_data_frame()
!

WRITE(*,*) 'Temps de la résolution : ', t2 - t1

END PROGRAM main 
