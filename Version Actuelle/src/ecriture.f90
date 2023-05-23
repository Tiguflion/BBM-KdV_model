MODULE mod_dat

USE numerics
USE mod_cond
USE mod_gnuplot

IMPLICIT NONE 


CONTAINS 


SUBROUTINE init_var(niter,date,dx)

INTEGER, INTENT(INOUT) :: niter
REAL(rp), INTENT(INOUT) :: date
REAL(rp), INTENT(INOUT) :: dx
 
OPEN(UNIT = 50, FILE = 'donnee', STATUS = 'OLD', ACTION = 'READ')

READ(50,*) xmin, xmax
READ(50,*) Nx
READ(50,*)Tfin
READ(50,*) c1, c2
READ(50,*) Cas_test
READ(50,*) Cas_topo
READ(50,*) H
READ(50,*) A
READ(50,*) x0
READ(50,*) dispersion
READ(50,*) alpha 


CLOSE(50)

niter = 0
date = 0.0_rp
dx = (xmax-xmin)/(Nx-1)

CALL init_tab()
CALL Cond_init(W(:,:),X(:),H)

END SUBROUTINE init_var 


SUBROUTINE init_tab()

INTEGER :: i



ALLOCATE(W(2,Nx), W_temp(2,Nx), Flux(2,Nx))
ALLOCATE(X(Nx))
ALLOCATE(S(2,Nx))
ALLOCATE(h_temp(Nx))

W(:,:) = 0._rp
Flux(:,:) = 0._rp
h_temp = 0._rp


DO i = 1, Nx
	X(i) = xmin + (i-1)*(xmax-xmin)/(Nx-1) 
END DO

END SUBROUTINE init_tab


SUBROUTINE ecrit_frame(t_sortie,dt_sortie,date)
	REAL(rp), INTENT(IN) :: dt_sortie
	REAL(rp), INTENT(IN) :: date
	REAL(rp), INTENT(INOUT) :: t_sortie

	if(t_sortie <= date) then 
		niter = niter+1
		CALL init_frame_data_names()
		CALL h_theo(x,date,h_temp)
		CALL Ecrits_data()
		CALL script_sauv_png("Hauteur d'eau à T = ")
	
		t_sortie = t_sortie + dt_sortie
	end if 	

END SUBROUTINE ecrit_frame



SUBROUTINE sol_dat()
REAL(rp) :: pos
INTEGER :: i
REAL(rp), DIMENSION(Nx) :: U
REAL(rp), DIMENSION(Nx) :: W_ex
REAL(rp), DIMENSION(Nx) :: sig_ex
REAL(rp) :: errmax, errL2
      write(6,*) 'fin de l''algorithme en',niter,' iterations'
      write(6,*) 'date de fin',date, Tfin
      
      CALL h_theo(x,0._rp,h_temp)
      CALL u_theo(X,Tfin,U)
      
      OPEN(unit=666, file = 'out/Solitaire/erreur',status = 'UNKNOWN')
      WRITE(666,*) "Nombez d'éléments :", Nx
      WRITE(666,*) 
      WRITE(666,*) "Variable U:"
      WRITE(666,*) "Erreur moyenne  :", sum(abs(h_temp(:)-W(1,:)))/Nx
      WRITE(666,*) "Erreur inf :", maxval(abs(h_temp(:)-W(1,:)))
      WRITE(666,*)
      WRITE(666,*) "Variable W:"
      WRITE(666,*) "Erreur moyenne  :", sum(abs(U(:)-W(2,:)))/Nx
      WRITE(666,*) "Erreur inf :", maxval(abs(U(:)-W(2,:)))
      close(666)
      open(unit=21,file='out/Solitaire/solutions',status = 'UNKNOWN')
      do i=1,Nx
        write(21,*) X(i),W(1,i),h_temp(i),W(2,i),U(i)
        
      enddo

	close(21)
	
END SUBROUTINE sol_dat 

END MODULE mod_dat 


