MODULE mod_dat

USE numerics
USE mod_cond
USE mod_gnuplot

IMPLICIT NONE 


CONTAINS 


SUBROUTINE init_var(niter,date,dx,dt_sortie)

INTEGER, INTENT(INOUT) :: niter
REAL(rp), INTENT(INOUT) :: date
REAL(rp), INTENT(INOUT) :: dx
REAL(rp), INTENT(INOUT) :: dt_sortie
LOGICAL :: fixe

INTEGER :: i
 
OPEN(UNIT = 50, FILE = 'donnee', STATUS = 'OLD', ACTION = 'READ')

READ(50,*) xmin, xmax
READ(50,*) Nx
READ(50,*)Tfin
READ(50,*) c1, c2
READ(50,*)
READ(50,*) Cas_test
READ(50,*) Cas_topo
READ(50,*) H
READ(50,*) A
READ(50,*) x0
READ(50,*)
READ(50,*) dispersion
READ(50,*) alpha 
READ(50,*)
READ(50,*) lB
READ(50,*) rB
READ(50,*)
READ(50,*) err_moy
CLOSE(50)


OPEN(UNIT = 50, FILE = 'sortie', STATUS = 'OLD', ACTION = 'READ')
READ(50,*) fixe
READ(50,*) 
READ(50,*)
IF(fixe) THEN 
	READ(50,*)dt_sortie
	READ(50,*)
	READ(50,*) 
	READ(50,*) 
	READ(50,*)
	nb_sortie = ceiling(Tfin/(dt_sortie))+1
	ALLOCATE(tab_sortie(nb_sortie))
	DO i = 1, nb_sortie
		tab_sortie(i) = (i-1)*dt_sortie
	END DO 
	tab_sortie(nb_sortie) = Tfin
ELSE
	READ(50,*)
	READ(50,*) 
	READ(50,*) 
	READ(50,*)  nb_sortie
	ALLOCATE(tab_sortie(nb_sortie))
	READ(50,*) tab_sortie(:)
END IF 

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


SUBROUTINE ecrit_frame()

	CALL init_frame_data_names()
	CALL Ecrits_data()
	niter = niter + 1


END SUBROUTINE ecrit_frame


SUBROUTINE ecrit_energie()
REAL(rp) :: Eu_n
REAL(rp) :: Ew_n
	Eu_n = sum(W(1,:)**2)
	Ew_n = sum(W(2,:)**2)
	WRITE(34,*) date, Eu_n/E_0,Ew_n/E_0

END SUBROUTINE ecrit_energie 


SUBROUTINE sol_dat()
REAL(rp) :: pos
INTEGER :: i
REAL(rp), DIMENSION(Nx) :: U
REAL(rp), DIMENSION(Nx) :: W_ex
REAL(rp), DIMENSION(Nx) :: sig_ex
REAL(rp) :: errmax, errL2
      !write(6,*) 'fin de l''algorithme en',niter,' iterations'
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
      OPEN(unit=21,file='out/Solitaire/solutions',status = 'UNKNOWN')
      DO i=1,Nx-1
           write(21,*) X(i),W(1,i),h_temp(i),W(2,i),U(i)!,X(D(1)),0._rp
      END DO 
      write(21,*) X(Nx),W(1,Nx),h_temp(Nx),W(2,Nx),U(Nx)!,X(D(2)),0._rp
      close(21)
	
END SUBROUTINE sol_dat 



END MODULE mod_dat 


