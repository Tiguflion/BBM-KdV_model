PROGRAM affiche



IMPLICIT NONE 


INTEGER, PARAMETER :: rp = 8

LOGICAL :: fixe 
INTEGER ::  nb_sortie, niter
REAL(rp) :: date, Tfin
REAL(rp), DIMENSION(:), ALLOCATABLE :: tab_sortie
LOGICAL :: video
CHARACTER(12) :: data_name
CHARACTER(13) :: frame_name

date = 0._rp

CALL init_var()


DO niter = 1,nb_sortie
	date = tab_sortie(niter)
	CALL init_frame_data_names()
	CALL script_sauv_png("Hauteur d'eau à T = ")
END DO 

	
CONTAINS 


SUBROUTINE init_var()

IMPLICIT NONE 

INTEGER :: i
REAL(rp) :: dt_sortie
 
OPEN(UNIT = 50, FILE = 'donnee', STATUS = 'OLD', ACTION = 'READ')

READ(50,*) 
READ(50,*) 
READ(50,*)Tfin


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
	nb_sortie = ceiling(Tfin/(dt_sortie)) + 1
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
END SUBROUTINE init_var 


SUBROUTINE init_frame_data_names()
	WRITE(frame_name,'(a,i4.4,a)') 'plot.',niter,'.png'
	WRITE(data_name,'(a,i4.4,a)') 'data',niter,'.dat'
END SUBROUTINE init_frame_data_names

SUBROUTINE script_sauv_png(title)
	CHARACTER(21), INTENT(IN) :: title
	CHARACTER(50) :: var_poubelle

	OPEN(UNIT = 111, FILE = "hauteur_eau.gnu", STATUS = 'UNKNOWN')
		WRITE(111,*) 'set title "',title,date, '"'
		WRITE(111,*) 'set terminal png'
		WRITE(111,*) 'set output "frame/',frame_name,'"'
		WRITE(111,*) 'plot "data/',data_name,'" using 1:2 title "u"w l, "data/'&
		&,data_name, '" using 1:3 title "w" w l, "data/',data_name, '" u 4:5 title "domaine" w l linetype 4 linewidth 3'
	CLOSE(111)
	
	CALL system("gnuplot hauteur_eau.gnu")

END SUBROUTINE script_sauv_png




END PROGRAM affiche
