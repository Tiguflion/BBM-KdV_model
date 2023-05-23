MODULE mod_correction

USE numerics 
USE mod_algebre
USE mod_operateur 
USE mod_cond
IMPLICIT NONE 

TYPE Penta_diag 
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Mat
END TYPE 


CONTAINS 

SUBROUTINE Correction(W_loc)

	IMPLICIT NONE 
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Mat
	INTEGER :: Nb_Vois
	INTEGER :: Ne ! Nombre d'éléments à corriger 
	INTEGER , DIMENSION(:)  , ALLOCATABLE :: Pos_diag !-- Position des coefficients diagonaux dans A_Creux
	INTEGER , DIMENSION(:)  , ALLOCATABLE :: Vois
	REAL(rp), DIMENSION(:,:)	      :: W_loc
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: A_penta
	REAL(rp), DIMENSION(:)  , ALLOCATABLE :: A_Creux
	REAL(rp), DIMENSION(:)  , ALLOCATABLE :: U ! Vitesse horizontale 
	REAL(rp), DIMENSION(:)  , ALLOCATABLE :: B ! Second membre
	REAL(rp), DIMENSION(:)  , ALLOCATABLE :: A_diag !- Coefficients diagonale de la matrice A

	INTEGER, DIMENSION(:)	, ALLOCATABLE :: vois_debug
	!REAL(rp), DIMENSION(Ne)    :: vec_debug
	!REAL(rp), DIMENSION(Ne,Ne) :: Mat_debug

	INTEGER i,j,k,l

	Ne = size(W_loc(1,:))

	CALL Init_Vec_Cor(Ne,B,U)
	 
	ALLOCATE(A_penta(5,Ne))
	A_Penta(:,:) = 0._rp
	
	CALL rempli_B(Ne,W_loc,dx,B)

	CALL rempli_A_penta(Ne,A_penta)
	
	
	CALL Conditions_Bord(Ne,W_loc,Nb_vois,0._rp,0._rp,0._rp,Pos_diag,vois,B,A_creux) !-- Utiliser des variables globales, obtenue à
	! partir d'un fichier texte pour obtenir a, b, c, tq u0 = au1 + bu2 + c
	!Faire une disjonction de cas, en fonction du bord, afin de permettre d'avoir deux conditions sur deux bords différents.

	CALL res_Chol_Penta(Ne, A_penta, U, B)

	W_loc(1,:) = U(:)

END SUBROUTINE Correction

SUBROUTINE Init_Vec_Cor(Ne,B,U)

IMPLICIT NONE 
	INTEGER,			       INTENT(IN)    :: Ne
	REAL(rp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: U
	REAL(rp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: B	

	INTEGER :: i 

	ALLOCATE(B(Ne),U(Ne))

	B(:) = 0._rp
	U(:) = 0._rp


END SUBROUTINE Init_Vec_Cor


SUBROUTINE rempli_A_penta(Ne,A_penta)
	INTEGER, 		   INTENT(IN)    :: Ne
	REAL(rp), DIMENSION(5,Ne), INTENT(INOUT) :: A_penta

	REAL(rp), DIMENSION(Ne) :: Vec,Vec2
	Vec2(:) = -(alpha**2)
	Vec(:) = 1._rp

	A_penta(:,:) = aD2X_Penta(vec2(:), dx, Ne) + Diag_Penta(vec(:),Ne,0)

	A_penta(5,1) = 0._rp
	A_penta(5,2) = 0._rp
	A_penta(4,1) = 0._rp
	A_penta(2,Ne) = 0._rp
	A_penta(1,Ne-1) = 0._rp
	A_penta(1,Ne) = 0._rp
	
END SUBROUTINE Rempli_A_penta

SUBROUTINE rempli_A_creux(Ne,pos_diag,vois,A_creux)
	INTEGER, INTENT(IN) :: Ne
	INTEGER, DIMENSION(:), INTENT(INOUT) :: pos_diag
	INTEGER, DIMENSION(:), INTENT(INOUT) :: vois
	REAL(rp), DIMENSION(:), INTENT(INOUT) :: A_creux 
	
	INTEGER :: i
	REAL(rp), DIMENSION(Ne) :: Vec
	
!	Vec(:) = 1._rp
!	!--- Si on utilise plusieurs opérateurs, ajouter une fonction de tri en dessous de diag_creux
!	DO i = 1,Ne
!		CALL aD2X_Creux(i,Ne,-(alpha**2),pos_diag,vois,A_creux)
!	END DO 
!	
!	CALL Diag_Creux(Vec,Ne, 0, Pos_diag,Vois, A_creux)
	
END SUBROUTINE rempli_A_creux 
	

SUBROUTINE rempli_B(Ne,W_loc,dx,B)
	INTEGER, INTENT(IN) :: Ne
	REAL(rp), INTENT(IN) :: dx
	REAL(rp), DIMENSION(:,:), INTENT(IN)   :: W_loc
	REAL(rp), DIMENSION(:), INTENT(INOUT) :: B
	
	B = 0._rp
	B = W_loc(1,1:Ne) + alpha*Da(W_loc(2,1:Ne),dx,Ne)
	B(1) = B(1) + alpha*W_loc(2,2)/(2*dx)
	B(Ne) = B(Ne) - alpha*W_loc(2,Ne-1)/(2*dx)

END SUBROUTINE rempli_B 


SUBROUTINE Conditions_Bord(Ne,W_loc,Nb_vois,alp,bet,gam,Pos_diag,vois,B,A_creux) !u_0 = alp*u1 + bet*u2 + gam
!Penser à changer le nom des constantes alp, bet, et gam, on peut faire l'amalgame entre alp et alpha 
	IMPLICIT NONE 
	
	INTEGER,		  INTENT(IN) 	:: Ne
	INTEGER, 		  INTENT(IN)    :: Nb_vois    
	REAL(rp), 		  INTENT(IN)    :: alp
	REAL(rp), 		  INTENT(IN)	:: bet
	REAL(rp), 		  INTENT(IN) 	:: gam
	INTEGER,  DIMENSION(:),   INTENT(IN) 	:: pos_diag
	REAL(rp), DIMENSION(:,:),   INTENT(IN)  :: W_loc
	INTEGER , DIMENSION(:),   INTENT(INOUT)	:: Vois 
	REAL(rp), DIMENSION(:),   INTENT(INOUT) :: B
	REAL(rp), DIMENSION(:),   INTENT(INOUT) :: A_creux 

	INTEGER, DIMENSION(4) :: int_bord 
	INTEGER, DIMENSION(3) :: Coeff	
	INTEGER, DIMENSION(3) :: d_e_f
	!--- Ajouter une variable gloabele Cas_cond, qui détermine la condition aux bords choisie 	
	INTEGER :: i,j,k

	
	Coeff(1) = alp
	Coeff(2) = bet
	Coeff(3) = gam

	d_e_f(1) = -Coeff(1)/(2*dx)
	d_e_f(2) = Coeff(1)/(2*dx)
	d_e_f(3) = -W_loc(2,2)

B(2) = B(2) - Coeff(3)*((alpha/(2*dx))**2)

B(Ne-1) = B(Ne-1) + Coeff(3)*((alpha/(2*dx))**2)

				
END SUBROUTINE Conditions_Bord




	
END MODULE mod_correction
