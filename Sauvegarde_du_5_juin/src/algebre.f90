MODULE mod_algebre

USE numerics

IMPLICIT NONE 

CONTAINS 

   SUBROUTINE gauss_seidel (n, A, X, B)
      INTEGER, INTENT(IN)  ::  n
      INTEGER :: i
      INTEGER :: j
      INTEGER :: k = 0
      REAL(rp), DIMENSION(n,n),  INTENT(IN)     ::  A
      REAL(rp), DIMENSION(n),    INTENT(INOUT)  ::  X
      REAL(rp), DIMENSION(n),    INTENT(IN)     ::  B
		REAL :: ERR = 1
		REAL :: S
		REAL :: R
		ERR = 1
		k = 0

		do while (ERR > 10E-8 .and. k < 10000)
		k = k+1
			ERR = 0._rp
			do i = 1,n
				S = dot_product(A(i,:),X(:))
				R = (B(i) - S)/A(i,i)
				ERR = ERR + R**2
				X(i) = X(i) + R
				ERR = sqrt(ERR)
			end do
			
		end do

!
      RETURN
      END SUBROUTINE gauss_seidel


! Fonction de la décomposition de Cholesky d'une matrice A
! ! @ Variable d'entrées 
! 		n  : Rang de la matrice 
!		A : La matrice à décomposer
! @ Sortie :
!		Cholesky : La solution du problème (Contenant les deux matrices triangulaires)

	FUNCTION Cholesky(n,A)

	IMPLICIT NONE 

	INTEGER :: n
	REAL(rp), DIMENSION(n,n) :: A
	
	REAL(rp), DIMENSION(n,n) :: Cholesky
	REAL(rp), DIMENSION(n,n) :: L
	REAL(rp) :: temp
	INTEGER :: i, j, k
	Cholesky = 0._rp
	L = 0._rp
	DO i = 1, n 
		L(i,i) = sqrt(A(i,i) - dot_product(L(i,1:i-1),L(i,1:i-1)))
		DO j = i+1,n
			L(j,i) = (A(j,i) - dot_product(L(i,1:j-1),L(j,1:j-1)))/L(i,i)
			L(i,j) = L(j,i) !-- Création de la parite inférieure
		END DO
	END DO 
	
	Cholesky = L	
	RETURN 
	END FUNCTION Cholesky	
	
	
		
! Fonction qui résoud le système triangulaire supérieur AX = B
! ! @ Variable d'entrées 
! 				n  : Le rang de la matrice 
!				A : Une matrice triangulaire supérieure
!				B : Le second membre du problème
! @ Sortie :
! 				X : La solution du problème		

	SUBROUTINE sys_trig_sup(n,A,X,B)
	
	IMPLICIT NONE 

	INTEGER, INTENT(IN)  ::  n
	INTEGER :: i
	REAL(rp), DIMENSION(n,n),  INTENT(IN)     ::  A
	REAL(rp), DIMENSION(n),  INTENT(INOUT)  :: X
	REAL(rp), DIMENSION(n),    INTENT(IN)     ::  B	
	
	DO i = n,1,-1
		X(i) = (B(i) - dot_product(A(i+1:n,i),X(i+1:n)))/A(i,i)
	END DO 


	END SUBROUTINE sys_trig_sup
	
	
		
! Fonction qui résoud le système triangulaire inférieur AX = B
! ! @ Variable d'entrées 
! 				n  : Le rang de la matrice 
!				A : Une matrice triangulaire inférieure
!				B : Le second membre du problème
! @ Sortie :
! 				X : La solution du problème			
	
	
	SUBROUTINE sys_trig_inf(n,A,X,B)
	
	IMPLICIT NONE 

	INTEGER, INTENT(IN)  ::  n
	INTEGER :: i
	INTEGER :: j
	REAL :: S
	REAL(rp), DIMENSION(n,n),  INTENT(IN)     ::  A
	REAL(rp), DIMENSION(n),    INTENT(INOUT)  ::  X
	REAL(rp), DIMENSION(n),    INTENT(IN)     ::  B	
	
	
	DO i = 1,n
		X(i) = (B(i) - dot_product(A(i,1:i-1),X(:)))/A(i,i)
	END DO 	
	
	END SUBROUTINE sys_trig_inf
	
! Fonction qui résoud le système matriciel AX = B, à l'aide de la méthode d'une décomposition de Cholesky
! ! @ Variable d'entrées 
! 				n  : Le rang de la matrice A
!				A : La matrice du problème
!				B : Le second membre du problème
! @ Sortie :
! 				X : La solution du problème

      SUBROUTINE res_Chol(n, A, X, B)
      
      IMPLICIT NONE 
      
      INTEGER, INTENT(IN)  ::  n
      INTEGER :: i
      INTEGER :: j
      REAL(rp), DIMENSION(n,n),  INTENT(IN)     ::  A
      REAL(rp), DIMENSION(n),    INTENT(INOUT)  ::  X
      REAL(rp), DIMENSION(n),    INTENT(IN)     ::  B
 
	REAL(rp), DIMENSION(n,n) :: LtL 
	REAL(rp), DIMENSION(n) :: Ly
 	Ly = 0._rp
      LtL = Cholesky(n,A)

      CALL sys_trig_inf(n,LtL,Ly,B)
      
      CALL sys_trig_sup(n,LtL,X,Ly)

!
      RETURN
      END SUBROUTINE res_Chol
      
	
	
	
	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!		OPERATION SUR DES MATRICES CREUSES 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!==============================================================================================================
!
!					STOCKAGE MORSE
!
!==============================================================================================================


		!-- Fonction qui résouds le problème matriciel AU = L , à l'aide de la matrice A_creux, en utilisant la méthode de Gauss-Seidel : 
! ! @ Variable d'entrées 
 
!				N_coeff : Taille du vecteur des poids A_creux
!				Ns : Nombre de points dans le maillage 
!				A_Creux : Vecteur des poids
!				A_cellule : Table des coefficients non nul de la matrice des poids 
!				Coeff_Diag : Tableau des coefficients diagonaux (= Pos_Creux)
!				Sec_mem : Second membre L
!	@ Variable de sortie
!				U : Solution du problème matriciel AU = L


      SUBROUTINE gauss_seidel_creux_Morse (N_Coeff,Ns, A_Creux,Vois,Coeff_Diag,U,B)
      INTEGER, INTENT(IN)  ::  Ns, N_Coeff
      INTEGER :: i
      INTEGER :: j
      INTEGER :: k 
      REAL(rp), DIMENSION(N_Coeff),  INTENT(IN)     ::  A_Creux
      INTEGER, DIMENSION(N_Coeff),  INTENT(IN)     ::  Vois
      INTEGER, DIMENSION(Ns+1), INTENT(IN) :: Coeff_Diag
      REAL(rp), DIMENSION(Ns),    INTENT(IN)::  B
      REAL(rp), DIMENSION(Ns),    INTENT(OUT)  ::  U
		REAL(rp) :: ERR
		REAL(rp) :: S
		REAL(rp) :: R
		k = 0
		ERR = 1._rp
		U(:) = 0._rp
		do while (ERR > 10E-10 .and. k < 20000)
			k = k+1
			ERR = 0._rp
			do i = 1,Ns
				R = B(i)/A_creux(Coeff_diag(i))
				do j = 1, (Coeff_diag(i+1) - Coeff_diag(i))
					R = R - A_Creux(j-1+Coeff_diag(i))*U(Vois(j-1 + Coeff_diag(i)))/A_creux(coeff_diag(i))
				end do 			
			
				! Code original 
				!S = 0._rp
				!
				!do j = 1, (Coeff_diag(i+1) - Coeff_diag(i))
				!	S = S+A_Creux(j-1+Coeff_diag(i))*U(Vois(j-1 + Coeff_diag(i)))
				!end do 
				!
				!R = (B(i) - S)/A_Creux(Coeff_diag(i))
				
				ERR = ERR + R*R
				U(i) = U(i) + R
				ERR = sqrt(ERR)
			end do
			if(k == 20000) then 
				write(*,*) ERR,k
			end if 
		end do

!
      RETURN
      END SUBROUTINE gauss_seidel_creux_Morse
!







  
!==============================================================================================================
!
!					STOCKAGE PENTA-DIAGONAL
!
!==============================================================================================================




SUBROUTINE res_Chol_Penta(n, A, X, B)
      
      IMPLICIT NONE 
      
      INTEGER, INTENT(IN)  ::  n
      INTEGER :: i
      INTEGER :: j
      REAL(rp), DIMENSION(5,n),  INTENT(IN)     ::  A
      REAL(rp), DIMENSION(n),    INTENT(INOUT)  ::  X
      REAL(rp), DIMENSION(n),    INTENT(IN)     ::  B
 
	REAL(rp), DIMENSION(5,n) :: LtL 
	REAL(rp), DIMENSION(n) :: Ly
 	
      LtL = Cholesky_penta(n,A)

      CALL sys_trig_inf_penta(n,LtL,Ly,B)
      
      CALL sys_trig_sup_penta(n,LtL,X,Ly)

!
      RETURN
      END SUBROUTINE res_Chol_Penta




	FUNCTION Cholesky_Penta(n,A)

	IMPLICIT NONE 

	INTEGER :: n
	REAL(rp), DIMENSION(5,n) :: A
	
	REAL(rp), DIMENSION(5,n) :: Cholesky_Penta
	REAL(rp), DIMENSION(5,n) :: L
	REAL(rp) :: temp
	INTEGER :: i, j, k
	Cholesky_Penta = 0._rp
	L = 0._rp
	
	! Dans la prise de note papier
	! A(3, i) = c_i (coeff diag)
	! A(2, i) = d_i (coeff diag + 1)
	! A(1, i) = e_i (coeff diag + 2)
	
	
	L(1,1) = 0._rp
	L(2,1) = 0._rp
	L(1,2) = 0._rp
	
	L(3,1) = sqrt(A(3,1))
	L(4,1) = A(4,1)/sqrt(A(3,1))
	L(5,1) = A(5,1)/sqrt(A(3,1))
	
	L(2,2) = L(4,1)
	L(3,2) = sqrt(A(3,2)- L(4,1)**2)
	L(4,2) = (A(4,2) - L(4,1)*L(5,1))/L(3,2)
	L(5,2) = A(5,2)/L(3,2)
	
	
	L(3,3) = sqrt(A(3,3) - L(5,1)**2 - L(4,2)**2)
	
	DO i = 3,n 

		L(1,i) = L(5,i-2)
		L(2,i) = L(4,i-1)
		L(3,i) = sqrt(A(3,i) - L(5,i-2)**2 - L(4,i-1)**2)
		L(4,i) = (A(4,i) - L(5,i-1)*L(4,i-1))/L(3,i)
		L(5,i) = A(5,i)/L(3,i)
	END DO  


	
	Cholesky_Penta = L	
	
	
	RETURN 
	END FUNCTION Cholesky_Penta		
	
	SUBROUTINE sys_trig_inf_penta(n,A,X,B) 
	
	IMPLICIT NONE 

	INTEGER, INTENT(IN)  ::  n
	INTEGER :: i
	INTEGER :: j
	REAL :: S
	REAL(rp), DIMENSION(5,n)    ::  A
	REAL(rp), DIMENSION(n)      ::  X
	REAL(rp), DIMENSION(n)      ::  B	
	
	X(1) = B(1)/A(3,1)
	X(2) = (B(2) - A(2,2)*X(1))/A(3,2)
	
	DO i = 3,n
		X(i) = (B(i) - A(1,i)*X(i-2) - A(2,i)*X(i-1))/A(3,i)
	END DO 	
	
	END SUBROUTINE sys_trig_inf_penta



	SUBROUTINE sys_trig_sup_penta(n,A,X,B)
	
	IMPLICIT NONE 

	INTEGER, INTENT(IN)  ::  n
	INTEGER :: i
	REAL(rp), DIMENSION(5,n),  INTENT(IN)     ::  A
	REAL(rp), DIMENSION(n),  INTENT(INOUT)  :: X
	REAL(rp), DIMENSION(n),    INTENT(IN)     ::  B	
	
	X(n) = B(n)/A(3,n)
	X(n-1) = (B(n-1) - A(4,n-1)*X(n))/A(3,n-1)
	
	DO i = n-2,1,-1
		X(i) = (B(i) - A(5,i)*X(i+2) - A(4,i)*X(i+1))/A(3,i)
	END DO 


	END SUBROUTINE sys_trig_sup_penta
END MODULE mod_algebre 
