MODULE mod_operateur

USE numerics


!===================================================================================================================
!
! Afin de comprendre de que font les fonctions, on utilise les accronymes suivants : 
!
! - D signifie "Dérivé", le chiffre qui suit désigne l'ordre de la dérivée (s'il n'y en a pas, on considère, une dérivation simple).

! - a Correspond à la variable que l'on mets en indice (exemple la fonction Da(u) corresponds à la fonction qui calcule la dérivée de u)
!
! -a1,a2, ... ont le même fonctionnement que a, elles sont utilisées dans le cadre où l'on utilise plusieurs variables. 
!
!
! Dans le cas d'un problème matriciel AX = B, il apparait également :
! 
! - X, l'inconnue du problème, (il apparait généralement dans l'énoncé des fonctionq qui fabriquent A)
!
! 
! Les fonctions qui n'utilisent pas ce jeu de variable, seront commentées
!
!
! Lexique des variables : 
! - Ne, le nombre d'éléments des vecteurs 
! - dx, l'ecart entre deux points 
!
!=================================================================================================================== 



CONTAINS 










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!					PARTIE 1 : FONCTIONS SUR DES VECTEURS
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











FUNCTION Da(a,dx,Ne)   ! ATTENTION : CETTE FONCTION NE PRENDS PAS EN COMPTE LES CONDITIONS AUX BORDS
	INTEGER :: Ne
	REAL(rp):: dx
	REAL(rp), DIMENSION(Ne) :: a
	
	REAL(rp), DIMENSION(Ne) :: Da
	
	INTEGER :: i
	Da = 0._rp

!Gradiant centré
	DO i = 2,Ne-1
		Da(i) = (a(i+1) - a(i-1))/(2._rp*dx)
	END DO 

	Da(1) = 0._rp
	Da(Ne) = 0._rp
END FUNCTION Da










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!					PARTIE 2 : FONCTIONS SUR DES MATRICES
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!==============================================================
! Fonction qui renvoies la matrice identitée pour une taille choisie 
!
!Entrée :
!	  Ne, Le nombre d'éléments de la matrice 
!
!==============================================================





FUNCTION Id(Ne)

	IMPLICIT NONE
	
	INTEGER :: Ne
	
	REAL(rp), DIMENSION(Ne,Ne) :: Id
	
	INTEGER :: i
	

	Id(:,:) = 0._rp
	
	DO i = 1,Ne
		Id(i,i) = 1._rp
	END DO 

END FUNCTION Id



!==============================================================
! Fonction qui remplie une diagonale de matrice 
!
!Entrée :
!	  Vec, un vecteur de valeur
!	  Ne, Le nombre d'éléments de la matrice 
!	  k, la position de la diagonale (positif si le décalage est la droite et négatif à droite, et vaut 0 s'il n'y a pas de décalage)	
!==============================================================


FUNCTION Diag(vec,Ne,k) 

	IMPLICIT NONE 
	INTEGER :: Ne 
	REAL(rp), DIMENSION(:) :: vec
	INTEGER :: k
	REAL(rp), DIMENSION(Ne,Ne) :: Diag
	INTEGER :: i
	
	Diag(:,:) = 0._rp
	IF(k <= 0) THEN
		DO i = 1, Ne+k
			Diag(i,i-k) = vec(i)
		END DO 
	ELSE IF(k > 0) THEN
		DO i = 1, Ne - k
			Diag(i+k,i) = vec(i)
		END DO 
	END IF 

END FUNCTION Diag


FUNCTION DaDX(Vec,dx,Ne)
	INTEGER  :: Ne
	REAL(rp) :: dx ! Dans le cadre ou l'on a un écart constant entre les points d'observations
	REAL(rp), DIMENSION(Ne) :: Vec
	
	REAL(rp), DIMENSION(Ne,Ne) :: DaDX
	
	DaDX(:,:) = 0._rp
	
	DaDX(2:Ne-1,2:Ne-1) = -Diag( (vec(1:Ne-2) + vec(3:Ne))/(4._rp*dx**2) ,Ne-2,0) &
	&+ Diag( vec(3:Ne)/(4._rp*dx**2),Ne-2,-2) + Diag( vec(3:Ne)/(4._rp*dx**2)  , Ne-2,2)


END FUNCTION DaDX



FUNCTION aD2X(vec, dx, Ne)
	INTEGER :: Ne
	REAL(rp) :: dx
	REAL(rp), DIMENSION(Ne) :: Vec
	
	REAL(rp), DIMENSION(Ne,Ne) :: aD2X
	
	aD2X(:,:) = 0._rp
	
	aD2X(1:Ne,1:Ne) = -Diag( 2*vec(1:Ne)/(4._rp*dx**2) ,Ne,0) &
	&+ Diag( vec(1:Ne)/(4._rp*dx**2),Ne,-2) + Diag( vec(1:Ne)/(4._rp*dx**2)  , Ne,2)	

END FUNCTION aD2X





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!					PARTIE 2 : FONCTIONS SUR DES MATRICES CREUSES
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










!==============================================================================================================
!
!					FONCTION DE RECHERCHE DE VOISINS EN FONCTION DES OPERATEURS
!
!==============================================================================================================

SUBROUTINE Vois_DaDX(i,Ne,Nb_Vois,pos_diag)

	INTEGER 	      , INTENT(IN)    :: i
	INTEGER		      , INTENT(IN)    :: Ne
	INTEGER, DIMENSION(:) , INTENT(INOUT) :: pos_diag
	INTEGER 	      , INTENT(INOUT) :: Nb_vois

	INTEGER :: j
	
	Nb_Vois = Nb_vois+2
	DO j = i+1,Ne+1
		pos_diag(j) = pos_diag(j)+2
	END DO 
	
END SUBROUTINE Vois_DaDX


SUBROUTINE Vois_aD2X(i,Ne,Nb_Vois,pos_diag)

	INTEGER 	      , INTENT(IN)    :: i
	INTEGER		      , INTENT(IN)    :: Ne
	INTEGER, DIMENSION(:) , INTENT(INOUT) :: pos_diag
	INTEGER 	      , INTENT(INOUT) :: Nb_vois

	INTEGER :: j
	
	Nb_Vois = Nb_vois+2
	DO j = i+1,Ne+1
		pos_diag(j) = pos_diag(j)+2
	END DO 
	
END SUBROUTINE Vois_aD2X



!==============================================================================================================
!
!					FONCTION D'OPERATEURS CREUX 
!
!==============================================================================================================










!==============================================================================================================
!
!					STOCKAGE MORSE
!
!==============================================================================================================




SUBROUTINE Id_Creux_Morse(Ne,Vec, Pos_diag, A_creux)
	INTEGER, INTENT(IN) :: Ne
	INTEGER(rp), DIMENSION(:) ,    INTENT(IN) :: Pos_diag
	REAL(rp),    DIMENSION(Ne),    INTENT(IN) :: Vec
	REAL(rp),    DIMENSION(:) , INTENT(INOUT) :: A_creux
	
	!Note : Contrairement, à l'appel plein qui est fait par une fonction, on a besoin ici d'instaurer le vecteur vec, car on est dans une subroutine, on ne peux donc pas faire de combinaison linéaire sur l'appel.
	
	DO i = 1,Ne
		A_creux(Pos_diag(i)) = A_creux(pos_diag(i)) + Vec(i)
		! Pas de nécéssitée d'ajouter les voisins ici, car on a dit qu'un point agit forcément sur lui même, lors de l'initialisation du vecteur voisin 
	END DO 
	
END SUBROUTINE Id_Creux_Morse


SUBROUTINE Diag_creux_Morse(Vec,Ne, k, Pos_diag, Vois, A_creux)
	INTEGER, INTENT(IN) :: Ne
	INTEGER, INTENT(IN) :: k
	INTEGER,  DIMENSION(:), INTENT(IN)    :: Pos_diag
	REAL(rp), DIMENSION(:), INTENT(IN)    :: vec
	INTEGER,  DIMENSION(:), INTENT(INOUT) :: Vois 
	REAL(rp), DIMENSION(:), INTENT(INOUT) :: A_creux 
	
	INTEGER :: i,j
	
	DO i = 1,Ne
		DO j = Pos_diag(i), Pos_diag(i+1) - 1
			IF(i+k == Vois(j) .or. (Vois(j) == 0 .and. i+k > 0)) then 
				Vois(j) = i+k
				A_Creux(j) = A_creux(j) + Vec(i)
				EXIT
			END IF 
		END DO 
	END DO 


END SUBROUTINE Diag_Creux_Morse






SUBROUTINE DaDX_Creux_Morse(i,Ne,Vec,pos_diag,vois,A_creux)

	INTEGER, DIMENSION(:),  INTENT(IN) :: pos_diag 
	INTEGER,	        INTENT(IN) :: i
	INTEGER, 	        INTENT(IN) :: Ne
	REAL(rp), DIMENSION(:),	INTENT(IN) :: Vec
	INTEGER, DIMENSION(:),  INTENT(INOUT) :: vois
	REAL(rp), DIMENSION(:), INTENT(INOUT) :: A_creux
	
	
	REAL(rp), DIMENSION(3) :: Coeff
	INTEGER, DIMENSION(3) :: Pos
	
	INTEGER k,j


	Coeff(:) = 0._rp
	
	Pos(2) = i-2
	Pos(1) = i
	Pos(3) = i+2
	
	!-----  Faire une séléction de cas, en fonction de la valeur de i 
		
	Coeff(2) = Vec(i-1)/(4._rp*dx**2)
	Coeff(1) = -(Vec(i+1) + Vec(i-1))/(4._rp*dx**2)
	Coeff(3) = Vec(i+1)/(4._rp*dx**2)



	!----- Shémas type de remplissage creux
	

	k = 1

	
	DO j = Pos_diag(i), Pos_diag(i+1)-1
		if((vois(j) == Pos(k) .or. vois(j) == 0) .and. (k <= 3 .and. k > 0) ) then 
			if(vois(j) == 0) then 
				A_Creux(j) = Coeff(k)
			else 
				A_creux(j) = A_creux(j) + Coeff(k)
			end if 
			vois(j) = pos(k)
			
			IF(i >= Ne-1 .and. k == 2) THEN 
				k = 4
			ELSE IF(i <= 2) THEN 
				k = k+2
			ELSE
				k = k+1
			END IF 

		end if 
	END DO 
 

END SUBROUTINE DaDX_Creux_Morse




SUBROUTINE aD2X_creux_Morse(i,Ne,vec,pos_diag,vois,A_creux)

	INTEGER, DIMENSION(:),  INTENT(IN) :: pos_diag 
	INTEGER,	        INTENT(IN) :: i
	INTEGER, 	        INTENT(IN) :: Ne
	REAL(rp), 		INTENT(IN) :: vec
	INTEGER, DIMENSION(:),  INTENT(INOUT) :: vois
	REAL(rp), DIMENSION(:), INTENT(INOUT) :: A_creux
	
	
	REAL(rp), DIMENSION(3) :: Coeff
	INTEGER, DIMENSION(3) :: Pos
	
	INTEGER k,j


	Coeff(:) = 0._rp
	
	Pos(1) = i
	Pos(2) = i-2
	Pos(3) = i+2
	
	Coeff(1) = -(vec)/(2._rp*(dx**2))
	Coeff(2) = vec/(4._rp*(dx**2))
	Coeff(3) = vec/(4._rp*(dx**2))


	!----- Shémas type de remplissage creux
	k = 1
	DO j = Pos_diag(i), Pos_diag(i+1)-1
		IF((vois(j) == Pos(k) .or. vois(j) == 0) .and. k <= 3) THEN 
			IF((i == 1 .OR. i == Ne) .and. k == 1) THEN !-- Calculer théoriquement le cas Ne
				Coeff(1) = Coeff(1)/2._rp
			END IF 
			IF(vois(j) == 0) THEN 
				A_Creux(j) = Coeff(k)
			ELSE 
				A_creux(j) = A_creux(j) + Coeff(k)
			END IF 
			vois(j) = pos(k)
			IF(i >= Ne-1 .and. k == 2) THEN 
				k = 4
			ELSE IF(i <= 2) THEN 
				k = k+2
			ELSE
				k = k+1
			END IF 

		end if 
	END DO 
 

END SUBROUTINE aD2X_creux_Morse





!==============================================================================================================
!
!					STOCKAGE PENTA-DIAGONAL
!
!==============================================================================================================


FUNCTION Diag_Penta(vec,Ne,k) 

	IMPLICIT NONE 
	INTEGER :: Ne 
	REAL(rp), DIMENSION(:) :: vec
	INTEGER :: k
	REAL(rp), DIMENSION(5,Ne) :: Diag_Penta
	INTEGER :: i
	
	Diag_Penta(:,:) = 0._rp
	
IF(k <= 0) THEN 
	Diag_penta(3+k,1:Ne+k) = vec(1:Ne+k)
ELSE
	Diag_penta(3+k,1+k:Ne) = vec(1+k:Ne)
END IF  	

END FUNCTION Diag_penta


FUNCTION aD2X_Penta(vec, dx, Ne)
	INTEGER :: Ne
	REAL(rp) :: dx
	REAL(rp), DIMENSION(Ne) :: Vec
	
	REAL(rp), DIMENSION(5,Ne) :: aD2X_Penta
	
	aD2X_Penta(:,:) = 0._rp
	
	aD2X_Penta(:,:) = -Diag_Penta( 2*vec(1:Ne)/(4._rp*dx**2) ,Ne,0) &
	&+ Diag_Penta( vec(1:Ne)/(4._rp*dx**2),Ne,-2) + Diag_Penta( vec(1:Ne)/(4._rp*dx**2)  , Ne,2)
	
END FUNCTION aD2X_Penta






END MODULE mod_operateur 
