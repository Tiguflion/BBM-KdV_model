MODULE numerics 

IMPLICIT NONE



INTEGER, PARAMETER :: rp = 8
INTEGER :: Nx

REAL(rp), PARAMETER :: g = 9.81

REAL(rp), DIMENSION(:,:), ALLOCATABLE :: W, W_temp !-- Vecteur des valeurs approchées 
REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Flux !-- Vecteur des différents Flux 
REAL(rp), DIMENSION(:,:), ALLOCATABLE :: S !-- Vecteur du second membre 
REAL(rp), DIMENSION(:),ALLOCATABLE :: evo_nu !-- On peut généraliserla dimension avec une variable par la suite que l'on mettra dans le jeu de donnée 
REAL(rp), DIMENSION(:), ALLOCATABLE :: h_temp

REAL(rp) :: dx,dt !-- Pas d'espace et de temps
REAL(rp) :: xmin,xmax !-- Extremums du domaine en 1D
REAL(rp), DIMENSION(:), ALLOCATABLE :: X !-- Vecteur de discrétisation du maillage en 1D

REAL(rp) :: date,Tfin !-- Temps actuel, et temps final

REAL(rp) :: Lambda

INTEGER :: Cas_test 
INTEGER :: Cas_topo
INTEGER :: Nobs
INTEGER :: niter
INTEGER :: compt

REAL(rp) :: H !-- Niveau d'eau moyen
REAL(rp) :: A !-- Amplitude de la vague
REAL(rp) :: x0 !-- Centre de la vague initiale

REAL(rp) :: c1
REAL(rp) :: c2

LOGICAL :: video
CHARACTER(12) :: data_name
CHARACTER(13) :: frame_name

LOGICAL :: dispersion
REAL(rp) :: alpha !-- Terme de liuen entre les deux coefficients de KdV

REAL(rp) :: lB
REAL(rp) :: rB

REAL(rp) :: E_0


REAL(rp) :: err_moy ! Erreur moyenne maximale, tolérée dans le domaine.

REAL(rp) :: err_tot
INTEGER, DIMENSION(2) :: D !-- Fronctières de la zone de correction


LOGICAL :: fixe 
INTEGER ::  nb_sortie
REAL(rp), DIMENSION(:), ALLOCATABLE :: tab_sortie


END MODULE numerics 
