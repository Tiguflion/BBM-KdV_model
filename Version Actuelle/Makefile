#
#
## Commandes shell
RM = rm -f

#
#
## Compilateur
FC = gfortran
LD = $(FC)

#
#
## Options compilateur
#FFLAGS  = -march=native -O3 -ffree-form -fbacktrace -g
FFLAGS  = -O2 #-ffixed-form -fbacktrace -g -fcheck=all -Wall -ffpe-trap=zero,overflow,invalid
LDFLAGS = -O3
#LDFLAGS = 

#
#
## Executable

SRCDIR = src
EXEC = exe

#
###
# objects Fortran
FOBJS = \
	numerics.o \
	gnuplot.o \
	algebre.o \
	operateur.o \
	vit_et_topo.o \
	cond.o \
	ecriture.o \
	flux.o \
	calcul.o \
	correction.o \
	erreur.o \
	main.o

FMODS = \
	numerics.mod \
	mod_gnuplot.mod \
	mod_algebre.mod \
	mod_operateur.mod \
	mod_vit_top.mod \
	mod_advection.mod \
	mod_cond.mod \
	mod_dat.mod  \
	mod_flux.mod \
	mod_calcul.mod \
	mod_correction.mod \
	mod_erreur.mod \
# sources Fortran
FSRCS = \
	$(SRCDIR)/numerics.f90 \
	$(SRCDIR)/gnuplot.f90 \
	$(SRCDIR)/algebre.f90 \
	$(SRCDIR)/operateur.f90 \
	$(SRCDIR)/vit_et_topo.f90 \
	$(SRCDIR)/advection.f90 \
	$(SRCDIR)/cond.f90 \
	$(SRCDIR)/ecriture.f90 \
	$(SRCDIR)/flux.f90 \
	$(SRCDIR)/calcul.f90 \
	$(SRCDIR)/correction.f90 \
	$(SRCDIR)/erreur.f90 \
	$(SRCDIR)/main.f90 
#
#.o: $(SRCDIR)/%.f90
#	$(FC) $(FFLAGS) $(FSRCS) -c $(FOBJS)
#
#
## Dependances compilation
all: $(FSRCS)
	$(LD) $(LDFLAGS) $(FSRCS) -o $(EXEC)

#
#
# nettoyage
#clean:
	   rm  *.mod

#
numerics.o: $(SRCDIR)/numerics.f90
algebre.o: $(SRCDIR)/algebre.f90 
gnuplot.o: $(SRCDIR)/gnuplot.f90 numerics.o
operateur.o: $(SRCDIR)/operateur.f90 numerics.o
vit_et_topo.o: $(SRCDIR)/vit_et_topo.f90 numerics.o
avection.o: $(SRCDIR)/advection.f90 vit_et_topo.o numerics.o
cond.o: $(SRCDIR)/cond.f90 numerics.o vit_et_topo.o
ecriture.o: $(SRCDIR)/ecriture.f90 numerics.o cond.o
flux.o: $(SRCDIR)/flux.f90 cond.o numerics.o vit_et_topo.o
calcul.o: $(SRCDIR)/calcul.f90 numerics.o cond.o flux.o vit_et_topo.o
correction.o: $(SRCDIR)/correction.f90 numerics.o
erreur.o: $(SRCDIR)/erreur.f90 numerics.o operateur.o correction.o
main.o: $(SRCDIR)/main.f90 flux.o cond.o ecriture.o correction.o gnuplot.o numerics.o 
