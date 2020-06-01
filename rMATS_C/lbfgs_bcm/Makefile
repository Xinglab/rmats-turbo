FC = gfortran

FFLAGS = -O

DRIVER1 = driver1.o
DRIVER2 = driver2.o
DRIVER3 = driver3.o

ROUTINES = routines.o

all :  lbfgsb1 lbfgsb2 lbfgsb3 

lbfgsb1 : $(DRIVER1) $(ROUTINES)
	$(FC) $(FFLAGS) $(DRIVER1) $(ROUTINES) -o x.lbfgsb1

lbfgsb2 : $(DRIVER2) $(ROUTINES)
	$(FC) $(FFLAGS) $(DRIVER2) $(ROUTINES) -o x.lbfgsb2

lbfgsb3 : $(DRIVER3) $(ROUTINES)
	$(FC) $(FFLAGS) $(DRIVER3) $(ROUTINES) -o x.lbfgsb3

