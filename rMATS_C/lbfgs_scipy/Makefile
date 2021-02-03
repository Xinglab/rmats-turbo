# Set compiler if not set by environment variable
ifndef FC
FC := gfortran
endif

FFLAGS = -c -O2

LBFGSB  = lbfgsb.f
LINPACK = linpack.f
TIMER   = timer.f

all :  lbfgsb.o linpack.o timer.o
