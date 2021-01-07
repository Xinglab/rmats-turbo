# Set compilers if not set by environment variables
ifndef CC
CC := gcc
endif
ifndef CXX
CXX := g++
endif
ifndef FC
FC := gfortran
endif

FFLAGS := -c
LDFLAGS := -lm -lgfortran -lgsl -lgslcblas -lgomp -lblas -llapack
CFLAGS := -Wall -O2 -msse2 -funroll-loops -fopenmp

# Allow gsl flags to be set by environment variable
ifdef GSL_LDFLAGS
LDFLAGS := $(GSL_LDFLAGS) $(LDFLAGS)
endif
ifdef GSL_CFLAGS
CFLAGS := $(GSL_CFLAGS) $(CFLAGS)
endif

INC_DIR = include
SRC_DIR = src
FOR_DIR = lbfgs_scipy

BIN = rMATSexe
OBJ = util.o myfunc.o cthreadpool.o
FOBJ = lbfgsb.o linpack.o blas.o timer.o
LBFGSB = lbfgsb.o
LINPACK = linpack.o
TIMER = timer.o
ROUTINES = routines.f

all: $(BIN)

$(FOR_DIR)/$(FOBJ):
	cd $(FOR_DIR) && make

rMATSexe: $(SRC_DIR)/main.c $(INC_DIR)/global.h $(INC_DIR)/myfunc.h $(SRC_DIR)/myfunc.c $(INC_DIR)/util.h $(SRC_DIR)/util.c $(INC_DIR)/cthreadpool.h $(SRC_DIR)/cthreadpool.c $(FOR_DIR)/$(LBFGSB) $(FOR_DIR)/$(LINPACK) $(FOR_DIR)/$(TIMER)
	# $(FC) $(FFLAGS) $(FOR_DIR)/$(ROUTINES) -o $(FOR_DIR)/$(FOBJ)
	$(CC) $(CFLAGS) -o $@ $(filter %.o %.c %.cc %.a, $^) $(LDFLAGS)
	make clean

# $(OBJ) :
# 	$(CXX) -c $(CFLAGS) -o $@ $(firstword $(filter %.cpp %.c %.cc, $^) )
# 
# $(FOBJ):
# 	$(FC) $(FFLAGS) -o $@ $(firstword $(filter %.f, $^) )

install:
	cp -f -r $(BIN)  $(INSTALL_PATH)

Rpack:
	make clean
	rm -rf rMATS rMATS*.tar.gz
	cp -r R-package rMATS
	rm -rf rMATS/src/*.o rMATS/src/*.so rMATS/src/*.dll
	rm -rf rMATS/src/*/*.o
	cp -r src/* rMATS/src
	rm rMATS/src/main.c
	cp -r include rMATS/include
	# mkdir rMATS/src/lbfgs_bcm
	# cp -r lbfgs_bcm/routines.f rMATS/src/lbfgs_bcm
	cp lbfgs_bcm/routines.f rMATS/src
	cp wrapper/rMATS_wrapper.h rMATS/include
	cp wrapper/rMATS_wrapper.c rMATS/src
	cp ./LICENSE rMATS
	# cat R-package/src/Makevars|sed '2s/.*/PKGROOT=./' > rMATS/src/Makevars
	cp rMATS/src/Makevars rMATS/src/Makevars.win
	# R CMD build --no-build-vignettes rMATS
	# R CMD build rMATS
	# rm -rf rMATS
	# R CMD check --as-cran rMATS*.tar.gz

Rbuild:
	make Rpack
	R CMD build rMATS
	rm -rf rMATS

Rcheck:
	make Rbuild
	R CMD check --as-cran rMATS*.tar.gz

clean:
	# rm -rf $(OBJ) *.o */*.o */*/*.o *~ */*~ */*/*~
	rm -rf $(OBJ) *.o src/*.o $(FOR_DIR)/*.o
	rm -rf rMATS.Rcheck
