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
OBJ = util.o myfunc.o
FOBJ = lbfgsb.o linpack.o blas.o timer.o
LBFGSB = lbfgsb.o
LINPACK = linpack.o
TIMER = timer.o

all: $(BIN)

$(FOR_DIR)/$(FOBJ):
	cd $(FOR_DIR) && make

rMATSexe: $(SRC_DIR)/main.c $(INC_DIR)/global.h $(INC_DIR)/myfunc.h $(SRC_DIR)/myfunc.c $(INC_DIR)/util.h $(SRC_DIR)/util.c $(FOR_DIR)/$(LBFGSB) $(FOR_DIR)/$(LINPACK) $(FOR_DIR)/$(TIMER)
	$(CC) $(CFLAGS) -o $@ $(filter %.o %.c %.cc %.a, $^) $(LDFLAGS)
	make clean

install:
	cp -f -r $(BIN)  $(INSTALL_PATH)

clean:
	rm -rf $(OBJ) *.o src/*.o $(FOR_DIR)/*.o
