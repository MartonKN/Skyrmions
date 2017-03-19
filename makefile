# First type LANG=en_US.UTF-8
PNAME=imagTimeGPESolver
CC=icpc
CFLAGS=-I. -openmp
LIBS=-mkl -lfftw3_omp -lfftw3
DEPS=GPE_dataStructures.h GPE_IO.h GPE_allocations.h GPE_solverFunctions.h GPE_spline.h GPE_monitoring.h GPE_random.h dcomplex.h matrix.h
OBJ=GPE_main.o GPE_IO.o GPE_allocations.o GPE_solverFunctions.o GPE_spline.o GPE_monitoring.o GPE_random.o dcomplex.o matrix.o

all: $(PNAME)

%.o: %.cpp $(DEPS)
	$(CC) -O3 -parallel -xT -c -o $@ $< $(CFLAGS) $(LIBS)
# WARNING: the -xT flag means optimalization for Core2Duo processor. On other kinds of processors the program may not run because of this.

$(PNAME): $(OBJ)
	$(CC) -O3 -parallel -xT -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm  *.o *~