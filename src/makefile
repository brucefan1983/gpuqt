# use these to build a CPU version
#CC = g++
#CFLAGS = -O3 -std=c++11 -x c++ -DDEBUG -DCPU_ONLY

# use these to build a GPU version
CC = nvcc
CFLAGS = -O3 -std=c++11 -arch=sm_35 -use_fast_math -DDEBUG


SOURCES = main.cu lsqt.cu sigma.cu vector.cu hamiltonian.cu model.cu model_general.cu model_lattice.cu
_OBJ = main.o lsqt.o sigma.o vector.o hamiltonian.o model.o model_general.o model_lattice.o
HEADERS = sigma.h vector.h lsqt.h hamiltonian.h model.h


ODIR = obj
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: lsqt

$(ODIR)/%.o: %.cu $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

lsqt: $(OBJ)
	$(CC) $^ -o $@

clean:
	rm $(ODIR)/*o lsqt
