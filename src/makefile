CC = nvcc
#CFLAGS = -arch=sm_20 -O3 -use_fast_math -DDEBUG -std=c++11
CFLAGS = -arch=sm_35 -O3 -use_fast_math -DDEBUG -std=c++11  
SOURCES = main.cu gpuqt.cu sigma.cu vector.cu hamiltonian.cu model.cu
_OBJ = main.o gpuqt.o sigma.o vector.o hamiltonian.o model.o
HEADERS = sigma.h vector.h gpuqt.h hamiltonian.h model.h


ODIR = obj
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: gpuqt

$(ODIR)/%.o: %.cu $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

gpuqt: $(OBJ)
	$(CC) $^ -o $@

clean:
	rm $(ODIR)/*o gpuqt
