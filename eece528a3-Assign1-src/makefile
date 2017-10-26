compiler=gcc
mpicompiler = mpicc
basicflag = -O3 -ftree-vectorize -funroll-loops -Wextra -mtune=intel -m64
mpiflag = -Wall -Wpedantic -Wextra -Wconversion -O3 -ftree-vectorize -m64 -march=native -funroll-loops 

all: compiling basic mpitransform

clean:
	rm -rf *.o *.exe

clean_all: clean
	rm -rf *.s *.csv

.PHONY: all clean clean_all

basic:
	$(compiler) $(basicflag) Assign1_Basic_Final.c -o Assign1_Basic_Final.exe 
	
compiling:
	echo "Compiling the programs"


mpitransform:

	$(mpicompiler) $(mpiflag) Assign1_MPI_Final.c -o Assign1_MPI_Final.exe
