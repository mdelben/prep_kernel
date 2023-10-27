DFLAGS=-DONE_API

FC = mpif90 
CC = mpicc
LINKER = ${FC}
CPP    = cpp -C -nostdinc -fopenmp

#Opt Flags
FFLAGS = -fc=ifx -free  -O2 -g -traceback -check shape -fp-model precise -no-ipo -align array64byte -fiopenmp -fopenmp-targets=spir64 -qmkl -lmkl_sycl -lsycl -lOpenCL

