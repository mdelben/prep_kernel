# module load oneapi/eng-compiler/2023.10.15.002  e4s/23.05/2023.05.15.006-2 hdf5/1.14.2-oneapi-mpich mpich/52.2
DFLAGS=-DONE_API

FC = mpif90 
CC = mpicc
LINKER = ${FC}
CPP    = cpp -C -nostdinc -fopenmp

#Opt Flags
FFLAGS = -fc=ifx -free  -O2 -g -traceback -check shape -fp-model precise -no-ipo -align array64byte -fiopenmp -fopenmp-targets=spir64 -qmkl -lmkl_sycl -lsycl -lOpenCL

