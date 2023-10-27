# module load cray-fftw  cray-hdf5-parallel  craype-accel-amd-gfx90a rocm  cray-python

DFLAGS=-DHIP_API

HIP_INC = -I${ROCM_PATH}/include/ -J/ccs/home/mdelben/frontier_BGW/LOCAL_libs/hipfort/build/include/hipfort/amdgcn/
HIP_LIB = /ccs/home/mdelben/frontier_BGW/LOCAL_libs/hipfort/build/lib/libhipfort-amdgcn.a -L${ROCM_PATH}/lib -lamdhip64  -lhipfft -lhipblas

FC = ftn 
CC = cc
LINKER = ${FC}  
CPP    = cpp  -C -E -P  -nostdinc

#Opt Flags
FFLAGS = -O1 -f free -fopenmp -g -ef -hlist=a -hacc_model=auto_async_none:no_fast_addr:deep_copy ${HIP_INC}
LIBS   = ${HIP_LIB}
