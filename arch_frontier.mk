# module load cray-fftw  cray-hdf5-parallel  craype-accel-amd-gfx90a rocm  cray-python

DFLAGS=-DHIP_API

HIP_INC = -J/ccs/home/mdelben/frontier_BGW/LOCAL_libs/hipfort/build/include/hipfort/amdgcn/ -I${ROCM_PATH}/include/
HIP_LIB = /ccs/home/mdelben/frontier_BGW/LOCAL_libs/hipfort/build/lib/libhipfort-amdgcn.a -L${ROCM_PATH}/lib -lamdhip64  -lhipfft -lhipblas

FC = ftn 
CC = cc
LINKER = ${FC}
CPP    = cpp  -C -E -P  -nostdinc

#Opt Flags
FFLAGS = -f free -h acc -homp -g -ef -hacc_model=auto_async_none:no_fast_addr:no_deep_copy  ${HIP_INC} ${HIP_LIB} 

