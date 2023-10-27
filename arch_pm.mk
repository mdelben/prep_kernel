DFLAGS=-DNVHPC

FC = ftn 
CC = cc
LINKER = ${FC}
CPP    = /usr/bin/cpp  -C   -nostdinc

#Opt Flags
FFLAGS = -Mfree -acc -mp=multicore,gpu -gpu=cc80  -Mcudalib=cublas,cufft -Mcuda=lineinfo -traceback -Minfo=all,mp,acc -gopt -traceback

