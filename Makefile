.SUFFIXES: .f90 .o
#-----------------------
# Basis
#-----------------------
include arch.mk

F90_SRC = aux_mod.f90 prep_mod.f90 prep.f90

OBJS = $(F90_SRC:.f90=.o)

GOAL = prep.x

%.p.f90 : %.f90
	${CPP} ${DFLAGS} $< -o $@

%.o : %.p.f90
	${FC} -c ${FFLAGS} $< -o $@

default :
	${MAKE} ${GOAL}

$(GOAL) : $(OBJS)
	$(LINKER) $(FFLAGS) -o $(GOAL) $(OBJS) $(LIBS)

clean :
	rm -f *.o *.mod *.modmic prep.x
