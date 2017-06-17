#nix makefile for fortran-file	

# put the name of the target program here
TARGET = Gel
SRC = globals.f nanopore.f mrrrr.f rota36.f creador.f cadenas72mr.f fkfun.f chains.definitions.f monomers.definitions.f free_energy.f mainsolver.f saveresults.f fkpsol.f fkpset.f call_kinsol.f savetodisk.f fitspl.f lookup.f splint.f splines.f geom.f graftpoints.f lookup_kai.f kai.f cadenas_mk.f checking_actual_config.f TUNING.f create_protein.f arrayparas.f binary_array.f solver.f presolver.f

# some definitions
SHELL = /bin/bash
FFLAGS= -O3 #-fbounds-check #-O3 #${F90FLAGS} # -O3 -traceback -check all ${F90FLAGS}

#LDFLAGS= -L/export/apps/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s

#LFFLAGS=$(LDFLAGS)

LDFLAGS=-L/home/khuang28/sundials-2.5.0-build/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s

LFLAGS=$(LDFLAGS)

#LFLAGS= -L/home/khuang/sundials-2.5.0-build/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath

#FF = /state/partition1/apps/openmpi-1.6.3/bin/mpif77
#FF = /share/apps/openmpi-1.10.0_no_ib/bin/mpif90
FF = /home/khuang28/mpich/bin/mpif90
#FF = mpif77 #${F90}


all:	$(TARGET)

$(TARGET): $(SRC:.f=.o)
	$(FF) -o $(TARGET) $(SRC:.f=.o) $(LFLAGS) $(LDFLAGS)

.f.o:
	${FF} -c ${FFLAGS}  $(SRC) $(LFLAGS) $(LDFLAGS)

install: all
	cp $(TARGET) ~/bin


clean:	
	@rm -f $(SRC:.f=.o) $(SRC:.f=.d) $(TARGET) *~

realclean: clean
	@rm -f .depend

depend dep:
	@$(FF)  $(CFLAGS) -MM $(SRC) > .depend 

ifeq (.depend, $(wildcard .depend))
include .depend
endif
