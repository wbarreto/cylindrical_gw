FC=gfortran
PGM     = cgw
SRC     = $(PGM).f base.f gaussj.f initial.f

OBJ     = $(PGM).o base.o gaussj.o initial.o

default: $(PGM)

$(PGM): $(OBJ)
	$(FC) -o cgw $(OBJ)

#       Clean up executables, object, data and prof files
clean:
	rm -f $(PGM) *.o core *~ *% mon.out fort.* *.dat *.gr *.log

#
# Dependencies.
# Object        Source          Include files
#
#$(PGM).o:       $(PGM).f        
