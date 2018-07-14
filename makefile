FC=gfortran
PGM     = cgw
SRC     = $(PGM).f base.f initial.f gaussj.f dynsys.f metric.f rk.f io.f

OBJ     = $(PGM).o base.o initial.o gaussj.o dynsys.o metric.o rk.o io.f

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
