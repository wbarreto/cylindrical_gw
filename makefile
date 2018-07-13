FC=gfortran
PGM     = ads
SRC     = $(PGM).f initial.f gaussj.f\
	frame.f field.f metric.f dynsys.f rk.f constraint.f io.f ko.f

OBJ     = $(PGM).o initial.o gaussj.o\
	frame.o field.o metric.o dynsys.o rk.o constraint.o io.o ko.o

default: $(PGM)

$(PGM): $(OBJ)
	$(FC) -o ads $(OBJ)

#       Clean up executables, object, data and prof files
clean:
	rm -f $(PGM) *.o core *~ *% mon.out fort.* *.dat *.gr *.log

#
# Dependencies.
# Object        Source          Include files
#
#$(PGM).o:       $(PGM).f        
