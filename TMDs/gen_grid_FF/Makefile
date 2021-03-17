COMP = gfortran
FFLAGS = -O2    # This makes the code run faster! Optimizes gfortran...
RUNFLAGS = -ffixed-line-length-none -cpp -fbacktrace -fcheck=all
DEBUGFLAGS = -g -Wall -Wextra -Wsurprising -Walign-commons -ffpe-trap=zero,invalid,overflow,underflow
SOURCE = gen_grid.f TMDFF.f DSS14_HESSIAN/fDSS_HESSIAN14.f
LHAPDFLIBS = `lhapdf-config --ldflags`

executables:
	$(COMP) $(FFLAGS) -c $(SOURCE)
	$(COMP) $(LIBPATH) $(LIB) *.o -o gen_grid.out $<  $(LHAPDFLIBS)

clean:
	-rm *.out
#	-rm *mod
	-rm *.o
#	-rm fort.*
	-rm *~
