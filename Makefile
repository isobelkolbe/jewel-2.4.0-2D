all: jewel-2.4.0-vac jewel-2.4.0-simple jewel-2.4.0-2D

# path to LHAPDF library
LHAPDF_PATH := /Users/isobelkolbe/software/lhapdf6/lib

FC := gfortran
FFLAGS := -g -static -std=legacy -O2


jewel-2.4.0-vac: jewel-2.4.0.o medium-vac.o pythia6425mod-lhapdf6.o meix.o
	$(FC) -o $@ -L$(LHAPDF_PATH) $^ -lLHAPDF -lstdc++

jewel-2.4.0-simple: jewel-2.4.0.o medium-simple.o pythia6425mod-lhapdf6.o meix.o
	$(FC) -o $@ -L$(LHAPDF_PATH) $^ -lLHAPDF -lstdc++


jewel-2.4.0-2D: medium-2D.o jewel-2.4.0.o pythia6425mod-lhapdf6.o meix.o
	$(FC) -o $@ -L$(LHAPDF_PATH) $^ -lLHAPDF -lstdc++


clean:
	rm -f medium-*.o 
	rm -f jewel*.o
	rm -f pythia6425mod.o meix.o
	rm -f *~

.PHONY: all
