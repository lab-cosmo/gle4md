.PHONY: all clean
FC := gfortran -x f95-cpp-input 
LD := gfortran
FFLAGS := -O3 -mtune=native -DUSELIBS
OBJS := $(patsubst %.f90,%.o,$(wildcard *.f90))

all: harmonic.x

harmonic.x: $(OBJS)
	$(LD) $(FFLAGS) -o harmonic.x $(OBJS) -lblas -llapack

md-nose.o: md-tools.o
md-gle.o: md-tools.o
harmonic.o: md-tools.o md-nose.o md-gle.o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm *~ *.o *.mod harmonic.x
