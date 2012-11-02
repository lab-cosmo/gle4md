-include ../make.in

SOURCES=$(shell ls *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
DEPENDS=$(SOURCES:.cpp=.d)
TOOLS=gle-fit gle-scaler gle-analyze oned-pimd oned-test
.PHONY: all obj clean libs $(TOOLS)

all: libs $(TOOLS)

clean:
	rm -f $(OBJECTS) $(DEPENDS)

-include $(DEPENDS)
%.d: %.cpp
	$(MAKEDEPEND) $< > $@.$$$$;                      \
           sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;     \
           rm -f $@.$$$$


# generates -- and targets -- executables in another directory
$(TOOLS): % : ../bin/%
$(TOOLS:%=../bin/%): ../bin/%: %.o libglefit.o libcolor.o
	$(CXX) $(CXXFLAGS) -o $@ $< libcolor.o  libglefit.o $(LDFLAGS)
	
../bin/oned-pimd: oned-pimd.o  libglefit.o libcolor.o
	$(CXX) $(CXXFLAGS) -o $@ $< libcolor.o  libglefit.o $(LDFLAGS) -lfftw3
../bin/oned-test: oned-test.o  libglefit.o libcolor.o
	$(CXX) $(CXXFLAGS) -o $@ $< libcolor.o  libglefit.o $(LDFLAGS) -lfftw3
