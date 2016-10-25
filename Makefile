### Makefile.tpl -
##
## Author: Karl Ljungkvist

CXXFLAGS=-Wall -O3

PROGS=multiply

all: $(PROGS)

multiply: multiply.cc
	$(CXX) $(CXXFLAGS) $^ -o $@


clean:
	$(RM) $(PROGS) *.o
