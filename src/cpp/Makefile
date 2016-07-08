PLATFORM = mac64
VERSION = 650
INC      = /Library/gurobi$(VERSION)/$(PLATFORM)/include/
CC       = gcc
CPP      = g++
CARGS    = -m64 -g -w
CPPLIB   = -L/Library/gurobi$(VERSION)/$(PLATFORM)/lib/ -lgurobi_c++ -lgurobi65

all: missing missing2

missing: missing.cpp
	$(CPP) $(CARGS) -o missing.o missing.cpp -I$(INC) $(CPPLIB) -lm

missing2: missing2.cpp
	$(CPP) $(CARGS) -o missing2.o missing2.cpp -I$(INC) $(CPPLIB) -lm
clean:
	rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp *.dSYM
