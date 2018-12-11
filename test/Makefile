OBJS = ../src/poly.o ../src/basis.o 

all:
	$(MAKE) -C ../src
	cxxtestgen --error-printer -o tests.cpp mandatory.h
	g++ -std=c++11 -larmadillo -o tests tests.cpp $(OBJS)
