PP = g++
CFLAGS = -Wall -std=c++11 -larmadillo
LDFLAGS = -larmadillo -std=c++11

BINDIR = ../bin
HEADERDIR = header

TARGET = $(BINDIR)/main
OBJS = poly.o basis.o 

all: $(TARGET)
	
$(TARGET): $(OBJS)
	$(PP) -o $@ $(LDFLAGS) main.cpp $^

%.o: %.cpp %.h
	$(PP) -c $(CFLAGS) $< -o $@
	

.PHONY: clean
clean:
	rm -f $(OBJS)
	rm -f $(TARGET)
