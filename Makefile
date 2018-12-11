MAKE = make

all:
	$(MAKE) src
	$(MAKE) doc
	$(MAKE) test

run:
	./bin/main

src: 
	$(MAKE) -C src

doc:
	$(MAKE) -C doc
.PHONY: clean src doc

clean:
	$(MAKE) -C src clean
	$(MAKE) -C doc clean
	$(MAKE) -C toPlot clean

tests:
	$(MAKE) -C test
