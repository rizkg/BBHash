# CC=/usr/bin/g++
CXX ?= g++
CFLAGS = -O3 -std=c++11 -lpthread
EXEC=Bootest example example_strings
al: $(EXEC)

ifeq ($(prof),1)
 CFLAGS+= -pg
endif
ifeq ($(deb),1)
 CFLAGS+= -O0 -DASSERTS -g
endif

.PHONY: test clean


test:
	./Bootest 10000 1 2 -check
all: $(EXEC)

example: example.cpp
	$(CXX) -o $@  $^ $(CFLAGS)

example_strings: example_strings.cpp
	$(CXX) -o $@  $^ $(CFLAGS)


Bootest:  bootest.cpp
	$(CXX) -o $@  $^ $(CFLAGS)

%.o: %.cpp %.h
	$(CXX) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm Bootest
	rm example
