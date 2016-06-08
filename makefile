# CC=/usr/bin/g++
CC=g++
CFLAGS = -O3 -std=c++11 -lpthread -march=native
EXEC=Bootest example
all: $(EXEC)

ifeq ($(prof),1)
 CFLAGS+= -pg
endif
ifeq ($(deb),1)
 CFLAGS+= -O0 -DASSERTS -g
endif

.PHONY: test clean


test:
	./Bootest 10000 1 -check
all: $(EXEC)

example: example.cpp
	$(CC) -o $@  $^ $(CFLAGS)

Bootest:  bootest.cpp
	$(CC) -o $@  $^ $(CFLAGS)

%.o: %.cpp %.h
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm Bootest
	rm example
