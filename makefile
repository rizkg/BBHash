# CC=/usr/bin/g++
CXX ?= g++
CFLAGS = -O3 -std=c++11 -lpthread -Wno-unused-result -Wno-format
EXEC=Bootest example example_custom_hash
all: $(EXEC)

ifeq ($(prof),1)
 CFLAGS+= -pg
endif
ifeq ($(deb),1)
 CFLAGS+= -O0 -DASSERTS -g
endif

ifeq ($(sani),1)
 CFLAGS= -std=c++11 -lpthread -fsanitize=address -fno-omit-frame-pointer -O1 -fno-optimize-sibling-calls -g
endif

.PHONY: test clean


test:
	./Bootest 10000 1 2 -check
all: $(EXEC)

example: example.cpp
	$(CXX) -o $@  $^ $(CFLAGS)

example_custom_hash: example_custom_hash.cpp
	$(CXX) -o $@  $^ $(CFLAGS)


Bootest:  bootest.cpp
	$(CXX) -o $@  $^ $(CFLAGS)

%.o: %.cpp %.h
	$(CXX) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm Bootest
	rm example
