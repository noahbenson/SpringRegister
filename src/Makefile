CC=g++
COPT=-Wall -O3
LIBS=-lm -pthread

default: all

all: springs

clean:
	\rm -rf springs params.o
springs: springs.hh dhash.hh vdw.hh schira.hh main.cc angle_springs.hh params.o
	$(CC) $(COPT) $(LIBS) -o springs main.cc params.o
params.o: params.hh params.cc
	$(CC) $(COPT) -o params.o -c params.cc
