#!/usr/bin/make
# Makefile: makefile for Revolver algorithms
# (c) Juanyun Luo, 2019

objects= main.o lib.o libgraph.o math.o sort_ways.o
all : $(objects)
	g++ -o main -pthread $(objects)

main.o : libgraph.h lib.h math.h sort_ways.hpp
lib.o: libgraph.h lib.h math.h sort_ways.hpp
libgraph.o : libgraph.h math.h lib.h
math.o : libgraph.h lib.h

clean:
	rm *.o
