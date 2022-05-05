###############################################################################
# Makefile for assignment 1, Parallel and Distributed Computing 2022.
###############################################################################

CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm

BIN = matmul

all: $(BIN)

matmul: qsort.c qsort.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
