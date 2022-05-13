###############################################################################
# Makefile for assignment 1, Parallel and Distributed Computing 2022.
###############################################################################

CC = mpicc
CFLAGS = -std=c99 -g
LIBS = -lm

BIN = qsort

all: $(BIN)

qsort: qsort.c qsort.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
