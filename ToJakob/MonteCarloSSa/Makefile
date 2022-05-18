CC = sudo mpicc
CFLAGS = -std=c99 -g
LIBS = -lm

BIN = mcssa

all: $(BIN)

mcssa: mcssa.c mcssa.h
	$(CC) $(CFLAGS) -o  $@ $< $(LIBS)

clean:
	sudo $(RM) $(BIN)