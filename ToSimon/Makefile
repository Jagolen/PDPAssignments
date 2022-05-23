CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm

BIN = montecarlo

all: $(BIN)

montecarlo: montecarlo.c montecarlo.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
