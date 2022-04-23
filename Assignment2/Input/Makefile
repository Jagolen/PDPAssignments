
CC = gcc
CFLAGS = -O3
LIBS = -lm

BIN = random_input

all: $(BIN)

random_input: random_input.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
