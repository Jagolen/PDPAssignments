CC = mpicc
CFLAGS = -std=c99 -g
LIBS = -lm

BIN = sendresv_replace
BIN2 = contiguous

all: $(BIN) $(BIN2)

sendresv_replace: sendresv_replace.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

contiguous: contiguous.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
	$(RM) $(BIN2)
