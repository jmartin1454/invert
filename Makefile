# Makefile for building C stuff with GSL

CFLAGS=-ansi -Wall -pedantic -g
LDFLAGS=-lgsl -lgslcblas -lm

# uncomment the following lines for a private GSL installation
GSLINCLUDES=/usr/local/include
GSLLIBS=/usr/local/lib
CFLAGS=-Wall -g -I$(GSLINCLUDES)
LDFLAGS=-L$(GSLLIBS) -lgsl -lgslcblas -lm
# you may also have to set LD_LIBRARY_PATH to GSLLIBS


CC=gcc

DEST=invert

$(DEST): $(DEST).o
	$(CC) $(DEST).o $(LDFLAGS) -o $(DEST)

$(DEST).o: $(DEST).c
	$(CC) $(CFLAGS) -c $(DEST).c

clean:
	rm -f *~ *.o core $(DEST)
