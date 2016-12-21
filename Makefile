#############################################################
# Makefile for the Elkies package
# Authors: Andreas-Stephan Elsenhans and Joerg Jahnel
# Date: 2007-12-03
#############################################################

CC = g++
LD = g++

CFLAGS = -O2 -funroll-all-loops -falign-functions=64 -Wall -pedantic
LFLAGS = -L/opt/local/lib -lm -lgmp -lgmpxx 

DEPS = fixedpoint.h

OBJ = elkies_alg.o elkies_alg_init.o

VPATH = .

all:    elkies_alg

%.o:%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

elkies_alg: $(OBJ)
	$(LD) -o $@ $^ $(LFLAGS)
#	strip elkies_alg

clean:	
	rm -rf *.o elkies_alg
