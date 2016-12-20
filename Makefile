#############################################################
# Makefile for the Elkies package
# Authors: Andreas-Stephan Elsenhans and Joerg Jahnel
# Date: 2007-12-03
#############################################################

CC = g++
LD = g++

CFLAGS = -O2 -funroll-all-loops -falign-functions=64 -Wall -pedantic
LFLAGS = -L/opt/local/lib -lm -lgmp -lgmpxx 

DEPS = festkomma.h

OBJ = elkies_allg.o elkies_allg_init.o

VPATH = .

all:    elkies_allg

%.o:%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# Zweistufiges Linken.
# GMP muss statisch gelinkt werden, weil es nur auf dem Hauptrechner gwdu105
# installiert ist. 
elkies_allg: $(OBJ)
	$(LD) -o $@ $^ $(LFLAGS)
#	strip elkies_allg

clean:	
	rm -rf *.o
