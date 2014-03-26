TARGET = spectrum
CC = gcc
CFLAGS = -O2 -lm -lfftw3 -g -pg
INC = -I. -I../plot/
LIB = -L.

spectrum = Spectrum
SRCspectrum = spectrum.c
OBJspectrum=$(SRCspectrum:.c=.o) \
	../plot/gnuplot_i.o

all: spectrum

spectrum: $(OBJspectrum)
	${CC} $(OBJspectrum) ${CFLAGS} ${INC} ${LIB} -o $(spectrum)

$(OBJspectrum): $(SRCspectrum)
	${CC} $(SRCspectrum) ${CFLAGS} ${INC} ${LIB} -lpthread -c

clean:
	@rm -rf *.o *.out $(Q2) $(spectrum)
