TARGET = spectrum
CC = gcc
CFLAGS = -O2 -lm -lfftw3 -g -pg
INC = -I. -I../plot/
LIB = -L.

Q2 = Gibbs
SRCQ2 = q2.c
OBJQ2=$(SRCQ2:.c=.o) \
	../plot/gnuplot_i.o

Q3 = Spectrum
SRCQ3 = q3.c
OBJQ3=$(SRCQ3:.c=.o) \
	../plot/gnuplot_i.o

q2: $(OBJQ2)
	${CC} $(OBJQ2) ${CFLAGS} ${INC} ${LIB} -o $(Q2)

q3: $(OBJQ3)
	${CC} $(OBJQ3) ${CFLAGS} ${INC} ${LIB} -o $(Q3)

$(OBJQ2): $(SRCQ2)
	${CC} $(SRCQ2) ${CFLAGS} ${INC} ${LIB} -lpthread -c

$(OBJQ3): $(SRCQ3)
	${CC} $(SRCQ3) ${CFLAGS} ${INC} ${LIB} -lpthread -c

clean:
	@rm -rf *.o *.out $(Q2) $(Q3)
