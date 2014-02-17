TARGET = dsp_fft
CC = gcc
CFLAGS = -O2 -lm -lfftw3
INC = -I.
LIB = -L.

SRC = main.c gnuplot_i.c
OBJ=$(SRC:.c=.o)

main: $(OBJ)
	${CC} $(OBJ) ${CFLAGS} ${INC} ${LIB} -o $(TARGET)

$(OBJ): $(SRC)
	${CC} $(SRC) ${CFLAGS} ${INC} ${LIB} -lpthread -c

clean:
	@rm -rf *.o $(TARGET)
