CFLAGS=-Wall -O2
CC=g++
LFLAGS=
INCLUDE=


PROGS = sensitivities bench demo

all:		${PROGS}

sensitivities:	sensitivities.o bs.o

bench:		bench.o bs.o

demo:		demo.o bs.o

bs.o:		bs.cpp
		${CC} -Wall -O2 -DNDEBUG -c bs.cpp

%.o: %.cpp
		${CC} ${CFLAGS} ${INCLUDE} -c $< -o $@

clean:
		rm -f *.o ${PROGS} a.out
