EXECS=par_diff ser_diff par_cart
MPICC=mpicc
GCC=gcc
CFLAGS = -g -Og -Wall -fdiagnostics-color=auto

all: ${EXECS}

ser_diff: ser_diff.c
	${GCC} ${CFLAGS} -o ser_diff ser_diff.c

par_diff: par_diff.c
	${MPICC} ${CFLAGS} -o par_diff par_diff.c

par_cart: par_cartesian.c
	${MPICC} ${CFLAGS} -o par_cart par_cartesian.c

clean:
	rm ${EXECS}


