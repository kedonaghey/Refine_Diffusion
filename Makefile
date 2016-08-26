EXECS=ser_refine par_refine finalpar_refine
MPICC=mpicc
GCC=gcc
CFLAGS = -g -Og -Wall

all: ${EXECS}

par_refine: par_gridrefine.c
	${MPICC} ${CFLAGS} -o par_refine par_gridrefine.c

finalpar_refine: finalpar_gridrefine.c
	${MPICC} ${CFLAGS} -o finalpar_refine finalpar_gridrefine.c

ser_refine: ser_gridrefine.c
	${GCC} -o ser_refine ser_gridrefine.c -lm

clean:
	rm ${EXECS}

