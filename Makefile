EXECS=par_refine
MPICC=mpicc
GCC=gcc
CFLAGS = -g -Og -Wall

all: ${EXECS}

par_refine: par_gridrefine.c
	${MPICC} ${CFLAGS} -o par_refine par_gridrefine.c


clean:
	rm ${EXECS}

