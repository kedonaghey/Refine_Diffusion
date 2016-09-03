CC=gcc
CFLAGS = -O2 -Wall
LDFLAGS=-lm

EXECS=ser_refine
OBJS = ser_gridrefine.o compute_timesteps.o ser_grid_utils.o interface_funcs.o

all: ${EXECS}

ser_gridrefine.c: compute_timesteps.h ser_grid_utils.h interface_funcs.h
compute_timesteps.c: compute_timesteps.h
ser_grid_utils.c: ser_grid_utils.h
interface_funcs.c: interface_funcs.h

ser_refine: ${OBJS}
	gcc ${LDFLAGS} -o ser_refine ${OBJS}

clean:
	${RM} ${OBJS} ${EXECS}



