CC	= gcc
VPATH	=ops:x-sec
CFLAGS= -I. -lgsl -lgslcblas -lm
DEPS	= core.h
ODIR	= obj
_OBJ	= main.o basis.o qcd.o amy_prepInt.o operators.o thermal.o htl.o
OBJ 	= $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

eta: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

