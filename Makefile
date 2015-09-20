CC	= gcc
VPATH	=integrand:x-sec
CFLAGS= -I. -lgsl -lgslcblas -lm
DEPS	= core.h
ODIR	= obj
_OBJ	= main.o basis.o qcd.o amy_prepInt.o operators.o thermal.o htl.o
OBJ 	= $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -g -o $@ $< $(CFLAGS)

eta: $(OBJ)
	gcc -g -o $@ $^ $(CFLAGS)

