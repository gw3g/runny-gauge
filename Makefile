CC	= gcc
CFLAGS	= -I. -lgsl -lgslcblas -lm
SDIR	= src
ODIR	= build
OUT	= out
TARGET	= bin/eta

# need to make a choice of ONE suffix
SRCEXT	= c

# find ALL *.SRCEXT files in ~/src directory
SRC	= $(wildcard $(SDIR)/*.$(SRCEXT))
OBJ	= $(patsubst $(SDIR)/%,$(ODIR)/%,$(SRC:.$(SRCEXT)=.o))
INC	= -I include

$(TARGET): $(OBJ)
	@mkdir -p bin
	@mkdir -p $(OUT)/data
	$(CC) $(CFLAGS) $^ -o $(TARGET)

$(ODIR)/%.o: $(SDIR)/%.$(SRCEXT)
	@mkdir -p $(ODIR)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	$(RM) -r $(ODIR) $(TARGET); $(RM) bin/*;

# auxiliary compiles go here:

temp:
	gle -o "out/eta(temp).pdf" -d pdf "out/plotter/temp.gle"
rel:
	gle -o "out/_REL.pdf" -d pdf "out/plotter/rel.gle"
u:
	gle -o "out/u(eta).pdf" -d pdf "out/plotter/u_eta.gle"
C:
	gle -o "out/C.pdf" -d pdf "out/plotter/C.gle"

# test
tester: $(OBJ)
	$(CC) $(CFLAGS) $(INC) build/matrix.o spike/mat_test.c $(BIN) -o bin/tester

.PHONY: clean

