CC	= gcc
CFLAGS	= -I. -lgsl -lgslcblas -lm
SDIR	= src
ODIR	= build
OUT	= out
TARGET	= bin/eta

# need to make a choice of ONE suffix
SRCEXT	= c

SRC	= $(wildcard $(SDIR)/*.$(SRCEXT))
OBJ	= $(patsubst $(SDIR)/%,$(ODIR)/%,$(SRC:.$(SRCEXT)=.o))
INC	= -I src

$(TARGET): $(OBJ)
	@mkdir -p bin
	@mkdir -p $(OUT)/data
	@echo "Compiling : $(CC) $(CFLAGS) $^ -o $(TARGET)"; $(CC) $(CFLAGS) $^ -o $(TARGET)

$(ODIR)/%.o: $(SDIR)/%.$(SRCEXT)
	@mkdir -p $(ODIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo "Cleaning : $(RM) -r $(ODIR) $(TARGET)"; $(RM) -r $(ODIR) $(TARGET); $(RM) bin/*;

# auxiliary compiles go here:

plotter1:
	gle -o "out/Z2_D4.pdf" -d pdf "out/plotter/Z2_D4.gle"

# test
tester: $(OBJ)
	$(CC) $(CFLAGS) $(INC) build/matrix.o spike/mat_test.c $(BIN) -o bin/tester

.PHONY: clean

