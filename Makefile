CC	= gcc
CFLAGS	= -std=c99 -I. -lgsl -lgslcblas -lm
SDIR	= src
ODIR	= build
OUT	= out
TARGET	= bin/eta
NF	= 0

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
compare:
	gnuplot "out/plotter/g4eta.gp"
temperature:
	gnuplot -c "out/plotter/running.gp" $(NF)
relative:
	gnuplot "out/plotter/relative.gp"

.PHONY: clean

