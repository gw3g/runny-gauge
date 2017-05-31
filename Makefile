CC	= gcc
CFLAGS	= -I. -lgsl -lgslcblas -lm
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

temperature:
	gle -o "out/eta_running.png" -tex -d png -dpi 100 "out/plotter/temp.gle" $(NF)
relative:
	gle -o "out/eta-to-NLL.png" -tex -d png -dpi 100 "out/plotter/relative.gle" $(NF)
g4eta:
	gnuplot "out/plotter/eta.gp"

.PHONY: clean

