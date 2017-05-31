CXX=g++
CFLAGS=-I.
DEPS = Grid.h

OBJ = Grid.o \
      Grid_init.o \
      Grid_utils.o \
      Grid_iteration.o \
      Grid_moveOneCell.o \
      Grid_calc_Scattering.o \
      Grid_saveStokes.o \
      Grid_muller_Matrix.o \
      utils.o \

OBJM = main.o

%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

imc : $(OBJM) $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS)

.PHONY: install clean

install : imc
	./install.perl

clean :
	rm -f *.o *~
