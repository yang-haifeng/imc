CXX=g++
CFLAGS=-I.
DEPS = Grid.h

OBJ = Grid.o \
      Grid_init.o \
      Grid_utils.o \
      Grid_iteration.o \
      Grid_moveOneCell.o \
      Grid_calcScattering.o \
      Grid_calcEmission.o \
      Grid_saveStokes.o \
      Grid_Image.o \
      Grid_Integrate.o \
      Grid_muller_Matrix_AlignES.o \
      utils.o \
      #Grid_muller_Matrix.o \

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
