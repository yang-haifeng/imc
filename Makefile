CXX=g++
#CXX=mpicxx
CFLAGS=-I.
DEPS = Grid.h \
       utils.h \
       typedef.h \
       parameter_input.h \
       io_wrapper.h

OBJ = Grid.o \
      Grid_init.o \
      Grid_utils.o \
      Grid_iteration.o \
      Grid_moveOneCell.o \
      Grid_calcScattering.o \
      Grid_calcEmission.o \
      Grid_calcExtinction.o \
      Grid_saveload.o \
      Grid_Image.o \
      Grid_Integrate_Align.o \
      Grid_muller_Matrix_AlignES.o \
      utils.o \
      parameter_input.o \
      io_wrapper.o \
      # \
      Grid_muller_Matrix.o \
      Grid_Integrate.o \

OBJM = main.o

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

imc : $(OBJM) $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS)

.PHONY: install clean

install : imc
	./install.perl

clean :
	rm -f *.o *~
