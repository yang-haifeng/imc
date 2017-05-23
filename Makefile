CXX=g++
CFLAGS=-I.
DEPS = Grid.h

OBJ = 

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
