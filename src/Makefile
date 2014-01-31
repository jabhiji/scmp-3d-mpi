# This file can be used to compile and link CARTGEN++ 

# compiler used
CC = mpic++

# optional compile time flags (-O2, -O3 etc)
CFLAGS = -O3

EXE = sc3d.x

## Build targets

# root target (builds the final executable)

$(EXE):	initialize.o \
	streaming.o \
	calc_dPdt.o \
	updateMacro.o \
	updateEquilibrium.o \
	writeMesh.o \
	sc3d.o
	$(CC) initialize.o streaming.o calc_dPdt.o updateMacro.o updateEquilibrium.o writeMesh.o sc3d.o -o $(EXE) -L /Users/jabhiji/MYLIBS/hdf5/lib -lm -lhdf5 -lz -liconv

# compile dependencies

initialize.o: initialize.h initialize.cpp
	$(CC) $(CFLAGS) -c initialize.cpp -o initialize.o

streaming.o: streaming.h streaming.cpp
	$(CC) $(CFLAGS) -c streaming.cpp -o streaming.o

calc_dPdt.o: calc_dPdt.h calc_dPdt.cpp
	$(CC) $(CFLAGS) -c calc_dPdt.cpp -o calc_dPdt.o

updateMacro.o: updateMacro.h updateMacro.cpp
	$(CC) $(CFLAGS) -c updateMacro.cpp -o updateMacro.o

updateEquilibrium.o: updateEquilibrium.h updateEquilibrium.cpp
	$(CC) $(CFLAGS) -c updateEquilibrium.cpp -o updateEquilibrium.o

writeMesh.o: writeMesh.h writeMesh.cpp
	$(CC) $(CFLAGS) -I /Users/jabhiji/MYLIBS/hdf5/include -c writeMesh.cpp -o writeMesh.o

sc3d.o: sc3d.h sc3d.cpp
	$(CC) $(CFLAGS) -c sc3d.cpp -o sc3d.o

clean:
	/bin/rm -f *.o

veryclean:
	/bin/rm -f *.o *.x

wipedata:
	/bin/rm -f ../out/*.xmf ../out/xdmf/*.xmf ../out/hdf5/*.h5
