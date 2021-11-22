CC   = gcc
OBJ  = imageio.o objfunc.o util.o fourier.o framework.o forcehistogram.o forcehistogram_w.o
LIBS =  -lm 

all: fh_new fh_std fh_wang


clean:	
	rm -f $(OBJ) fh_new fh_std fh_wang

fh_new:	imageio.o objfunc.o util.o fourier.o framework.o
	$(CC) imageio.o objfunc.o util.o fourier.o  framework.o -o "fh_new" $(LIBS)

fh_std: imageio.o objfunc.o util.o forcehistogram.o
	$(CC) imageio.o objfunc.o util.o forcehistogram.o -o "fh_std" $(LIBS)

fh_wang: imageio.o objfunc.o util.o forcehistogram_w.o
	$(CC) imageio.o objfunc.o util.o forcehistogram_w.o -o "fh_wang" $(LIBS)

imageio.o: imageio.c
	$(CC) -c imageio.c -o imageio.o $(LIBS)

objfunc.o: objfunc.c
	$(CC) -c objfunc.c -o objfunc.o $(LIBS)

fourier.o: fourier.c
	$(CC) -c fourier.c -o fourier.o $(LIBS)

framework.o: framework.c
	$(CC) -c framework.c -o framework.o $(LIBS)

forcehistogram.o: forcehistogram.c
	$(CC) -c forcehistogram.c -o forcehistogram.o $(LIBS)

forcehistogram_w.o: forcehistogram_w.c
	$(CC) -c forcehistogram_w.c -o forcehistogram_w.o $(LIBS)

util.o: util.c
	$(CC) -c util.c -o util.o $(LIBS)
