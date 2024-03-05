GSLDIR=../gsl-2.7.1/LIBRARY
CFLAGS=-g -O3 -march=native -fgnu89-inline -pedantic -I $(GSLDIR)/include/ -Warray-bounds -Wall -fopenmp
LDFLAGS= -L $(GSLDIR)/lib/ -lm -lgsl -lgslcblas
FPIC=-fpic
SHARED=-shared

all:
	$(CC) $(CFLAGS) $(FPIC) -c src/LinearAlgebra.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/DataStatistics.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/MathConversions.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/NumberLists.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/Quaternions.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/SphericalCoords.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/CompareFunctions.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/RotationMatrices.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/Helpers.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/VSH.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/LegendreP.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/Cholesky.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/time.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/io.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/bitmasks.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/sorting.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/Statistics.c -I include/
	$(CC) $(CFLAGS) $(FPIC) -c src/SparseLinearAlgebra.c -I include/
	$(CC) $(CFLAGS) $(FPIC) $(SHARED) -o libGWHelper.so *.o $(LDFLAGS)

clean:
	rm -rf *.o *.exe *.so
