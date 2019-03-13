MPICC=mpicc
CC=clang
CFLAGS=-g -O3 -fgnu89-inline -pedantic
GSLDIR=/PATH/TO/YOUR/GSL-2.4/...
LDFLAGS=-L $(GSLDIR)/lib/ -lm -lmpi -lgsl -lgslcblas
FPIC=-fpic
SHARED=-shared

all:
	$(MPICC) $(CFLAGS) $(FPIC) -c src/MPIChunkReadFile.c -I include/
	$(MPICC) $(CFLAGS) $(FPIC) -c src/LinearAlgebra.c -I include/
	$(MPICC) $(CFLAGS) $(FPIC) -c src/DataStatistics.c -I include/
	$(MPICC) $(CFLAGS) $(FPIC) -c src/MPIJobDistribution.c -I include/
	$(CC) $(CFLAGS)    $(FPIC) -c src/MathConversions.c -I include/
	$(CC) $(CFLAGS)    $(FPIC) -c src/NumberLists.c -I include/
	$(CC) $(CFLAGS)    $(FPIC) -c src/Quaternions.c -I include/
	$(CC) $(CFLAGS)    $(FPIC) -c src/SphericalCoords.c -I include/
	$(CC) $(CFLAGS)    $(FPIC) -c src/CompareFunctions.c -I include/
	$(CC) $(CFLAGS)    $(FPIC) -c src/RotationMatrices.c -I include/
	$(CC) $(CFLAGS)    $(FPIC) -c src/Helpers.c -I include/
	$(MPICC) $(CFLAGS)    $(FPIC) -c src/VSH.c -I include/
	$(CC) $(CFLAGS)    $(FPIC) -c src/LegendreP.c -I include/
	$(CC) $(CFLAGS)    $(FPIC) -c src/Cholesky.c -I include/
	$(MPICC) $(CFLAGS) $(FPIC) $(SHARED) -o libGWHelper.so *.o $(LDFLAGS)

clean:
	rm -rf *.o *.exe *.so
	rm -rf tests/*.o tests/*.exe

check:
	$(MPICC) $(CFLAGS) -o tests/testMPIChunkReadFile.exe *.o tests/testMPIChunkReadFile.c -I include/ $(LDFLAGS)
	$(MPICC) $(CFLAGS) -o tests/testDataStatistics.exe   *.o tests/testDataStatistics.c   -I include/ $(LDFLAGS)
	$(MPICC) $(CFLAGS) -o tests/testMathConversions.exe  *.o tests/testMathConversions.c  -I include/ $(LDFLAGS)
	$(MPICC) $(CFLAGS) -o tests/testLinearAlgebra.exe  *.o tests/testLinearAlgebra.c  -I include/ $(LDFLAGS)
	$(MPICC) $(CFLAGS) -o tests/testVSH.exe  *.o tests/testVSH.c  -I include/ $(LDFLAGS)
	$(MPICC) $(CFLAGS) -o tests/testVSH2.exe  *.o tests/testVSH2.c  -I include/ $(LDFLAGS)
	$(MPICC) $(CFLAGS) -o tests/testVSH3.exe  *.o tests/testVSH3.c  -I include/ $(LDFLAGS)
	$(MPICC) $(CFLAGS) -o tests/testVSHFit.exe  *.o tests/testVSHFit.c  -I include/ $(LDFLAGS) -L $(GSLDIR)
	$(MPICC) $(CFLAGS) -o tests/testCholesky.exe  *.o tests/testCholesky.c  -I include/ $(LDFLAGS) -L $(GSLDIR)
	$(MPICC) $(CFLAGS) -o tests/testLinAlg.exe  *.o tests/testLinAlg.c  -I include/ $(LDFLAGS) -L $(GSLDIR)
	$(MPICC) $(CFLAGS) -o tests/testVSHDirection.exe  *.o tests/testVSHDirection.c  -I include/ $(LDFLAGS) -L $(GSLDIR)
