/*
 * This file is part of the GWHL (Gravitational Wave Helper Lib).
 *
 * GWHL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GWHL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GWHL.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>

/**
 * slices numWorkUnitsIN, into equal sized chunks giving out the
 * first (inclusive) and last (inclusive) index which the process has to handle.
 * Indices start to count at zero.
 * 
 */
void GWHL_MPIJobDist_getMyWorkChunk(long long int numWorkUnitsIN, long long int *firstIdxOUT, long long int *lastIdxOUT){
  int mpisize;
  int mpirank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  
  long long int perP = numWorkUnitsIN / mpisize;
  *firstIdxOUT = mpirank*perP;
  *lastIdxOUT = (mpirank+1)*perP - 1;
  if(*lastIdxOUT >= numWorkUnitsIN){
    *lastIdxOUT = numWorkUnitsIN - 1;
  }
}
