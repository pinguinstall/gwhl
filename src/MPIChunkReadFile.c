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

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

/**
 * This function takes an file pointer to a text file and creates a return string which is a globally
 * non-overlapping piece of the input file. Furthermore it starts with exactly the beginning of a
 * line and ends exactly at the end of a line. In this way lines of a text input file can be equally
 * distributed among the ranks. Every line of the input file is exactly returned once by one single rank.
 *
 * If the number of lines is smaller than the number of ranks empty strings are returned.
 *
 * This function takes approximately double of the amount of memory needed for the return sting
 * during its runtime. It is very convenient but not very fast.
 *
 * @todo performance and memory improvement
 *
 * @param fp - MPI file pointer
 * @param overlap - should be more than the maximum expected line width
 * @param com - an MPI inter-communicator
 * @return A char pointer to a string which contains all whole (!) lines read by the calling MPI rank
 */
char* GWHL_MPIChunkReadFile(MPI_File *fp, int overlap, MPI_Comm com) {
    int myrank, mpisize, ierr;
    long long int mychunksize;
    char* readchunk = NULL;
    MPI_Comm_rank(com, &myrank);
    MPI_Comm_size(com, &mpisize);
    MPI_Offset filesize;
    MPI_Offset mystart, myend;

    long long int i;

    // get the total size of the file
    ierr = MPI_File_get_size(*fp, &filesize);
    if (ierr) {
        if (myrank == 0) fprintf(stderr, "Couldn't get file size of file\n");
    }

    // delete the EOF of the text file end
    filesize--;

    mychunksize = filesize / mpisize;
    mystart = myrank * mychunksize;
    myend   = mystart + mychunksize - 1;

    if (myrank == mpisize-1) myend = filesize-1;

    // add overlap to the end of everyone's chunk except last rank
    if (myrank != mpisize-1)
        myend += overlap;

    mychunksize =  myend - mystart + 1;

    // allocate memory
    readchunk = malloc( (mychunksize + 1) * sizeof(char) );

    // every rank reads it chunk (including the overlap)
    MPI_File_read_at_all(*fp, mystart, readchunk, mychunksize, MPI_CHAR, MPI_STATUS_IGNORE);
    readchunk[mychunksize] = '\0';

    //search for last "\n" and delete it including the rest of the file
    long long int locend = mychunksize - 1;
    if (myrank != mpisize-1) {
        while(readchunk[locend] != '\n') locend--;
        locend--; // not include the newline sign itself
    }

    // if me is not rank 0 (i.e. I start at the beginning)
    long long int locstart = 0;
    if (myrank != 0) {
        while(readchunk[locstart] != '\n') locstart++;
        locstart++;
    }

//     fprintf(stderr, "%i %lli %lli\n", myrank, locstart, locend);

    long long int mystringsize = locend - locstart + 1;
    char *returnchunk = calloc( (mystringsize + 1), sizeof(char) );

    // copy in the right string (memcpy is not working?)
    for(i = 0; i<mystringsize; i++) {
        returnchunk[i] = readchunk[locstart + i];
    }

    returnchunk[mystringsize] = '\0';

    free(readchunk);

    return returnchunk;
}

/**
 * returns the positions at which to cut for equally distribution of lines of the input file under
 * the MPI ranks.
 */
void GWHL_MPIgetFileChunkBorders(MPI_File *fp, int overlap, MPI_Comm com, long long int *start, long long int *end, long long int *size) {
    int myrank, mpisize, ierr;
    long long int mychunksize;
    char* readbegin = NULL;
    char* readend = NULL;
    MPI_Comm_rank(com, &myrank);
    MPI_Comm_size(com, &mpisize);
    MPI_Offset filesize;
    MPI_Offset mystart, myend;

    long long int i;

    // get the total size of the file
    ierr = MPI_File_get_size(*fp, &filesize);
    if (ierr) {
        if (myrank == 0) fprintf(stderr, "Couldn't get file size of file\n");
    }

    // delete the EOF of the text file end
    filesize--;

    mychunksize = filesize / mpisize;
    mystart = myrank * mychunksize;
    myend   = mystart + mychunksize - 1;

    if (myrank == mpisize-1) myend = filesize-overlap;

    mychunksize =  myend - mystart + 1;

    // set it to mystart, first
    *start = mystart;
    // get some data for searching first '\n'
    // allocate memory
    readbegin = malloc( (overlap + 1) * sizeof(char) );    
    // every rank reads the beginning of its chunk
    MPI_File_read_at_all(*fp, mystart, readbegin, overlap, MPI_CHAR, MPI_STATUS_IGNORE);
    
    readbegin[overlap] = '\0';
  
    // search for first '\n' and returns the next possition
    if(myrank != 0) {
        i = 0;
        while((readbegin[i] != '\n') && (i < overlap)) i++;
        if(readbegin[i] == '\0') {
            *start = -1;
        } else {
            *start = mystart + i + 1;
        }
    } else { // start rank
        *start = 0;
    }
    
    // last rank reads last overlap bytes, kind of stupid, but all have to read
    readend = malloc( (overlap + 1) * sizeof(char) );
    MPI_File_read_at_all(*fp, myend, readend, overlap, MPI_CHAR, MPI_STATUS_IGNORE);
        
    // treading the end
    if (myrank != mpisize-1) {
        readend[overlap] = '\0';       
        i = 0;
        while((readend[i] != '\n') && (i < overlap)) i++;
        if(readend[i] == '\0') {
            *end = -1;
        } else {
            *end = myend + i - 1;
        }
    } else {
        // I am the last mpi process
        *end = filesize;
    }

    *size = *end - *start;

    free(readend);
    free(readbegin);
}
