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
#include <stdio.h>
#include <string.h>
#include <math.h>
#define GWHLHELPERBUFSIZE 100

void GWHL_printMatrix(FILE *stream, const char *fmtString, const double *mat, const long long unsigned int sizeM, const long long unsigned int sizeN, int endWithNewline){
  char last[GWHLHELPERBUFSIZE];
  char first[GWHLHELPERBUFSIZE];
  long long unsigned int i;
  long long unsigned int j;
  
  snprintf(last, sizeof(last), "%s%s", fmtString, "\n");
  snprintf(first, sizeof(first), "%s%s", fmtString, " ");
  
  for(i = 0; i < sizeM; i++){
      for(j = 0; j < sizeN; j++){
	if(j < sizeN-1){
	  fprintf(stream, first, mat[i*sizeN+j]);
	} else {
	  if(i < sizeM-1){
	    fprintf(stream, last, mat[i*sizeN+j]);
	  } else if (i == sizeM -1 && endWithNewline){
	    fprintf(stream, last, mat[i*sizeN+j]);
	  } else {
	    fprintf(stream, first, mat[i*sizeN+j]);
	  }
	}
      }
  }
}


void GWHL_printVector(FILE *stream, const char *fmtString, const double *vec, const long long unsigned int length, int endWithNewline){
    char last[GWHLHELPERBUFSIZE];
    char first[GWHLHELPERBUFSIZE];
    long long unsigned int i;
    snprintf(last, sizeof(last), "%s%s", fmtString, "\n");
    snprintf(first, sizeof(first), "%s%s", fmtString, " ");
    for(i = 0; i < length - 1; i++){
	    fprintf(stream, first, vec[i]);
    }
    
    if(endWithNewline){
	if(endWithNewline > 1){
	    fprintf(stream, last, vec[i]);
	    for(i = 0; i < endWithNewline - 1; i++){
	      fprintf(stream, "\n");
	    }
	} else {
	    fprintf(stream, last, vec[i]);
	}
    } else {
	fprintf(stream, first, vec[i]);
    }
}


/**
 * get the value from datain which has max(abs(datain[i]))
 * @return datain[i] where max(abs(datain[i])) over all i
 */
double GWHL_Helpers_getMaxAbsFromData(const double *datain, long long int len){
	long long int i;
	double max = 0.0;
	for(i = 0; i < len; i++){
		if(fabs(datain[i]) > max){
			max = datain[i];
		}
	}
	return max;
}



/**
 * get the value from upper right triangle of the NxN matrix
 * matIn which has max(abs(matIn[i])) without the diagonal
 * @return matin[i] where max(abs(datain[i])) over all i in upper right triangle w/o diagonal
 */
double GWHL_Helpers_getMaxAbsFromUpperRT(const double *matIn, long long int N){
	long long int i, j, k = 1;
	double max = 0.0;
	for (i = 0; i < N; i++) {
		for (j = k; j < N; j++) {  
			if(fabs(matIn[N * i + j]) > max){
				max = matIn[N * i + j];
			}
		}
		k++;
	}
	return max;
}
