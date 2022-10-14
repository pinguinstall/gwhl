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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>

extern void GWHL_printMatrix(FILE *stream, const char *fmtString, const double *mat, const long long unsigned int sizeN, const long long unsigned int sizeM, int endWithNewline);
extern void GWHL_printVector(FILE *stream, const char *fmtString, const double *vec, const long long unsigned int length, int endWithNewline);

/**
 * decompose A using LL-Cholesky
 */
double GWHL_Cholesky_decompLL(const double *A, long long int n, double *R){
  long long d = 0;
  long long j = 0; // col
  long long i = 0; // row
  long long int k;

  double s = 0;
  memcpy(R, A, n*n*sizeof(double));
  
  for (j = 0; j < n; j++) {
    for(i = 0; i <= j; i++){
      s = 0;
      for(k = 0; k < i; k++){
	s += R[k*n+i] * R[k*n+j];
      }
      if(i < j){
	if (R[i*n+i] > 0){
	  R[i*n+j] = (R[i*n+j]-s) / R[i*n+i];
	} else {
	  R[i*n+j] = 0.0;
	}
      } else {
	if(R[j*n+j] > s){
	  R[j*n+j] = sqrt(R[j*n+j] - s);
	} else {
	  R[j*n+j] = 0;
	  d++;
	}
      }
    }
    for(i = j+1; i < n; i++){
      R[i*n+j] = 0.0;
    }
  }
  
  return d;
}


void GWHL_Cholesky_reduceRhsLL(const double *R, const double *b, long long int n, double *yRes){
  long long int i, k;
  double s = 0;
  memcpy(yRes, b, n*sizeof(double));
  
  for(i = 0; i < n; i++){
    if(R[i*n+i] > 0){
      s = 0.0;
      for(k = 0; k < i; k++){
	s += R[k*n+i] * yRes[k];
      }
      yRes[i] = (yRes[i]-s)/R[i*n+i];
    } else {
      yRes[i] = 0.0;
    }
  }
}


void GWHL_Cholesky_solveLL(const double *R, const double *y, long long int n, double *xRes){
  long long int i, k;
  double s = 0;
  
  memcpy(xRes, y, n*sizeof(double));
  
  for(i = n-1; i >= 0; i--){
    if(R[i*n+i] > 0){
      s = 0.0;
      for(k = i+1; k < n; k++){
	s += R[i*n+k] * xRes[k];
      }
      xRes[i] = (xRes[i]-s) / R[i*n+i];
    } else {
      xRes[i] = 0.0;
    }
  }
}


double GWHL_Cholesky_SolveSystemLL(double *A, double *b, long long int nSize, double *xRes){
    double *R = malloc(nSize*nSize*sizeof(double));
    double *y = malloc(nSize*sizeof(double));
    double d;
    d = GWHL_Cholesky_decompLL(A, nSize, R);
    
    GWHL_Cholesky_reduceRhsLL(R, b, nSize, y);
    GWHL_Cholesky_solveLL(R, y, nSize, xRes);
    
    free(R);
    free(y);
    return d;
}
