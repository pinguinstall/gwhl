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

/**
 * @file   LinearAlgebra.c
 * @author Robin Geyer (robin.geyer@tu-dresden.de)
 * @date   October, 2015
 * @brief  Helper routines for small scale linear algebra
 * @see    LinearAlgebra.h
 *
 * Different helper routines to handle statistics of the least-squares fit
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "Helpers.h"

/**
 * @brief  matR = matA * a
 * 
 * matrix size = 3x3
 * 
 */
void GWHL_matScale3d(double *matA, double a, double *matR) {
    int i;
    for(i = 0; i < 9; i++) {
        matR[i] = matA[i] * a;
    }
}


/**
 * @brief  matR = matA + matB
 * 
 * matrix size = 3x3
 * 
 */
void GWHL_matMatAdd3d(double *matA, double *matB, double *matR) {
    int i;
    for(i = 0; i < 9; i++) {
        matR[i] = matA[i] + matB[i];
    }
}


/**
 * @brief  vecR = matA . vecB
 * 
 * matrix size = 3 x 3, vector size = 3
 * 
 */
void GWHL_matVecMul3d(double *matA, double *vecB, double *vecR) {
    double tmp[3];
    int i;
    
    tmp[0] = (matA[0] * vecB[0]) + (matA[1] * vecB[1]) + (matA[2] * vecB[2]);
    tmp[1] = (matA[3] * vecB[0]) + (matA[4] * vecB[1]) + (matA[5] * vecB[2]);
    tmp[2] = (matA[6] * vecB[0]) + (matA[7] * vecB[1]) + (matA[8] * vecB[2]);
    
    // Wee need this tmp-var to assure that we can use matA OR matB == matR
    for(i = 0; i < 3; i++){
        vecR[i] = tmp[i];
    }
}


void GWHL_matVecMulNd(double *matA, double *vecB, long long int nColsNVec, long long int nRows, double *vecR){
    long long int i, j;
    for(i = 0; i < nRows; i++){
	vecR[i] = 0.0;
	for(j = 0; j < nColsNVec; j++){
	    vecR[i] += matA[i*nRows + j] * vecB[j];
	}
    }
}

/**
 * @brief  result = vecA . vecB
 * @return vector dot product
 * 
 * vector size = 3
 * 
 */
double GWHL_vecVecDot3d(double *vecA, double *vecB) {
    return vecA[0] * vecB[0] + vecA[1] * vecB[1] + vecA[2] * vecB[2];
}


/**
 * @brief  result = vecA . vecB
 * @return vector dot product
 * 
 * vector size = len
 * 
 */
double GWHL_vecVecDotNd(double *vecA, double *vecB, long long int len) {
    long long int i = 0;
    double ret = 0;
    for(i = 0; i < len; i++){
      ret += vecA[i] * vecB[i];
    }
    return ret;
}


/**
 * @brief  matR = matA . matB
 * 
 * matrix size = 3 x 3
 * 
 */
void GWHL_matMatMul3d(double *matA, double *matB, double *matR){
    double tmp[9];
    int i;
    
    tmp[0] = (matA[0] * matB[0]) + (matA[1] * matB[3]) + (matA[2] * matB[6]);
    tmp[1] = (matA[0] * matB[1]) + (matA[1] * matB[4]) + (matA[2] * matB[7]);
    tmp[2] = (matA[0] * matB[2]) + (matA[1] * matB[5]) + (matA[2] * matB[8]);
    tmp[3] = (matA[3] * matB[0]) + (matA[4] * matB[3]) + (matA[5] * matB[6]);
    tmp[4] = (matA[3] * matB[1]) + (matA[4] * matB[4]) + (matA[5] * matB[7]);
    tmp[5] = (matA[3] * matB[2]) + (matA[4] * matB[5]) + (matA[5] * matB[8]);
    tmp[6] = (matA[6] * matB[0]) + (matA[7] * matB[3]) + (matA[8] * matB[6]);
    tmp[7] = (matA[6] * matB[1]) + (matA[7] * matB[4]) + (matA[8] * matB[7]);
    tmp[8] = (matA[6] * matB[2]) + (matA[7] * matB[5]) + (matA[8] * matB[8]);
    
    // Wee need this tmp-var to assure that we can use matA OR matB == matR
    for(i = 0; i < 9; i++){
        matR[i] = tmp[i];
    }
}


/**
 * @brief  matT = matA^T
 * 
 * matrix size = 3 x 3
 * 
 */
void GWHL_matTranspose3d(double *matA, double *matT){
    double tmp[9];
    int i;
    
    tmp[0] = matA[0];
    tmp[1] = matA[3];
    tmp[2] = matA[6];
    tmp[3] = matA[1];
    tmp[4] = matA[4];
    tmp[5] = matA[7];
    tmp[6] = matA[2];
    tmp[7] = matA[5];
    tmp[8] = matA[8];
    
    // Wee need this tmp-var to assure that we can use matA == matT
    for(i=0; i<9; i++){
        matT[i] = tmp[i];
    }
}


void  GWHL_matTransposeNd(double *matA, long long int m, long long int n, double *matT){
    long long int i, j;
    long long int r = 0;
    for(i = 0; i < n; i++){ // for every coloumn
	for(j = 0; j < m; j++){ // for every row/entry in coloumn
	    matT[r] = matA[j*n+i];
	    r++;
	}
    }
}

/**
 * @brief mat = v . v^T
 */
void GWHL_vecVecTMul(double *vecAIn, double *vecBIn, double *matOut, long long unsigned int size){
  long long unsigned int i,j;
  for(i = 0; i < size; i++){
    for(j = 0; j < size; j++){
      matOut[i*size+j] = vecAIn[i] * vecBIn[j];
    }
  }
}


/**
 * @brief mat += v . v^T
 */
void GWHL_vecVecTMulUpdate(double *vecAIn, double *vecBIn, double *matInOut, long long unsigned int size){
  long long unsigned int i,j;
  for(i = 0; i < size; i++){
    for(j = 0; j < size; j++){
      matInOut[i*size+j] += vecAIn[i] * vecBIn[j];
    }
  }
}


/**
 * @brief mat += (v . v^T) * weight
 */
void GWHL_vecVecTMulScaleUpdate(double *vecAIn, double *vecBIn, double weight, double *matInOut, long long unsigned int size){
  long long unsigned int i,j;
  for(i = 0; i < size; i++){
    for(j = 0; j < size; j++){
      matInOut[i*size+j] += (vecAIn[i] * vecBIn[j] * weight);
    }
  }
}


/**
 * @brief  vecR = vecA + vecB
 * 
 * vector size = 3
 * 
 */
void GWHL_vecVecAdd3d(double *vecA, double *vecB, double *vecR){
    vecR[0] = vecA[0] + vecB[0];
    vecR[1] = vecA[1] + vecB[1];
    vecR[2] = vecA[2] + vecB[2];
}


/**
 * @brief  vecR = vecA * a
 * 
 * vector size = 3
 * 
 */
void GWHL_vecScale3d(double *vecA, double a, double *vecR){
    vecR[0] = vecA[0] * a;
    vecR[1] = vecA[1] * a;
    vecR[2] = vecA[2] * a;
}


/**
 * @brief  vecR += vecA * a
 * 
 * vector size = length
 * this routine is not aliasing safe
 * 
 */
void GWHL_vecScaleNd(double *vecA, double a, unsigned long long int length, double *vecR){
  unsigned long long int i;
  for(i = 0; i < length; i++){
    vecR[i] = vecA[i] * a;
  }
}


/**
 * @brief  vecR += vecA * a
 * 
 * vector size = length
 * this routine is not aliasing safe
 */
void GWHL_vecScaleNdUpdate(double *vecA, double a, unsigned long long int length, double *vecRupdated){
  unsigned long long int i;
  for(i = 0; i < length; i++){
    vecRupdated[i] += vecA[i] * a;
  }
}


/**
* @brief normalize a vector of length "vl" in place
* 
* @param vec vector to normalize (in/out)
* @param vl vector length
*/
void GWHL_normalizeVectorN(double *vec, unsigned long long int vl){
    double s = 0;
    unsigned long long int  i = 0;
    
    for(i = 0; i < vl; i++){
        s += vec[i]*vec[i];
    }
    
    s = 1 / sqrt(s);
    
    for(i = 0; i < vl; i++){
        vec[i] = vec[i] * s;
    }
}

/**
* @brief normalize a vector of length "vl" in place using the max norm v / |max_i(v)|
* 
* @param vec vector to normalize (in/out)
* @param vl vector length
*/
void GWHL_normalizeVectorNMaxNorm(double *vec, unsigned long long int vl){
    double s = 0;
    unsigned long long int  i = 0;
    s = 1 / fabs(GWHL_Helpers_getMaxAbsFromData(vec, vl));
    
    for(i = 0; i < vl; i++){
        vec[i] = vec[i] * s;
    }
}

/**
* @brief normalize a 3d vector in place
* 
* @param vec vector to normalize (in/out)
*/
void GWHL_normalizeVector3d(double *vec){
    double s = 1 / sqrt(vec[0]*vec[0] + vec[1]*vec[1] +vec[2]*vec[2]);
    
    vec[0] = vec[0] * s;
    vec[1] = vec[1] * s;
    vec[2] = vec[2] * s;
}

/**
 * resets an array or matrix to all 0.0 values using efficient memset. array must be initialized
*/
void GWHL_zeroArray(double *array, uint64_t size){
    memset(array, 0, size*sizeof(double));
}

/**
 * adds two arbitrary long arrays, result will be in "b"
 */
void GWHL_arrayAddN(double *aIn, double *bInOut, uint64_t size){
    uint64_t i;
    for(i = 0; i < size; i++){
        bInOut[i] += aIn[i];
    }
}

/**
 * adds two arbitrary long arrays, result will be in "b"
 */
void GWHL_vecVecAddN(double *aIn, double *bInOut, uint64_t size){
    uint64_t i;
    for(i = 0; i < size; i++){
        bInOut[i] += aIn[i];
    }
}

/**
 * add two arbitrary arrays, results will be in c and overwritten
 */
void GWHL_vecVecAddNCopy(double *aIn, double *bIn, double *cOut, uint64_t size){
    uint64_t i;
    for(i = 0; i < size; i++){
        cOut[i] = aIn[i] + bIn[i];
    }
}

void GWHL_vecVecAddScaleN(double *aIn, double *bInOut, double s, uint64_t size){
    uint64_t i;
    for(i = 0; i < size; i++){
        bInOut[i] = aIn[i] + s * bInOut[i];
    }
}

/**
 * a = a * b
 */
void GWHL_vecVecMulN(double *aInOut, double *bIn, uint64_t len){
    uint64_t i;
    for(i = 0; i < len; i++){
        aInOut[i] = aInOut[i] * bIn[i];
    }
}

/**
 * c = a * b
 */
void   GWHL_vecVecMulN_copy(double *a, double *b, double *cOut, uint64_t len){
    uint64_t i;
    for(i = 0; i < len; i++){
        cOut[i] = a[i] * b[i];
    }
}

/**
 * a = a / b
 */
void GWHL_vecVecDivN(double *a, double *b, uint64_t len){
    uint64_t i;
    for(i = 0; i < len; i++){
        a[i] = a[i] / b[i];
    }
}

/**
 * a[i] = a[i] / b
 */
void GWHL_vecDivN(double *a, double b, uint64_t len){
    uint64_t i;
    for(i = 0; i < len; i++){
        a[i] = a[i] / b;
    }
}


/**
* @brief standard right handed cross product for 3d vectors
* 
* @param vecA input vector a
* @param vecB input vector b
* @param cross output a x b
*/
void GWHL_crossProductVec3d(double *vecA, double *vecB, double *cross){
    double out[3] = {0.0, 0.0, 0.0};
    
    out[0] = vecA[1] * vecB[2] - vecA[2] * vecB[1];
    out[1] = vecA[2] * vecB[0] - vecA[0] * vecB[2];
    out[2] = vecA[0] * vecB[1] - vecA[1] * vecB[0];
    
    cross[0] = out[0];
    cross[1] = out[1];
    cross[2] = out[2];
}



/**
* @brief computes the euclidean norm (length) of the input vector
* 
* @param vIn p_vIn : input vector
* @param len p_len : length of input vector
* @return double - euclidean norm of vIn
*/
double GWHL_vecNormNd(double *vIn, long long int len){
    double sqrtSum = GWHL_vecVecDotNd(vIn, vIn, len);
    return sqrt(sqrtSum);
}


/**
* @brief computes the sqrt of each element in vInOut
* 
* @param vInOut p_vInOut:input and output array, is overwritten
* @param len p_len:length of array
*/
void GWHL_vecSqrtElemWise(double *vInOut, long long unsigned int len){
    long long unsigned int i;
    for(i=0; i<len; i++){
        vInOut[i] = sqrt(vInOut[i]);
    }
}


/**
* @brief returns the row sum norm of a MxN matrix
* 
* @param A p_A:matrix
* @param M p_M:dimension M
* @param N p_N:dimension N
* @return double - row sum norm
*/
double GWHL_LinAlg_getRowSumNorm(double *A, long long unsigned int M, long long unsigned int N){
    long long unsigned int i, j;
    double sums[M];
    
    GWHL_zeroArray(sums, M);
    
    for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
            sums[i] += fabs(A[N*i + j]);
        }
    }
    
    return GWHL_Helpers_getMaxFromData(sums, M);
}

/**
* @brief returns the coloumn sum norm of a MxN matrix
* 
* @param A p_A:matrix
* @param M p_M:dimension M
* @param N p_N:dimension N
* @return double - coloumn sum norm
*/
double GWHL_LinAlg_getColSumNorm(double *A, long long unsigned int M, long long unsigned int N){
    long long unsigned int i, j;
    double sums[N];
    
    GWHL_zeroArray(sums, N);
    
    for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
            sums[j] += fabs(A[N*i + j]);
        }
    }
    
    return GWHL_Helpers_getMaxFromData(sums, N);
}



void GWHL_LinAlg_mirrorMatrixUpperToLower(double *A, long long unsigned int M, long long unsigned int N){
    long long unsigned int i,j;
    
    for(i = 0; i < M; i++){
        for(j = i; j < N; j++){
            A[N*j + i] = A[N*i + j];
        }
    }
}
