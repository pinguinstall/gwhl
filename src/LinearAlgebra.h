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

#ifndef ____GWHL_LINEARALGEBRA__
#define ____GWHL_LINEARALGEBRA__

void GWHL_matScale3d(double *matA, double a, double *matR);
void GWHL_matMatAdd3d(double *matA, double *matB, double *matR);
void GWHL_matVecMul3d(double *matA, double *vecB, double *vecR);
double GWHL_vecVecDot3d(double *vecA, double *vecB);
void GWHL_matMatMul3d(double *matA, double *matB, double *matR);
void GWHL_matTranspose(double *matA, double *matT);
void GWHL_vecVecTMul(double *vecAIn, double *vecBIn, double *matOut, long long unsigned int size);
void GWHL_vecVecTMulUpdate(double *vecAIn, double *vecBIn, double *matInOut, long long unsigned int size);
void GWHL_vecVecAdd3d(double *vecA, double *vecB, double *vecR);
void GWHL_vecScale3d(double *vecA, double a, double *vecR);
void GWHL_vecScaleNd(double *vecA, double a, unsigned long long int length, double *vecR);
void GWHL_vecScaleNdUpdate(double *vecA, double a, unsigned long long int length, double *vecRupdated);
void GWHL_normalizeVectorN(double *vec, int vl);
void GWHL_zeroArray(double *array, uint64_t size);
void GWHL_arrayAddN(double *a, double *b, uint64_t size);
void GWHL_vecVecMulN(double *a, double *b, uint64_t len);
void GWHL_vecVecDivN(double *a, double *b, uint64_t len);
void GWHL_vecDivN(double *a, double b, uint64_t len);
void GWHL_crossProductVec3d(double *vecA, double *vecB, double *cross);
void GWHL_vecVecTMulScaleUpdate(double *vecAIn, double *vecBIn, double weight, double *matInOut, long long unsigned int size);
double GWHL_vecNormNd(double *vIn, long long int len);


#endif
