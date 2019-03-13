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

#ifndef ____NUMBERLISTS__H
#define ____NUMBERLISTS__H

#include <inttypes.h>

void GWHL_getLogscaleList(double start, double end, long long int numPoints, double *result);
long long int GWHL_popcnt(int *vec, long long unsigned int n);
void GWHL_getLinearList(double start, double end, long long int numPoints, double *result);
void GWHL_getLinearCenterList(double start, double end, long long int numBuckets, double *result);
void GWHL_getLinearList_int(long long int start, long long int end, long long int numPoints, long long int *result);
void GWHL_addWithMask(double *vec, int vecl, double *result, int reslen, int *mask);
void GWHL_setCompressedWithMask(double *vec, int vecl, double *result, int reslen, int *mask);
void GWHL_setUncompressedWithMask(double *vec, int vecl, double *result, int reslen, int *mask);
void GWHL_setArrayToXd(double *array, double x, uint64_t arraylen);
void GWHL_setArrayToXi(int *array, int x, uint64_t arraylen);
void GWHL_copyArray(double *src, double *dest, uint64_t len);

#endif
