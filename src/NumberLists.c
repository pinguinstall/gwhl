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
 * @file   NumberLists.c
 * @author Robin Geyer (robin.geyer@tu-dresden.de)
 * @date   October, 2015
 * @brief  Functions for generating lists of numbers
 *
 */

#include <stdio.h>
#include <math.h>
#include <inttypes.h>
#include <omp.h>

/**
 * @brief   generates list with exponetial growing spacing with double values
 *
 * @param [in] start - start value
 * @param [in] end - last value
 * @param [in] numPoints - number of entries
 * @param [out] result - resulting vector
 *
 */
void GWHL_getLogscaleList(double start, double end, long long int numPoints, double *result){
    double k;
    unsigned long int i;
    
    k = pow((end/start),(1/(((double) numPoints) - 1)));

    for(i=0; i<numPoints; i++){
        result[i] = start * pow(k,i+1);
    }
    
    if(numPoints == 1){
        result[0] = start;
    }
}


/**
 * @brief   population count
 *
 * @param [in] vec - start value
 * @param [in] n - last value
 * @return number of nonzero entries
 *
 */
long long int GWHL_popcnt(int *vec, long long unsigned int n){
    long long int i;
    long long int cnt = 0;
    for(i=0; i<n; i++){
        if(vec[i] != 0){
            cnt++;
        }
    }
    
    return cnt;
}


/**
 * @brief   generates a linear list with double entries
 *
 * @param [in] start - start value (inclusive)
 * @param [in] end - last value (exclusive)
 * @param [in] numPoints - number of entries
 * @param [out] result - result list
 * 
 */
void GWHL_getLinearList(double start, double end, long long int numPoints, double *result){
    long long int i;
    
    long double inc = ((long double) end - (long double) start) / (long double) numPoints;
    long double val = (long double) start;
    
    if(numPoints == 1){
        result[0] = start;
        return;
    }
    
    for(i=0; i<numPoints; i++){
        result[i] = (double) val;
        val += inc;
    }
}

/**
 * @brief   generates a linear list with double entries
 *
 * @param [in] start - start value
 * @param [in] end - last value
 * @param [in] numPoints - number of entries
 * @param [out] result - result list
 * 
 */
void GWHL_getLinearCenterList(double start, double end, long long int numBuckets, double *result){
    long long int i;
    
    long double inc = ((long double) end - (long double) start) / (long double) numBuckets;
    long double val = ((long double) start) + inc / 2.0;

    if(numBuckets == 1){
        result[0] = (end - start) / 2.0;
        return;
    }
    
    for(i=0; i<numBuckets; i++){
        result[i] = (double) val;
        val += inc;
    }
    

}

/**
 * @brief   generates a linear list with integer entries
 *
 * @param [in] start - start value
 * @param [in] end - last value
 * @param [in] numPoints - number of entries
 * @param [out] result - result list
 * 
 */
void GWHL_getLinearList_int(long long int start, long long int end, long long int numPoints, long long int *result){
    long long int i;
    
    long long int inc = (end - start) / numPoints;
    long long int val = start;
    
    if(inc == 0) inc = 1;

    for(i=0; i<numPoints; i++){
        result[i] = val;
        val += inc;
    }
    
    if(numPoints == 1){
        result[0] = start;
    }
}


/**
 * adds a vector of length n to one of m >= n, with respect to mask which is also of length m and
 * has n elements set to 1.
 * 
 * vec    = (1, 2, 3)
 * result = (3, 3, 3, 3, 3)
 * mask   = (1, 0, 0, 1, 1)
 * 
 * result after call =
 *          (4, 3, 3, 5, 6)
 */
void GWHL_addWithMask(double *vec, long long int vecl, double *result, long long int reslen, int *mask){
    long long int i, m = 0;
    for(i = 0; i < reslen; i++){
        if(mask[i] != 0){
            result[i] += vec[m];
            m++;
        }
        // that should not happen but is a sanity check
        if(m >= vecl){
            return;
        }
    }
}

/**
 * set the entries of a compressed vector "result" using the "mask", from a uncompressed vector "vec"
 */
void GWHL_setCompressedWithMask(double *vec, long long int vecl, double *result, long long int reslen, int *mask){
    long long int i, m = 0;
    for(i = 0; i < vecl; i++){
        if(mask[i] != 0){
            result[m] = vec[i];
            m++;
        }
        if(m >= reslen){
            return;
        }
    }
}

/**
 * set the entries of a uncompressed vector "result" using the "mask", from a compressed vector "vec"
 */
void GWHL_setUncompressedWithMask(double *vec, long long int vecl, double *result, long long int reslen, int *mask){
   
    long long int i, m = 0;
    for(i = 0; i<reslen; i++){
        if(mask[i] != 0){
            result[i] = vec[m];
            m++;
        }
        if(m >= vecl){
            return;
        }
    }
}

/**
 * initialize whole array with number x
 */
void GWHL_setArrayToXd(double *array, double x, uint64_t arraylen){
    uint64_t i;
    #pragma omp parallel for schedule(static) if(arraylen > 100000)
    for(i = 0; i < arraylen; i++){
        array[i] = x;        
    }
}

/**
 * initialize whole array with number x
 */
void GWHL_setArrayToXi(int *array, int x, uint64_t arraylen){
    uint64_t i;
    #pragma omp parallel for schedule(static) if(arraylen > 100000)
    for(i = 0; i < arraylen; i++){
        array[i] = x;        
    }
}

void GWHL_setArrayToXllu(long long unsigned int *array, long long unsigned int x, uint64_t arraylen){
    uint64_t i;
    #pragma omp parallel for schedule(static) if(arraylen > 100000)
    for(i = 0; i < arraylen; i++){
        array[i] = x;        
    }
}

void GWHL_setArrayToXll(long long int *array, long long int x, uint64_t arraylen){
    uint64_t i;
    #pragma omp parallel for schedule(static) if(arraylen > 100000)
    for(i = 0; i < arraylen; i++){
        array[i] = x;        
    }
}

/**
 * copy array
 */
void GWHL_copyArray(double *src, double *dest, uint64_t len){
    uint64_t i;
    
    #pragma omp parallel for schedule(static) if(len > 100000)
    for(i=0; i<len; i++){
        dest[i] = src[i];
    }
}

/**
 * @brief set all entries to its absolute value (in-place replacement)
 * 
 * @param vec : input vector
 * @param len : vector length
 */
void GWHL_Lists_setAbsoluteEntries(double *vec, long long unsigned int len){
    long long unsigned int i = 0;
    
    #pragma omp parallel for schedule(static) if(len > 100000)
    for(i = 0; i < len; i++){
        vec[i] = fabs(vec[i]);
    }
}

/**
 * @brief set all entries to its squared value (in-place replacement)
 * 
 * @param vec : input vector
 * @param len : vector length
 */
void GWHL_Lists_setSquaredEntries(double *vec, long long unsigned int len){
    long long unsigned int i = 0;
    
    #pragma omp parallel for schedule(static) if(len > 100000)
    for(i = 0; i < len; i++){
        vec[i] = vec[i] * vec[i];
    }
}
