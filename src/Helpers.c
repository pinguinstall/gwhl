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
#include <limits.h>
#include <float.h>

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
                if(i < sizeM - 1){
                    fprintf(stream, last, mat[i*sizeN+j]);
                } else if (i == sizeM - 1 && endWithNewline){
                    fprintf(stream, last, mat[i*sizeN+j]);
                } else {
                    fprintf(stream, first, mat[i*sizeN+j]);
                }
            }
        }
    }

    if(endWithNewline > 1){
        for(i = 0; i < endWithNewline - 1; i++){
            fprintf(stream, "\n");
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
* @brief fprintf to two files with one call
* 
* @param f2 p_f1:first file
* @param f1 p_f2:second file
* @param fmtstr p_fmtstr:format string
*/
void GWHL_Helpers_print_2files(FILE *f1, FILE *f2, char const *fmtstr, ...){ 
    va_list valist;
    if(f1 != NULL){
        va_start(valist, fmtstr);
        vfprintf(f1, fmtstr, valist);
        va_end(valist);
    }
    
    if(f2 != NULL){
        va_start(valist, fmtstr);
        vfprintf(f2, fmtstr, valist);
        va_end(valist);
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
                if(fabs(datain[i]) > fabs(max)){
                    max = datain[i];
                }
        }
        return max;
}


/**
 * get the value from datain which has min(abs(datain[i]))
 * @return datain[i] where min(abs(datain[i])) over all i
 */
double GWHL_Helpers_getMinAbsFromData(const double *datain, long long int len){
        long long int i;
        double min = DBL_MAX;
        for(i = 0; i < len; i++){
                if(fabs(datain[i]) < fabs(min)){
                        min = datain[i];
                }
        }
        return min;
}


/**
 * get the value from datain which has max(datain[i])
 * @return datain[i] where max(datain[i]) over all i
 */
double GWHL_Helpers_getMaxFromData(const double *datain, long long int len){
        long long int i;
        double max = -DBL_MAX;
        
        for(i = 0; i < len; i++){
                if(datain[i] > max){
                        max = datain[i];
                }
        }
        return max;
}


/**
 * get the value from datain which has min(datain[i])
 * @return datain[i] where max(datain[i]) over all i
 */
double GWHL_Helpers_getMinFromData(const double *datain, long long int len){
        long long int i;
        double min = DBL_MAX;
        
        for(i = 0; i < len; i++){
                if(datain[i] < min){
                        min = datain[i];
                }
        }
        return min;
}


/**
* @brief returns the value where max(abs(matIn[i])) from upper right triangle of the NxN matrix (excluding the diagonal)
* 
* @param matIn - input matrix
* @param N - matrix size in both dimensions
* @return matIn[i] where max(abs(matIn[i])) in upper RT
*/
double GWHL_Helpers_getMaxAbsFromUpperRT(const double *matIn, long long int N){
        long long int i, j, k = 1;
        double max = 0.0;
        for (i = 0; i < N; i++) {
                for (j = k; j < N; j++) {  
                        if(fabs(matIn[N * i + j]) > fabs(max)){
                                max = matIn[N * i + j];
                        }
                }
                k++;
        }
        return max;
}


/**
* @brief returns maximum from the upper right triangle (excluding the diagonal) of matrix matIn
* 
* @param matIn - input matrix
* @param N - matrix size in both dimensions
* @return maximum upper RT value
*/
double GWHL_Helpers_getMaxFromUpperRT(const double *matIn, long long int N){
        long long int i, j, k = 1;
        double max = -DBL_MAX;
        for (i = 0; i < N; i++) {
                for (j = k; j < N; j++) {  
                        if(matIn[N * i + j] > max){
                                max = matIn[N * i + j];
                        }
                }
                k++;
        }
        return max;
}


/**
* @brief returns minimum from the upper right triangle (excluding the diagonal) of matrix matIn
* 
* @param matIn - input matrix
* @param N - matrix size in both dimensions
* @return minimum upper RT value
*/
double GWHL_Helpers_getMinFromUpperRT(const double *matIn, long long int N){
        long long int i, j, k = 1;
        double min = DBL_MAX;
        for (i = 0; i < N; i++) {
                for (j = k; j < N; j++) {  
                        if(matIn[N * i + j] < min){
                                min = matIn[N * i + j];
                        }
                }
                k++;
        }
        return min;
}
