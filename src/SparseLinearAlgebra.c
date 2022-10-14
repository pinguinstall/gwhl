#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include <omp.h>

#include "../include/GWHL_SparseTypes.h"

/**
 * @brief mat += (a . b^T) * w in sparse format, omp locked
 */
void GWHL_vecVecTMulScaleUpdate_sparse_ompsafe(GWHL_SparseVector_double *a, GWHL_SparseVector_double *b, double w, double *MOut, omp_lock_t *locks){
   
    assert(a->len == b->len);
    
    long long unsigned int i, j, n, m;
    
    n = a->len; // matrix dimensions for output matrix
    m = b->len;
    
    for(i = 0; i < a->numEntries; i++){
        omp_set_lock(&(locks[a->entryIdxs[i]])); // lock line in output matrix
        for(j = 0; j < b->numEntries; j++){
            MOut[a->entryIdxs[i] * m + b->entryIdxs[j]] += a->entries[i] * b->entries[j] * w;
        }
        omp_unset_lock(&(locks[a->entryIdxs[i]]));
    }
}


/**
 * @brief mat += (a . b^T) * w in sparse format, not thread safe
 */
void GWHL_vecVecTMulScaleUpdate_sparse(GWHL_SparseVector_double *a, GWHL_SparseVector_double *b, double w, double *MOut){
    
    assert(a->len == b->len);
    
    long long unsigned int i, j, n, m;
    
    n = a->len; // matrix dimensions for output matrix
    m = b->len;
    
    for(i = 0; i < a->numEntries; i++){
        for(j = 0; j < b->numEntries; j++){
            MOut[a->entryIdxs[i] * m + b->entryIdxs[j]] += a->entries[i] * b->entries[j] * w;
        }
    }
}


void GWHL_vecScaleNdUpdate_sparse(GWHL_SparseVector_double *a, double b, long long unsigned int len, double *Vout){
    long long unsigned int i;
    
    for(i = 0; i < a->numEntries; i++){
            assert(a->entryIdxs[i] < len);
            Vout[a->entryIdxs[i]] += b * a->entries[i];
    }
}
