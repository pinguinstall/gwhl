#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

// will not engage OpenMP if less than this length
#define SEQ_THRESHOLD 10000

//BUG: here the list pointers get swapped and cant be free()'d : "munmap_chunk(): invalid pointer
void swap(double *a, double *b){
    double tmp = 0;
    tmp = *a;
    *a = *b;
    *b = tmp;
}

long long int partition(double *d, long long int p, long long int ri){
    double pivot = d[p];
    long long int left = p;
    long long int right = ri;
    
    
    while (left < right){
        // increase until you get to a point where the element is greater that the pivot
        while (left < ri && d[left] <= pivot) ++left;
        // increase until you get to a point where the element is less or equal to the pivot
        while (right >= 0 && d[right] > pivot) --right;
        // swap if within bounds
        if (left < right && left < ri && right >= 0){
            swap(&d[left], &d[right]);
        }
    }

    // swap at the end
    if (left < right && left <ri && right >= 0){
        swap(&d[left], &d[right]);
    }
    d[p] = d[right];
    d[right] = pivot;
    return right;
}

void p_quicksort(double *data, long long int p, long long int low){
    long long int div;
    
    if(p < low){ 
        div = partition(data, p, low);
        #pragma omp task shared(data) if(low - p > SEQ_THRESHOLD) 
        p_quicksort(data, p, div - 1); 
        #pragma omp task shared(data) if(low - p > SEQ_THRESHOLD) 
        p_quicksort(data, div + 1, low); 
        #pragma omp taskwait
    }
}


/**
 * @brief parallel quicksort implementation
 * 
 * @param data : array of doubles to sort (this will happen in-place)
 * @param len : length of array
 * @todo make workable for arbitrary data types
 */
void GWHL_Sorting_parallel_quicksort(double *data, long long unsigned int len){
    omp_set_dynamic(0);
//     fprintf(stderr, "outer num omp threads = %i\n", omp_get_num_threads());
//     fprintf(stderr, "outer data = %p\n", data);
        #pragma omp parallel shared(data, len)
        {
            #pragma omp master
            {
                p_quicksort(data, 0, len-1);
            #pragma omp taskwait
            }
        }
}
