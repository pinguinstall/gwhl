/* GWHelperLib - A helper Library for GaiaGW (and general tasks)
*  Copyright (C) 2015-2022  Robin Geyer, et al.
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include "../include/GWHL_Bitmasks.h"
#include <stdint.h>

//TODO: check if the long masks are not better implemented using void* (consecutive memory arrays)

void GWHL_BitMasks_setAllZero(uint64_t *bitmask){
    *bitmask = 0;
}

void GWHL_BitMasks_setAllOne(uint64_t *bitmask){
    *bitmask = -1; // uint64_t set to -1 is 0xFFFFFFFF
}

void GWHL_BitMasks_setNthToZero(uint64_t *bitmask, int nth){
    *bitmask &= ~(1ULL << nth);
}

void GWHL_BitMasks_setNthToOne(uint64_t *bitmask, int nth){
    *bitmask |= 1ULL << nth; 
}

void GWHL_BitMasks_flipNth(uint64_t *bitmask, int nth){
    *bitmask ^= 1ULL << nth;
}

int GWHL_BitMasks_getNth(uint64_t *bitmask, int nth){
    return (*bitmask >> nth) & 1ULL;

}

int GWHL_BitMasks_checkAllZero(uint64_t *bitmask){
    if(*bitmask == 0){
        return 1;
    }
    
    return 0;
}

int GWHL_BitMasks_checkAllOne(uint64_t *bitmask){
    if(*bitmask == -1){
        return 1;
    }
    
    return 0;
}


/**
 * @brief compare function (checks if both numbers are equal)
 * 
 * @param bm1 bitmask 1
 * @param bm2 bitmask 2
 * @return 1 if equal, 0 if not equal
 */
int GWHL_BitMasks_compare(uint64_t *bm1, uint64_t *bm2){
    if(bm1 == bm2) return 1;
    return 0;
}


// implementation from wikipedia
int popcount64c(uint64_t x){
    //types and constants used in the functions below
    //uint64_t is an unsigned 64-bit integer variable type (defined in C99 version of C language)
    const uint64_t m1  = 0x5555555555555555; //binary: 0101...
    const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
    const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
    const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...
    
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
    return (x * h01) >> 56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
}


int GWHL_BitMasks_popcount(uint64_t *bitmask){
    return popcount64c((uint64_t) *bitmask);
}

int GWHL_BitMasks_large_map_init(uint64_t numBits, gwhl_bitmask_long *bm){
    uint64_t numLongs = numBits / (8 * sizeof(uint64_t));
    
    if(numLongs * (8 * sizeof(uint64_t)) < numBits){
        numLongs++;
    }
    
    bm->numbits = numBits;
    bm->numlongbins = numLongs;
    bm->mask = NULL;
    bm->mask = calloc(numLongs, sizeof(uint64_t));
    if(bm->mask == NULL) return -1;
    return 0;
}

void GWHL_BitMasks_large_map_setAllZero(gwhl_bitmask_long *bm){
    uint64_t i;
    for(i=0; i<bm->numlongbins; i++){
        bm->mask[i] = 0;
    }
}

void GWHL_BitMasks_large_map_setAllOne(gwhl_bitmask_long *bm){
    uint64_t i;
    for(i=0; i<bm->numlongbins; i++){
        bm->mask[i] = -1;
    }
}

int check_in_map(gwhl_bitmask_long *bm, uint64_t nth){
    if(bm->numbits > nth) return 1;
    fprintf(stderr, "try to select bit %lu, map only has %lu (counting starts at 0)\n",
            nth, bm->numbits);
    exit(EXIT_FAILURE);
}

void GWHL_BitMasks_large_map_setNthToZero(gwhl_bitmask_long *bm, uint64_t nth){
    check_in_map(bm, nth);
    uint64_t bin = nth / (8 * sizeof(uint64_t));
    GWHL_BitMasks_setNthToZero(&(bm->mask[bin]), nth % (8 * sizeof(uint64_t)));
}

void GWHL_BitMasks_large_map_setNthToOne(gwhl_bitmask_long *bm, uint64_t nth){
    check_in_map(bm, nth);
    uint64_t bin = nth / (8 * sizeof(uint64_t));
    GWHL_BitMasks_setNthToOne(&(bm->mask[bin]), nth % (8 * sizeof(uint64_t)));
}

void GWHL_BitMasks_large_map_flipNth(gwhl_bitmask_long *bm, uint64_t nth){
    check_in_map(bm, nth);
    uint64_t bin = nth / (8 * sizeof(uint64_t));
    GWHL_BitMasks_flipNth(&(bm->mask[bin]), nth % (8 * sizeof(uint64_t)));
}

int GWHL_BitMasks_large_map_getNth(gwhl_bitmask_long *bm, uint64_t nth){
    check_in_map(bm, nth);
    uint64_t bin = nth / (8 * sizeof(uint64_t));
    return GWHL_BitMasks_getNth(&(bm->mask[bin]), nth % (8 * sizeof(uint64_t)));
}

int GWHL_BitMasks_large_map_compare(gwhl_bitmask_long *bm1, gwhl_bitmask_long *bm2){
    uint64_t i = 0;
    if(bm1->numbits != bm2->numbits || bm1->numlongbins != bm2->numlongbins) return 0;
    for(i = 0; i < bm1->numlongbins; i++){
        if(bm1->mask[i] != bm2->mask[i]) return 0;
    }
    return 1;
}

uint64_t GWHL_BitMasks_large_map_popcnt(gwhl_bitmask_long *bm1){
    uint64_t i = 0;
    uint64_t cnt = 0;
    
    for(i = 0; i < bm1->numlongbins; i++){
        cnt += popcount64c(bm1->mask[i]);
    }
    return cnt;
}

void GWHL_BitMasks_large_map_print(gwhl_bitmask_long *bm, FILE *fp, int newline){
    long long int i = 0;
    for(i = bm->numbits - 1; i >= 0; i--){
        if(i % 4 == 0 && i != bm->numbits){
            fprintf(fp, "%i ", GWHL_BitMasks_large_map_getNth(bm, i));
        } else {
            fprintf(fp, "%i", GWHL_BitMasks_large_map_getNth(bm, i));
        }
    }
    
    for(i = newline; i > 0; i--){
        fprintf(fp, "\n");
    }
}
