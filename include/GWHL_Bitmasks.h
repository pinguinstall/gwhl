#ifndef ___GWHL_BITMASKS__H__
#define ___GWHL_BITMASKS__H__

#include <stdint.h>

typedef struct gwhl_bitmask_long {
    uint64_t numbits;
    uint64_t numlongbins;
    uint64_t *mask;
} gwhl_bitmask_long;

#endif
