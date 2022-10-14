#ifndef ____GWHL_SPARSETYPES__
#define ____GWHL_SPARSETYPES__

typedef struct GWHL_SparseVector_double {
    long long unsigned int len;
    long long unsigned int numEntries;
    long long unsigned int *entryIdxs;
    double *entries;
} GWHL_SparseVector_double;

#endif
