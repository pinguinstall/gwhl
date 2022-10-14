#include <stdlib.h>
#include <math.h>
#include "../include/GWHelperLib.h"

//TODO
double GWHL_Stat_quantile(double *vec, long long unsigned int len, double quantile){
    if(quantile > 1.0 || quantile < 0.0){
        fprintf(stderr, "GWHL_Stat_quantile: quantile must be between 0 and 1");
    }
    
    return 0;
}

//TODO
double GWHL_Stat_RSE(double *vec, long long unsigned int len){
    return 0;
}
