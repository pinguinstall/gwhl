#ifndef ___HELPERS__H_
#define ___HELPERS__H_

double GWHL_Helpers_getMaxAbsFromData(const double *datain, long long int len);
double GWHL_Helpers_getMinAbsFromData(const double *datain, long long int len);

double GWHL_Helpers_getMaxFromData(const double *datain, long long int len);
double GWHL_Helpers_getMinFromData(const double *datain, long long int len);
double GWHL_Helpers_getMaxFromUpperRT(const double *matIn, long long int N);
double GWHL_Helpers_getMinFromUpperRT(const double *matIn, long long int N);

#endif
