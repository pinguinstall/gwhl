#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/GWHelperLib.h"


/**
 * @brief computes the quantile of the given data in the vector
 * 
 * @param vec array of doubles
 * @param len length of array
 * @param quantile quantile
 * @param inPlace proces in-place (will sort data in the array), meaningless if isSortedData = 1 
 * @param isSortedData indicated whether data is already sorted (will not be checked internaly)
 * @return the quantile asked for
 */
double GWHL_Stat_quantile(double *vec, long long unsigned int len, double quantile, int inPlace, int isSortedData){
    if(quantile > 1.0 || quantile < 0.0){
        fprintf(stderr, "GWHL_Stat_quantile: quantile must be between 0 and 1");
    }
    
    if(len < 2) {
        return vec[0];
    }

    // sort
    double *tmp;
    if(!isSortedData){
        if(inPlace){
            qsort(vec, len, sizeof(double), GWHL_compareDouble);
            tmp = vec;
        } else {
            tmp = calloc(len, sizeof(double));
            memcpy(tmp, vec, len * sizeof(double));
            qsort(tmp, len, sizeof(double), GWHL_compareDouble);
        }
    } else { // if data is sorted, inplace is meaningless
        tmp = vec;
    }
    
    // if trivial quantile
    if(quantile == 1) {
        return tmp[len-1];
    }

    if(quantile == 0) {
        return tmp[0];
    }

    // Calculate quantile
    double t = quantile * (len - 1);
    long long int i = (long long int) t; // floor
    double p1 = tmp[i];
    double p2 = tmp[i+1];

    return ((((i + 1) - t) * p1) + ((t - i) * p2));
}

/**
 * @brief returns the RSE (robust scatter estimate) of the data in te array
 * 
 * @param vec array with the data
 * @param len length of data
 * @param inPlace proces in-place (will sort data in the array)
 * @return RSE
 */
double GWHL_Stat_RSE(double *vec, long long unsigned int len, int inPlace, int isSortedData){
    double q1 = 0.0;
    double q9 = 0.0;
    double *tmp;
    
    if(!inPlace && !isSortedData){
        tmp = calloc(len, sizeof(double));
        memcpy(tmp, vec, len * sizeof(double));
        qsort(tmp, len, sizeof(double), GWHL_compareDouble);
        
        q1 = GWHL_Stat_quantile(tmp, len, 0.1, 1, 1);    
        q9 = GWHL_Stat_quantile(tmp, len, 0.9, 1, 1);
    } else {
        tmp = vec;
        q1 = GWHL_Stat_quantile(tmp, len, 0.1, inPlace, isSortedData);    
        q9 = GWHL_Stat_quantile(tmp, len, 0.9, inPlace, 1);
    }
    
    return GWHLMATHCONSTANTS_RSE_FACTOR * (q9 - q1);
}


/**
 * @brief returns the RMS of data, weighted if necessary
 * 
 * @param vec data input
 * @param weights equally long array of weights (must be set to NULL if unweighted)
 * @param len length of arrays
 * @return RMS or weighted RMS value
 */
double GWHL_Stat_RMS(double *vec, double *weights, long long unsigned int len) {
    long long unsigned int i = 0;
    double s = 0.0;
    double ws = 0.0;
    if(weights != NULL){
        for(i = 0; i < len; i++){
            s += weights[i] * (vec[i] * vec[i]);
            ws += weights[i];
        }
        
        return sqrt(s/ws);
    } else {
        for(i = 0; i < len; i++){
            s += (vec[i] * vec[i]);
        }
        
        return sqrt(s/((double) len));
    }
}


/**
 * @brief returns the weighted average
 * 
 * @param dat array of data
 * @param weights array of weights
 * @param len len of arrays
 * @return weighted mean
 */
double GWHL_Stat_wMean(double *dat, double *weights, long long unsigned len) {
    long long unsigned int i = 0;
    double s = 0.0;
    double sw = 0.0;
    for(i = 0; i < len; i++){
        s += weights[i] * dat[i];
        sw += weights[i];
    }
    
    return s / sw;
}


/**
 * @brief computes the weighted or un-weighted R2 value (Coefficient of determination)
 * 
 * @param dat data (observed)
 * @param mod prediction (values for data-points given by model)
 * @param weights weights (un-weighted if this is NULL)
 * @param len length of arrays
 * @return R2
 */
double GWHL_Stat_R2(double *dat, double *mod, double *weights, long long unsigned int len){
    double mean = 0.0;
    double sr = 0.0;
    double sv = 0.0;
    long long unsigned int i = 0;
    if(weights != NULL){
        mean = GWHL_Stat_wMean(dat, weights, len);
        for(i = 0; i < len; i++){
            sr += weights[i] * pow(dat[i] - mod[i], 2.0);
            sv += weights[i] * pow(dat[i] - mean ,2.0);
        }
    } else {
        // simple mean
        for(i = 0; i < len; i++){
            mean += dat[i];
        }
        mean = mean / (double) len;
        
        for(i = 0; i < len; i++){
            sr += pow(dat[i] - mod[i], 2.0);
            sv += pow(dat[i] - mean, 2.0);
        }
    }
    
    return 1.0 - (sr/sv);
}

/**
 * @brief computes the weighted or un-weighted R2 value (Coefficient of determination)
 * 
 * @param dat data (observed)
 * @param mod prediction (values for data-points given by model)
 * @param weights weights (un-weighted if this is NULL)
 * @param len length of arrays
 * @return R2
 */


/**
 * @brief computes the weighted or un-weighted R2 value (Coefficient of determination) from residuals
 * 
 * @param dat data (observed)
 * @param residuals residuals from the fit
 * @param weights weights of observations
 * @param len length of arrays
 * @return double
 */
double GWHL_Stat_R2_Residuals(double *dat, double *residuals, double *weights, long long unsigned int len){
    double mean = 0.0;
    double sr = 0.0;
    double sv = 0.0;
    long long unsigned int i = 0;
    if(weights != NULL){
        mean = GWHL_Stat_wMean(dat, weights, len);
        for(i = 0; i < len; i++){
            sr += weights[i] * (residuals[i] * residuals[i]);
            sv += weights[i] * pow(dat[i] - mean, 2.0);
        }
    } else {
        // simple mean
        for(i = 0; i < len; i++){
            mean += dat[i];
        }
        mean = mean / (double) len;
        
        for(i = 0; i < len; i++){
            sr += (residuals[i] * residuals[i]);
            sv += pow(dat[i] - mean, 2.0);
        }
    }
    
    return 1.0 - (sr/sv);
}
