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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <float.h>
#include "../include/GWHelperLib.h"

GWHL_StatsContainer* GWHL_getStandardStatsFromDoubleData(double *datVec, unsigned long long int size, int inPlace){
  
  double *myData;
  double sumxi = 0;
  double sumxi2 = 0;
  double sumxi3 = 0;
  double sumxi4 = 0;
  
  GWHL_StatsContainer *cont = calloc(1, sizeof(GWHL_StatsContainer));
  cont->min = DBL_MAX;
  cont->max = -DBL_MAX;
  
  // sort (and if wished copy) the data
  if(inPlace){
    qsort(datVec, size, sizeof(double), GWHL_compareDouble);
    myData = datVec;
  } else {
    myData = malloc(size * sizeof(double));
    memcpy(myData, datVec, size*sizeof(double));
    qsort(myData, size, sizeof(double), GWHL_compareDouble);
  }
  
  // loop over data in order to get the stat moments
  unsigned long long int k;
  for(k = 0; k < size; k++){
    if(myData[k] > cont->max) cont->max = myData[k];
    if(myData[k] < cont->min) cont->min = myData[k];
    
    // build the sums from which we can calculate all the statistics
    sumxi += myData[k];
    sumxi2 += pow(myData[k],2);
    sumxi3 += pow(myData[k],3);
    sumxi4 += pow(myData[k],4);
  
    // increment the counter for valid values
    (cont->numPoints)++;
  }
  
  // compute the statistical moments
  cont->mean = sumxi / (double) cont->numPoints;
  cont->variance = (sumxi2 - cont->numPoints * pow(cont->mean, 2)) / (cont->numPoints - 1);
  cont->stddev = sqrt(cont->variance);
  
  double tmp = 0.0;
  tmp = sumxi3 - ((3.0 * cont->mean * ( sumxi2 - (cont->mean * sumxi))) + pow(cont->mean, 3) * cont->numPoints );
  tmp *= sqrt(cont->numPoints);
  tmp /= pow(sqrt(sumxi2 - pow(cont->mean, 2) * cont->numPoints), 3);
  cont->skewness = tmp;
  
  double s = sqrt(sumxi2 - pow(cont->mean, 2) * cont->numPoints);
  tmp = sumxi4 - (((4. * cont->mean * sumxi3)
    - (6. * pow(cont->mean,2) * sumxi2) + (4. * pow(cont->mean,3) * sumxi)) - (pow(cont->mean,4) * cont->numPoints));
  tmp *= (cont->numPoints / s / s / s / s);
  cont->kurtosis = tmp;
  
  // compute the quantiles
  // Weibull method, equals mathematica: {{0, 1}, {0, 1}} mean-based estimate (Weibull method)
  double p[] = {0.05, 0.1, 0.5, 0.90, 0.95};
  
  double t = p[0] * ((double) size + 1.0);
  long long unsigned int i = (long long unsigned int) t;
  double di = (double) i;
  if(i == 0){
    cont->q5 = myData[i];
  } else {
    cont->q5 = ((di + 1) - t) * myData[i-1] + (t - di) * myData[i];
  }
 
  t = p[1] * ((double) size + 1.0);
  i = (long long unsigned int) t;
  di = (double) i;
  if(i == 0){
    cont->q10 = myData[i];
  } else {
    cont->q10 = ((di + 1) - t) * myData[i-1] + (t - di) * myData[i];
  }
  
  t = p[2] * ((double) size + 1.0);
  i = (long long unsigned int) t;
  di = (double) i;
  cont->median = ((di + 1) - t) * myData[i-1] + (t - di) * myData[i];
  
  t = p[3] * ((double) size + 1.0);
  i = (long long unsigned int) t;
  di = (double) i;
  cont->q90 = ((di + 1) - t) * myData[i-1] + (t - di) * myData[i];
  
  t = p[4] * ((double) size + 1.0);
  i = (long long unsigned int) t;
  di = (double) i;
  cont->q95 = ((di + 1) - t) * myData[i-1] + (t - di) * myData[i];
  
  if(!inPlace){
    free(myData);
  }
  
  return cont;
}
