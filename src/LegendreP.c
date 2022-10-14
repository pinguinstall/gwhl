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
#include <math.h>
#include "../include/GWHL_Constants.h"

inline long long int GWHL_LegendreP_lm2idx(long long int l, long long int m, long long int N){
    return (l*(1 + l))/2 + m;
}


inline long long int GWHL_LegendreP_idx2m(long long int idx, long long int N){
    long long int l = ((-1 + sqrt(1 + 8*idx))/2);
    return idx - (l*(1 + l))/2;
}

inline long long int GWHL_LegendreP_idx2l(long long int idx, long long int N){
    return (-1 + sqrt(1 + 8*idx))/2;
}

void GWHL_LegendreP_idx2lm(long long int idx, long long int N, long long int *l, long long int *m){
    *l = (-1 + sqrt(1 + 8*idx))/2;
    *m = idx - ((*l)*(1 + (*l)))/2;
}

inline long long int GWHL_LegendreP_getOutputArraySize(long long int maxl){
    return (maxl + 1) * (maxl + 2) / 2;
}


/**
 * matrix layout is l first m second
 * (0,0 1,0 1,1 2,0 2,1 2,2 ...)
 */
// void GWHL_LegendreP_getLegendreP(double x, long long int lmax, double *P){
//     double sx = sqrt(1 - x * x);
//     P[0] = 1;
//     long long int l, m;
//     // diagonal elements
//     for(m = 0; m < lmax; m++){
//       P[GWHL_LegendreP_lm2idx(m+1, m+1, lmax)] = (2 * m + 1) * sx * P[GWHL_LegendreP_lm2idx(m, m, lmax)];
//       P[GWHL_LegendreP_lm2idx(m+1, m, lmax)] = (2 * m + 1) * x * P[GWHL_LegendreP_lm2idx(m, m, lmax)];
//     }
//     
//     for(m = 0; m < lmax - 1; m++){
//       for(l = m + 2; l <= lmax; l++){
//         P[GWHL_LegendreP_lm2idx(l, m, lmax)] = ((2*l-1) * x * P[GWHL_LegendreP_lm2idx(l-1, m, lmax)] - (l-1+m) * P[GWHL_LegendreP_lm2idx(l-2, m, lmax)])/(l-m);
//       }  
//     }
// }


void GWHL_LegendreP_getAlphas(long long int lmax, double *alphas){
    double fpi = 4.0 * GWHLMATHCONSTANTS_PI;
    long long int l, m;
    for(l = 0; l <= lmax; l++){
      alphas[GWHL_LegendreP_lm2idx(l, 0, lmax)] = sqrt((2*l+1)/fpi);
      for(m = 1; m <= lmax; m++){
        alphas[GWHL_LegendreP_lm2idx(l, m, lmax)] = -alphas[GWHL_LegendreP_lm2idx(l, m-1, lmax)]/sqrt((l+m)*(l+1-m));
      }
    }
}


void GWHL_LegendreP_getCosSinMap(double lambda, long long int lmax, double *cosmap, double *sinmap){
    cosmap[0] = 1.0;
    sinmap[0] = 0.0;
    double c1, s1;
    c1 = cos(lambda);
    s1 = sin(lambda);
    sinmap[1] = s1;
    cosmap[1] = c1;
    long long int m;
    
    for(m = 2; m <= lmax; m++){
      cosmap[m] = c1 * cosmap[m-1] - s1 * sinmap[m-1];
      sinmap[m] = s1 * cosmap[m-1] + c1 * sinmap[m-1];
    }
}


void GWHL_LegendreP_getBetas(long long int lmax, double *betas){
    double fpi = 4.0 * GWHLMATHCONSTANTS_PI;
    long long int l, m;
    betas[0] = 0; //TODO check if formaly correct
    for(l = 1; l <= lmax; l++){
      betas[GWHL_LegendreP_lm2idx(l, 0, lmax)] = sqrt((2*l+1)/fpi/(l*(l+1)));
      for(m = 1; m <= l; m++){
         betas[GWHL_LegendreP_lm2idx(l, m, lmax)] = -betas[GWHL_LegendreP_lm2idx(l, m-1, lmax)]/sqrt((l+m)*(l+1-m));
      }
    }
}


void GWHL_LegendreP_getAB(double x, long long int lmax, double *A, double *B){
    double sx = sqrt(1 - x * x);
    long long int l, m;
    
    A[0] = 0.0;
    B[0] = 0.0;

    for(l = 1; l <= lmax; l++){
	B[GWHL_LegendreP_lm2idx(l, 0, lmax)] = 0.0;
    }

    B[GWHL_LegendreP_lm2idx(1, 1, lmax)] = 1.0;
    for(m = 1; m < lmax; m++){
	B[GWHL_LegendreP_lm2idx(m+1, m+1, lmax)] = ((2. * m + 1.) * (m + 1.) / m) * sx * B[GWHL_LegendreP_lm2idx(m, m, lmax)];
	B[GWHL_LegendreP_lm2idx(m+1, m, lmax)] = (2. * m + 1.) * x * B[GWHL_LegendreP_lm2idx(m, m, lmax)];
    }
    
    for(m = 1; m < lmax-1; m++){
	for(l = m + 2; l <= lmax; l++){
	  B[GWHL_LegendreP_lm2idx(l, m, lmax)] = ((2.*l-1.)*x*B[GWHL_LegendreP_lm2idx(l-1, m, lmax)]-(l-1+m)*B[GWHL_LegendreP_lm2idx(l-2, m, lmax)])/(l-m);
	}
    }
    
    for(l = 1; l <= lmax; l++){
      A[GWHL_LegendreP_lm2idx(l, 0, lmax)] = sx * B[GWHL_LegendreP_lm2idx(l, 1, lmax)];
    }
    
    
    for(l = 1; l <= lmax; l++){
      for(m = 1; m <= l; m++){
	if(l == m){
	  A[GWHL_LegendreP_lm2idx(l, m, lmax)] = -l*x*B[GWHL_LegendreP_lm2idx(l, m, lmax)]/m;
	} else {
	  A[GWHL_LegendreP_lm2idx(l, m, lmax)] = (-l*x*B[GWHL_LegendreP_lm2idx(l, m, lmax)]+(l+m)*B[GWHL_LegendreP_lm2idx(l-1, m, lmax)])/m;
	}
      }
    }
}

