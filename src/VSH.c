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
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "NumberLists.h"
#include "../include/GWHelperLib.h"
#include <gsl/gsl_matrix.h>
#include <float.h>

unsigned long long int GWHL_VSH_getNumCoefficients(unsigned long long int maxl){
    return 2 * (maxl * (maxl + 2));
}


void GWHL_VSH_getGlmFlm(const int maxl, const double a, const double d, double *glmRI, double *flmRI){
    double *A = calloc(GWHL_LegendreP_getOutputArraySize(maxl), sizeof(double));
    double *B = calloc(GWHL_LegendreP_getOutputArraySize(maxl), sizeof(double));
    double *betas = calloc(GWHL_LegendreP_getOutputArraySize(maxl), sizeof(double));
    
    double *sinmap = calloc((maxl + 1), sizeof(double));
    double *cosmap = calloc((maxl + 1), sizeof(double));
    long long int l, m, i;
    i = 0;
    
    GWHL_LegendreP_getBetas(maxl, betas);
    GWHL_LegendreP_getAB(sin(d), maxl, A, B);
    
    GWHL_LegendreP_getCosSinMap(a, maxl, cosmap, sinmap);
    
    for(l = 1; l <= maxl; l++){
        flmRI[i] = betas[GWHL_LegendreP_lm2idx(l, 0, maxl)] * A[GWHL_LegendreP_lm2idx(l, 0, maxl)] * cosmap[0];
        glmRI[i] = -betas[GWHL_LegendreP_lm2idx(l, 0, maxl)] * B[GWHL_LegendreP_lm2idx(l, 0, maxl)] * sinmap[0];
        i++;
        for(m = 1; m <= l; m++){
                flmRI[i] = 2 * betas[GWHL_LegendreP_lm2idx(l, m, maxl)] * A[GWHL_LegendreP_lm2idx(l, m, maxl)] * cosmap[m];
                glmRI[i] = -2 * betas[GWHL_LegendreP_lm2idx(l, m, maxl)] * B[GWHL_LegendreP_lm2idx(l, m, maxl)] * sinmap[m];
                i++;

                flmRI[i] = -2 * betas[GWHL_LegendreP_lm2idx(l, m, maxl)] * A[GWHL_LegendreP_lm2idx(l, m, maxl)] * sinmap[m];
                glmRI[i] = -2 * betas[GWHL_LegendreP_lm2idx(l, m, maxl)] * B[GWHL_LegendreP_lm2idx(l, m, maxl)] * cosmap[m];
                i++;
        }  
    }
    
    free(sinmap);
    free(cosmap);
    free(A);
    free(B);
    free(betas);
}


void GWHL_VSH_getPowerOfOrdersUnNormalized(const int maxl, const double *vshCoefficients, double *result){
    long long unsigned int i = 0;
    long long int l, m;
    
    for(l = 1; l <= maxl; l++){
        result[l-1] = 0.0;
        result[l-1] += pow(vshCoefficients[i], 2.0);
        i++;
        for(m = 1; m <= l; m++){
            // real
            result[l-1] += 2 * pow(vshCoefficients[i], 2.0);
            i++;
            //imaginary
            result[l-1] += 2 * pow(vshCoefficients[i], 2.0);
            i++;
        }
    }
//     fprintf(stderr, "last i %llu\n", i);
}



void GWHL_VSH_getPowerOfOrdersNormalizedWTilde(const int maxl, const double *vshCoefficients, const double *vshCoefficientsSigmas, double *result){
    long long unsigned int i = 0;
    long long int l, m;
    
    for(l = 1; l <= maxl; l++){
        result[l-1] = 0.0;
        for(m = 0; m <= l; m++){
            if(m == 0) {
                result[l-1] += pow(vshCoefficients[i] / vshCoefficientsSigmas[i], 2.0);
                i++;
            } else {
                // real
                result[l-1] += pow(vshCoefficients[i] / vshCoefficientsSigmas[i], 2.0);
                i++;
                //imaginary
                result[l-1] += pow(vshCoefficients[i] / vshCoefficientsSigmas[i], 2.0);
                i++;
            }
        }
    }
}


void GWHL_VSH_getPowerOfOrdersNormalizedW(const int maxl, const double *vshCoefficients, const double *vshCoefficientsSigmas, double *result){
    long long unsigned int i = 0;
    long long int l, m;
    long long int nul;
    
    for(l = 1; l <= maxl; l++){
        result[l-1] = 0.0;
        for(m = 0; m <= l; m++){
            if(m == 0) {
                nul = i; // we only need sigma_l0
                result[l-1] += pow(vshCoefficients[i] / vshCoefficientsSigmas[nul], 2.0);
                i++;
            } else {
                // real
                result[l-1] += pow(vshCoefficients[i],2.0) / (pow(vshCoefficientsSigmas[nul], 2.0) / 2.0);
                i++;
                //imaginary
                result[l-1] += pow(vshCoefficients[i],2.0) / (pow(vshCoefficientsSigmas[nul], 2.0) / 2.0);
                i++;
            }
        }
    }
}


void GWHL_VSH_getPowerOfOrdersNormalizedWOver(const int maxl, const double *vshCoefficients, const double *vshCoefficientsSigmas, double *result){
    long long unsigned int i = 0;
    long long int l, m;
    double avgl[maxl];
    
    GWHL_setArrayToXd(avgl, 0.0, maxl);
    GWHL_setArrayToXd(result, 0.0, maxl);
    GWHL_VSH_getPowerOfOrdersUnNormalized(maxl, vshCoefficients, result);
    
    
    for(l = 1; l <= maxl; l++){
      for(m = 0; m <= l; m++){
	if(m == 0) {
	    avgl[l-1] += pow(vshCoefficientsSigmas[i], 2.0);
	    i++;
	} else {
	    // real
	    avgl[l-1] += 2.0 * pow(vshCoefficientsSigmas[i], 2.0);
	    i++;
	    //imaginary
	    avgl[l-1] += 2.0 * pow(vshCoefficientsSigmas[i], 2.0);
	    i++;
	}
      }
    }
    
    for(l = 1; l <= maxl; l++){
      result[l-1] = result[l-1] / (avgl[l-1] / (2*l + 1));
    }
    
}



void GWHL_VSH_printVSHCoefficients(FILE *fp, double *coeffs, int maxl){
  long long int numc = GWHL_VSH_getNumCoefficients(maxl);
  long long int i;
  
  for(i = 0; i < numc; i++){
    fprintf(fp, "%+9.5e ", coeffs[i]);
    if(i == numc/2 - 1) fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
}

double GWHL_VSH_getZ_scaled(double Wl, double l, double N){
  double n = N*(2 * l + 1);
  double Z = sqrt(4.5*n)*(pow((Wl/n),1.0/3.0)-(1.0-2.0/(9*n)));
  
  return Z;
}

double GWHL_VSH_getZ(double Wl, double l){  
  return GWHL_VSH_getZ_scaled(Wl, l, 1);
}

/**
 * returns the noise probability as described in Mignard Klioner 2012 Equation (85)
 * we return the probability here as Z is from a normal distribution with zero mean and unit variance
 */
double GWHL_VSH_getZConfidence(double Wl, double l){
  double Z = GWHL_VSH_getZ(Wl, l);  
  return 1.0 - gsl_cdf_gaussian_P(Z, 1.0);
}

double GWHL_VSH_getZConfidence_scaled(double Wl, double l, double N){
  double Z = GWHL_VSH_getZ_scaled(Wl, l, N);
  return 1.0 - gsl_cdf_gaussian_P(Z, 1.0);
}



double GWHL_VSH_getSurvivalForPower(double Wl, double l){
  return gsl_sf_gamma_inc_Q(l+0.5, Wl/2);
}

double GWHL_VSH_getSurvivalForPower_scaled(double Wl, double l, double n){
  return gsl_sf_gamma_inc_Q(n*(l+0.5), Wl/2);
}

double GWHL_VSH_getSurvivalForPowerSum(double Wl1, double Wl2, double l){
  return gsl_sf_gamma_inc_Q(2.0*l + 1, 0.5*(Wl1 + Wl2));
}


// internal function, converting vector to spherical coordinates
void GWHL_vec2bc(const double x, const double y, const double z, double *b, double *c){
  double no = sqrt(x*x + y*y + z*z);
  *b = GWHLMATHCONSTANTS_PI/2 - asin(z/no);
  *c = GWHLMATHCONSTANTS_PI - atan2(y, x);
  if(*c < 0){
    *c += 2 * GWHLMATHCONSTANTS_PI;
  }
}

/*
 * check if VSH expansion for l=2 is consistent with GW signal
 * for a gw signal the number returned by this function should be near zero
 */
double GWHL_VSH_checkVSHGWConsistency(const int maxl, double *vshcoffs) {
  if(maxl < 2) return DBL_MAX;
  //             0   1   2      3   4   5   6   7
  // vshcoffs = 10R 11R 11I    20R 21R 21I 22R 22I

  // temp array
  //0    1    2    3    4
  //fR20 fR21 fI21 fR22 fI22
  double f[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  
  // copy array
  memcpy(f,  &(vshcoffs[3]), 5);
  //OLD normalize array
  //GWHL_normalizeVectorNMaxNorm(f, 5);
  
  // new convert to uas
  GWHL_radToMicroasVec(f, 5);
  
  double one = pow(f[0], 3.0)/3.0;
  double two = sqrt(6) * (2*f[2]*f[4]*f[1] + f[3] * (f[1]*f[1] - f[2]*f[2]));
  double three = f[0] * (f[1]*f[1] + f[2]*f[2] - 2.0*(f[3]*f[3] + f[4]*f[4]));
  
  return one + two + three;
}

/**
 * returns the two points in celestial coordinates which define the symetry axis of the signal in l=2
 * return value is the smalles eigenvalue
 */
double GWHL_VSH_getSymetryAxisPointsOfL2(const int maxl, double *vshcoffs, double *ad, double *evsOut){
  if(maxl < 2) return DBL_MAX;
  //             0   1   2      3   4   5   6   7
  // vshcoffs = 10R 11R 11I    20R 21R 21I 22R 22I

  double b = 0, c = 0;
  double zptr[] = {0, 0, 1};
  double rotY[9];
  double rotZ[9];
  double rotYZ[9];
  double smallestEval = 0;
  
  double aCel = 0, bCel = 0;
  
  double cofsnew[8];
  // 3 from l = 1, 5 from l = 2  == 8
  memcpy(cofsnew, vshcoffs, 8 * sizeof(double));
  //NEW convert cofs to uas
  GWHL_radToMicroasVec(cofsnew, 8);
  
//     mat = {{-fR20 + Sqrt[6] fR22, -Sqrt[6] fI22, -Sqrt[6] fR21},
//            {-Sqrt[6] fI22, -fR20 - Sqrt[6] fR22, Sqrt[6] fI21},
//            {-Sqrt[6] fR21, Sqrt[6] fI21, 2 fR20}};
  double mat[] = {-cofsnew[3] + sqrt(6) * cofsnew[6], -sqrt(6) * cofsnew[7],               -sqrt(6) * cofsnew[4],
		   -sqrt(6) * cofsnew[7],               -cofsnew[3] - sqrt(6) * cofsnew[6], sqrt(6) * cofsnew[5],
		   -sqrt(6) * cofsnew[4],               sqrt(6) * cofsnew[5],                2.0 * cofsnew[3]};

    // OLD normalize the matrix with max norm
    // GWHL_normalizeVectorNMaxNorm(mat, 9);
 
  // set up GSL for eigenvalue computation, compute eigenvalues
  gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(3);
  gsl_vector *eval = gsl_vector_alloc(3);
  gsl_matrix *evec = gsl_matrix_alloc(3, 3);
  gsl_matrix_view gslmat = gsl_matrix_view_array(mat, 3, 3);
  gsl_eigen_symmv(&gslmat.matrix, eval, evec, ws);
  gsl_eigen_symmv_free(ws);
  
  // get smallest absolute eigenvalue
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
  smallestEval = gsl_vector_get(eval, 0);
  
  if(evsOut != NULL){
      evsOut[0] = gsl_vector_get(eval, 0);
      evsOut[1] = gsl_vector_get(eval, 1);
      evsOut[2] = gsl_vector_get(eval, 2);
      //fprintf(stderr, "all evals are: %+25.16e %+25.16e %+25.16e\n", evsOut[0],evsOut[1],evsOut[2]);
  }

  // convert vector associated to smallest eigenvalue to euler angles
  GWHL_vec2bc(gsl_matrix_get(evec, 0, 0), gsl_matrix_get(evec, 1, 0), gsl_matrix_get(evec, 2, 0), &b, &c);
  
  // use euler angles to get two rotational matrices
  GWHL_RotationMatrixY(b, rotY);
  GWHL_RotationMatrixZ(c, rotZ);
  
  // rotate the Z axis acording to the rotation matrices
  GWHL_matMatMul3d(rotZ, rotY, rotYZ);
  GWHL_matVecMul3d(rotYZ, zptr, zptr);
  
  // convert to celestial coordinates
  GWHL_getCelestCoordFromVec(zptr, &aCel, &bCel);

  ad[0] = fmod(aCel + 2 * GWHLMATHCONSTANTS_PI, 2 * GWHLMATHCONSTANTS_PI);
  ad[1] = bCel;
  
  ad[2] = fmod(ad[0] + GWHLMATHCONSTANTS_PI, 2 * GWHLMATHCONSTANTS_PI);
  ad[3] = -ad[1];
  
  // free structs
  gsl_vector_free(eval);
  gsl_matrix_free(evec);
  
  return smallestEval;
}


