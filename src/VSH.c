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

static double dLast = 9999.123456789;
// do not address this arrays directly, use always gsl_sf_legendre_array_index(l,m)
static double *betas = NULL;
static double *A = NULL;
static double *B = NULL;


unsigned long long int GWHL_VSH_getNumCoefficients(unsigned long long int maxl){
    return 2 * (maxl*(maxl + 2));
}


// void GWHL_VSH_getBetaFactors(int maxl, double *factors){
//     int l;
//     for(l = 1; l <= maxl; l++){
//       factors[l-1] = 1.0/sqrt(l*(l+1.0));
//     }
// }
// 
// 
// void GWHL_VSH_getcm(double alpha, int maxl, double *cms){
//     int m;
//     for(m = 0; m <= maxl; m++){
//       cms[m] = cos(alpha * m);
//     }
// }
// 
// 
// /**
//  * get the factors s_m
//  */
// void GWHL_VSH_getsm(double alpha, int maxl, double *sms){
//     int m;
//     for(m = 0; m <= maxl; m++){
//       sms[m] = sin(alpha * m);
//     }
// }
// 
// 
// /**
//  * returns factor A(x)
//  * @param x - parameter of function, e.g. sin (delta)
//  */
// double GWHL_VSH_getFactorA(double x){
//   return sqrt(1.0 - pow(x, 2.0));
// }
// 
// 
// /**
//  * compute all factors B_lm(x)
//  */
// void GWHL_VSH_getFactorsB(double x, int maxl, double *bs){
//     int m;
//     for(m = 0; m <= maxl; m++){
//       bs[m] = m * (1.0 / sqrt(1.0 - pow(x, 2.0)));
//     }
// }
// 
// 
// /**
//  * Function which pastes together everything for the G part of the design equation
//  * @param maxl - l_max in VSH expansion
//  * @param Dplm - derivative of spherical harmonics normalized ((-1)^m \sqrt((2l + 1) * (l-m)! / (4 \pi) / (l+m)!) P_l^m(x))
//  * @param sms - s_m = sin m*alpha
//  * @param cms - c_m = cos m*alpha
//  * @param blm - beta_lm
//  * @param Afactors - A_lm(x)
//  * @param flmRI - result array, must be pre-allocated
//  */
// void GWHL_VSH_getFlmRI(int maxl, double *Dplm, double *sms, double *cms, double *blm, double Afactor, double *flmRI){
//   int i = 0;
//   int l, m;
//   for(l = 1; l <= maxl; l++){
//     for(m = 0; m <= l; m++){
//       if(m == 0) {
// // 	fprintf(stderr, "1/sqrt(%i(%i+1)) %+25.16e\n", l, l, blm[l-1]);
// // 	fprintf(stderr, "cos(%i*alpha) %+25.16e\n", m, cms[m]);
// // 	fprintf(stderr, "sqrt(1-sind^2) %+25.16e\n", Afactor);
// 	  flmRI[i] =        Dplm[gsl_sf_legendre_array_index(l,m)] * blm[l-1] * cms[m] * Afactor;
// 	  i++;
//       } else {
// 	  // real
// 	  flmRI[i] =  2.0 * Dplm[gsl_sf_legendre_array_index(l,m)] * blm[l-1] * cms[m] * Afactor;
// 	  i++;
// 	  //imaginary
// 	  flmRI[i] = -2.0 * Dplm[gsl_sf_legendre_array_index(l,m)] * blm[l-1] * sms[m] * Afactor;
// 	  i++;
//       }
//     }
//   }
// }
// 
// 
// /**
//  * Function which pastes together everything for the G part of the design equation
//  * @param maxl - l_max in VSH expansion
//  * @param plm - spherical harmonics normalized ((-1)^m \sqrt((2l + 1) * (l-m)! / (4 \pi) / (l+m)!) P_l^m(x))
//  * @param sms - s_m = sin m*alpha
//  * @param cms - c_m = cos m*alpha
//  * @param blm - beta_lm
//  * @param Bfactors - B_lm(x)
//  * @param glmRI - result array, must be pre-allocated
//  */
// void GWHL_VSH_getGlmRI(int maxl, double *plm, double *sms, double *cms, double *blm, double *Bfactors, double *glmRI){
// //   int gsl_sf_legendre_array_e(const gsl_sf_legendre_t norm, const size_t lmax, const double x, const double csphase, double result_array[]);
// 
//   int i = 0;
//   int l, m;
//   for(l = 1; l <= maxl; l++){
//     for(m = 0; m <= l; m++){
//       if(m == 0) {
// 	  glmRI[i] =  -1.0 *  plm[gsl_sf_legendre_array_index(l, m)] * blm[l-1] * sms[m] * Bfactors[m];
// 	  i++;
//       } else {
// 	  // real
// 	  glmRI[i] = -2.0 * plm[gsl_sf_legendre_array_index(l, m)] * blm[l-1] * sms[m] * Bfactors[m];
// 	  i++;
// 	  //imaginary
// 	  glmRI[i] = -2.0 * plm[gsl_sf_legendre_array_index(l, m)] * blm[l-1] * cms[m] * Bfactors[m];
// 	  i++;
//       }
//     }
//   }
// }
// 
// 
// 
// void GWHL_VSH_getGlmFlm(const int maxl, const double a, const double d, double *glmRI, double *flmRI){
//     //   double *glmRI = calloc(maxl*(maxl+2), sizeof(double));
//     //   double *flmRI = calloc(maxl*(maxl+2), sizeof(double));
//     double Dplm[gsl_sf_legendre_array_n(maxl)];
//     double plm[gsl_sf_legendre_array_n(maxl)];
//     // init static arrays if not done before
//     if(blm == NULL){     
// 	blm = calloc(maxl, sizeof(double));
// 	cms = calloc(maxl+1, sizeof(double));
// 	sms = calloc(maxl+1, sizeof(double));
// 	Bfactors = calloc(maxl+1, sizeof(double));
// 	Afactor = GWHL_VSH_getFactorA(sin(d));
//     }
//     
//     // don't recompute if delta is not changed, e.g. if HEALPIX properties are exploited
// //     if(dLast != d){
// 	gsl_sf_legendre_deriv_array_e(GSL_SF_LEGENDRE_SPHARM, maxl, sin(d), -1, plm, Dplm);
// 	GWHL_VSH_getFactorsB(sin(d), maxl, Bfactors);
// 	Afactor = GWHL_VSH_getFactorA(sin(d));
// 	dLast = d;
// //     }
//     
//     GWHL_VSH_getBetaFactors(maxl, blm);
//     GWHL_VSH_getcm(a, maxl, cms);
//     GWHL_VSH_getsm(a, maxl, sms);
//     
//     GWHL_VSH_getFlmRI(maxl, Dplm,  sms, cms, blm, Afactor,  flmRI);
//     GWHL_VSH_getGlmRI(maxl, plm, sms, cms, blm, Bfactors, glmRI);
// }


void GWHL_VSH_getGlmFlm(const int maxl, const double a, const double d, double *glmRI, double *flmRI){
    if(A == NULL){
	A = calloc(GWHL_LegendreP_getOutputArraySize(maxl), sizeof(double));
	B = calloc(GWHL_LegendreP_getOutputArraySize(maxl), sizeof(double));
	betas = calloc(GWHL_LegendreP_getOutputArraySize(maxl), sizeof(double));
    }
    
    double *sinmap = malloc((maxl + 1) * sizeof(double));
    double *cosmap = malloc((maxl + 1) * sizeof(double));
    long long int l, m, i;
    i = 0;
    
    if(dLast != d){
	GWHL_LegendreP_getBetas(maxl, betas);
	GWHL_LegendreP_getAB(sin(d), maxl, A, B);
	dLast = d;
    }
    
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
}



void GWHL_VSH_freeBuffers(){
    // free all structures if wished
    if(betas != NULL){
      free(betas);
      betas = NULL;
    }
    if(A != NULL){
      free(A);
      A = NULL;
      
    }
    if(B != NULL){
      free(B);
      B = NULL;
    }
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
	    nul = i;
	    result[l-1] += pow(vshCoefficients[i] / vshCoefficientsSigmas[nul], 2.0);
	    i++;
	} else {
	    // real
	    result[l-1] += pow(vshCoefficients[i],2.0) / (pow(vshCoefficientsSigmas[nul], 2.0)/2.0);
	    i++;
	    //imaginary
	    result[l-1] += pow(vshCoefficients[i],2.0)/ (pow(vshCoefficientsSigmas[nul], 2.0)/2.0);
	    i++;
	}
      }
    }
}


void GWHL_VSH_getPowerOfOrdersNormalizedWOver(const int maxl, const double *vshCoefficients, const double *vshCoefficientsSigmas, double *result){
    long long unsigned int i = 0;
    long long int l, m;
    long long int nul;
    double avgl[maxl];
    
    GWHL_setArrayToXd(avgl, 0.0, maxl);
    GWHL_setArrayToXd(result, 0.0, maxl);
    GWHL_VSH_getPowerOfOrdersUnNormalized(maxl, vshCoefficients, result);
    
    
    for(l = 1; l <= maxl; l++){
      for(m = 0; m <= l; m++){
	if(m == 0) {
	    nul = i;
	    avgl[l-1] += pow(vshCoefficientsSigmas[nul], 2.0);
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


/**
 * returns the noise probability as described in Mignard Klioner 2012 Equation (85)
 * we return the probability here as Z is from a normal distribution with zero mean and unit variance
 */
double GWHL_VSH_getZConfidence(double Wl, double l){
  double n = 2 * l + 1;
  double Z = sqrt(4.5*n)*(pow((Wl/n),1.0/3.0)-(1.0-2.0/(9*n)));
  
  return 1.0 - gsl_cdf_gaussian_P(Z, 1.0);
}

/**
 * returns the noise probability as described in Mignard Klioner 2012 Equation (85)
 * we return only Z itself here 
 */
double GWHL_VSH_getZ(double Wl, double l){
  double n = 2 * l + 1;
  double Z = sqrt(4.5*n)*(pow((Wl/n),1.0/3.0)-(1.0-2.0/(9*n)));
  
  return Z;
}

double GWHL_VSH_getSurvivalForPower(double Wl, double l){
  return gsl_sf_gamma_inc_Q(l+0.5, Wl/2);
}

double GWHL_VSH_getSurvivalForPowerSum(double Wl1, double Wl2, double l){
  return gsl_sf_gamma_inc_Q(1+2*l, (Wl1+Wl2)/2);
}


// internal function
void GWHL_vec2bc(const double x, const double y, const double z, double *b, double *c){
  double no = sqrt(x*x + y*y + z*z);
  *b = GWHLMATHCONSTANTS_PI/2 - asin(z/no);
  *c = GWHLMATHCONSTANTS_PI - atan2(y, x);
  if(*c < 0){
    *c += 2 * GWHLMATHCONSTANTS_PI;
  }
}


/**
 * returns the two points in celestial coordinates which define the symetry axis of the signal in l=2
 */
void GWHL_VSH_getSymetryAxisPointsOfL2(const int maxl, double *vshcoffs, double *ad){
  if(maxl < 2) return;
  //             0   1   2    3   4   5   6   7
  // vshcoffs = 10R 11R 11I  20R 21R 21I 22R 22I
  double sqrt6 = sqrt(6);
  double b = 0, c = 0;
  double zptr[] = {0, 0, 1};
  double rotY[9];
  double rotZ[9];
  double rotYZ[9];
  
  double aCel = 0, bCel = 0;
  double mat[] = {-vshcoffs[3] + sqrt6 * vshcoffs[6], -sqrt6 * vshcoffs[7],               -sqrt6 * vshcoffs[4],
		   -sqrt6 * vshcoffs[7],               -vshcoffs[3] - sqrt6 * vshcoffs[6], sqrt6 * vshcoffs[5],
		   -sqrt6 * vshcoffs[4],               sqrt6 * vshcoffs[5],                2.0 * vshcoffs[3]};
		   
		   
//  GWHL_printMatrix(stderr, "%+25.16e", mat, 3, 3, 1);
 
  gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(3);
  gsl_vector *eval = gsl_vector_alloc(3);
  gsl_matrix *evec = gsl_matrix_alloc(3, 3);
  gsl_matrix_view gslmat = gsl_matrix_view_array(mat, 3, 3);
  
  gsl_eigen_symmv(&gslmat.matrix, eval, evec, ws);
  gsl_eigen_symmv_free(ws);
  
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
  
  GWHL_vec2bc(gsl_matrix_get(evec, 0, 0), gsl_matrix_get(evec, 1, 0), gsl_matrix_get(evec, 2, 0), &b, &c);
  
//   fprintf(stderr, "\n%+25.16e  %+25.16e  %+25.16e\n\n%+25.16e  %+25.16e  %+25.16e\n%+25.16e  %+25.16e  %+25.16e\n%+25.16e  %+25.16e  %+25.16e\n\n",
// 	  gsl_vector_get(eval, 0), gsl_vector_get(eval, 1), gsl_vector_get(eval, 2),
// 	  gsl_matrix_get(evec, 0, 0), gsl_matrix_get(evec, 0, 1), gsl_matrix_get(evec, 0, 2),
// 	  gsl_matrix_get(evec, 1, 0), gsl_matrix_get(evec, 1, 1), gsl_matrix_get(evec, 1, 2),
// 	  gsl_matrix_get(evec, 2, 0), gsl_matrix_get(evec, 2, 1), gsl_matrix_get(evec, 2, 2)
// 	);
  
  gsl_vector_free(eval);
  gsl_matrix_free(evec);
  
  GWHL_RotationMatrixY(b, rotY);
  GWHL_RotationMatrixZ(c, rotZ);
  
//   fprintf(stderr, "b = %+25.16e\t c = %+25.16e\n\n", b, c);
  
  GWHL_matMatMul3d(rotZ, rotY, rotYZ);
  GWHL_matVecMul3d(rotYZ, zptr, zptr);
  
  GWHL_getCelestCoordFromVec(zptr, &aCel, &bCel);

  ad[0] = fmod(aCel + 2 * GWHLMATHCONSTANTS_PI, 2 * GWHLMATHCONSTANTS_PI);
  ad[1] = bCel;
  
  ad[2] = fmod(ad[0] + GWHLMATHCONSTANTS_PI, 2 * GWHLMATHCONSTANTS_PI);
  ad[3] = -ad[1];
}

