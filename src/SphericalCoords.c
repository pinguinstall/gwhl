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
#include <inttypes.h>
#include "LinearAlgebra.h"
#include "../include/GWHL_Constants.h"

/**
 * compute the local triad (p q r) in spherical coordinates
 */
void GWHL_localTriad(double alpha, double delta, double *p, double *q, double *r) {
    double ca = cos(alpha);
    double sa = sin(alpha);
    double cd = cos(delta);
    double sd = sin(delta);

    p[0] = -sa;
    p[1] = ca;
    p[2] = 0.0;

    q[0] = -sd * ca;
    q[1] = -sd * sa;
    q[2] = cd;

    r[0] = cd * ca;
    r[1] = cd * sa;
    r[2] = sd;

}

/**
 * compute the local triad (p q 0) in spherical coordinates
 */
void GWHL_localTriadPQOnly(double alpha, double delta, double *p, double *q) {
    double ca = cos(alpha);
    double sa = sin(alpha);
    double cd = cos(delta);
    double sd = sin(delta);

    p[0] = -sa;
    p[1] = ca;
    p[2] = 0.0;

    q[0] = -sd * ca;
    q[1] = -sd * sa;
    q[2] = cd;
}

/**
 * @brief   compute unit vector in the celestial alpha direction
 * @param a alpha (right ascension)
 * @param d delta (declination)
 * @param [out] EHVa unit vector in alpha
 */
void GWHL_getEHValpha(double a, double d, double *EHVa) {
    EHVa[0] = -sin(a);
    EHVa[1] = cos(a);
    EHVa[2] = 0.0;
}


/**
 * @brief   compute unit vector in the celestial delta direction
 * @param a alpha (right ascension)
 * @param d delta (declination)
 * @param [out] EHVd unit vector in delta
 */
void GWHL_getEHVdelta(double a, double d, double *EHVd) {
    EHVd[0] = -cos(a) * sin(d);
    EHVd[1] = -sin(a) * sin(d);
    EHVd[2] = cos(d);
}

/**
 * @brief   compute unit vector in the celestial u direction
 * @param a alpha (right ascension)
 * @param d delta (declination)
 * @param [out] EHVu unit vector in u
 */
void GWHL_getEHVu(double a, double d, double *EHVu) {
    EHVu[0] = cos(a) * cos(d);
    EHVu[1] = sin(a) * cos(d);
    EHVu[2] = sin(d);
}

/**
 * @brief compute unit vector in the celestial alpha and delta direction
 * @param a alpha (right ascension)
 * @param d delta (declination)
 * @param [out] EHVa unit vector in alpha
 * @param [out] EHVd unit vector in alpha
 *
 */
void GWHL_getEHVFromAD(double a, double d, double *EHVa, double *EHVd) {
    EHVa[0] = -sin(a);
    EHVa[1] = cos(a);
    EHVa[2] = 0.0;

    EHVd[0] = -cos(a) * sin(d);
    EHVd[1] = -sin(a) * sin(d);
    EHVd[2] = cos(d);
}

/**
 * @brief compute unit vector in the celestial alpha direction
 * @param ad vector with alpha = ad[0] and delta = ad[1]
 * @param [out] EHVa unit vector in alpha
 * @param [out] EHVd unit vector in delta
 *
 */
void GWHL_getEHVFromAD2(const double *ad, double *EHVa, double *EHVd) {
    EHVa[0] = -sin(ad[0]);
    EHVa[1] = cos(ad[0]);
    EHVa[2] = 0.0;

    EHVd[0] = -cos(ad[0]) * sin(ad[1]);
    EHVd[1] = -sin(ad[0]) * sin(ad[1]);
    EHVd[2] = cos(ad[1]);
}


/**
 * @brief  get celestial coordinates from pointing vector
 * @param vec pointing vector
 * @param [out] a right ascension (alpha)
 * @param [out] d declination (delta)
 */
void GWHL_getCelestCoordFromVec(const double *vec, double *a, double *d) {
    // if pole
    double xy = sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
    if(xy <= 1.0e-16){
	*a = 0.0;
	*d = 0.5 * GWHLMATHCONSTANTS_PI * ((vec[2] > 0) - (vec[2] < 0)); // signum function
	return;
    }
    
    *a = atan2(vec[1], vec[0]);
    
    // normalize negative values
    if(*a < 0.0){
	*a += (2.0 * GWHLMATHCONSTANTS_PI);
    }
    
    *d = atan2(vec[2], xy);
}


/**
 * @brief get celestial coordinates from pointing vector
 * @param vec - pointing vector
 * @param [out] ad - 2D vector [0] = alpha, [1] = delta
 */
void GWHL_getCelestCoordFromVec2(const double *vec, double *ad) {   
    // if pole
    double xy = sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
    if(xy <= 1.0e-16){
	ad[0] = 0.0;
	ad[1] = 0.5 * GWHLMATHCONSTANTS_PI * ((vec[2] > 0) - (vec[2] < 0)); // signum function
	return;
    }
    
    ad[0] = atan2(vec[1], vec[0]);
    
    // normalize negative values
    if(ad[0] < 0.0){
	ad[0] += (2.0 * GWHLMATHCONSTANTS_PI);
    }
    
    ad[1] = atan2(vec[2], xy);
}


/**
 * @brief given the local shift in alpha* and delta, this function restores the euclidean tangend vector representing the shift
 */
void GWHL_restoreShiftVecFromComponents(const double alpha, const double delta, const double dA, const double dD, double *deltaU){
  double EHVa[3];
  double EHVd[3];
  
  GWHL_getEHValpha(alpha, delta, EHVa);
  GWHL_getEHVdelta(alpha, delta, EHVd);
  
  GWHL_vecScale3d(EHVa, dA, EHVa);
  GWHL_vecScale3d(EHVd, dD, EHVd);
  
  deltaU[0] = EHVa[0] + EHVd[0];
  deltaU[1] = EHVa[1] + EHVd[1];
  deltaU[2] = EHVa[2] + EHVd[2];
}



double GWHL_SphCo_getAngularDistanceCelest(double a1, double d1, double a2, double d2){
  return acos(sin(d1)*sin(d2) + cos(d1)*cos(d2)*cos(a1-a2));
}


double GWHL_SphCo_getAngularDistanceVec(const double *vec1, const double *vec2){
  double v1n[3] = {vec1[0], vec1[1], vec1[2]};
  double v2n[3] = {vec2[0], vec2[1], vec2[2]};
  GWHL_normalizeVectorN(v1n,3);
  GWHL_normalizeVectorN(v2n,3);
  
  return acos(GWHL_vecVecDot3d(v1n, v2n));
}
