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

/**
 * @file   MathConversions.c
 * @author Robin Geyer (robin.geyer@tu-dresden.de)
 * @date   October, 2015
 * @brief  mathematical conversion helper functions
 * @see    MathConversions.h
 */

#include <math.h>
#include "../include/GWHL_Constants.h"

/**
 * @brief   Degree to Radian conversion
 * @return  radian
 */
double GWHL_degToRad(double deg) {
    return (GWHLMATHCONSTANTS_PI * deg) / 180.0;
}

/**
* @brief inplace converts all values in vec from deg to rad
* 
* @param vec - input and output vector
* @param len - vector length
*/
void GWHL_degToRadVec(double *vec, unsigned long long int len){
    unsigned long long int i;
    for(i = 0; i < len; i++){
        vec[i] = GWHL_degToRad(vec[i]);
    }
}

/**
 * @brief   micro arc seconds to radian conversion
 * @return  radian
 */
double GWHL_microasToRad(double uas) {
    return (GWHLMATHCONSTANTS_PI * uas) / 648000000000.0;
}

double GWHL_uasToRad(double uas) {
    return (GWHLMATHCONSTANTS_PI * uas) / 648000000000.0;
}

/**
* @brief inplace converts all values in vec from uas to rad
* 
* @param vec - input and output vector
* @param len - vector length
*/
void GWHL_microasToRadVec(double *vec, unsigned long long int len){
    unsigned long long int i;
    for(i = 0; i < len; i++){
        vec[i] = GWHL_microasToRad(vec[i]);
    }
}


/**
 * @brief   Radian to micro arcsecond conversion
 * @return  micro arc seconds
 */
double GWHL_radToMicroas(double rad) {
    return (rad * 648000000000.0) / GWHLMATHCONSTANTS_PI;
}


/**
* @brief inplace converts all values in vec from rad to uas
* 
* @param vec - input and output vector
* @param len - vector length
*/
void GWHL_radToMicroasVec(double *vec, unsigned long long int len){
    unsigned long long int i;
    for(i = 0; i < len; i++){
        vec[i] = GWHL_radToMicroas(vec[i]);
    }
}


/**
 * @brief   radian to degree conversion
 * @return  degree
 */
double GWHL_radToDeg(double rad) {
    return (rad * 180.0) / GWHLMATHCONSTANTS_PI;
}

/**
* @brief inplace converts all values in vec from rad to deg
* 
* @param vec - input and output vector
* @param len - vector length
*/
void GWHL_radToDegVec(double *vec, unsigned long long int len){
    unsigned long long int i;
    for(i = 0; i < len; i++){
        vec[i] = GWHL_radToDeg(vec[i]);
    }
}
