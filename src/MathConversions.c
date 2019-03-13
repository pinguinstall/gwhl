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
 * @brief   micro arc seconds to radian conversion
 * @return  radian
 */
double GWHL_microasToRad(double mas) {
    return (GWHLMATHCONSTANTS_PI * mas) / 648000000000.0;
}


/**
 * @brief   Radian to micro arcsecond conversion
 * @return  micro arc seconds
 */
double GWHL_radToMicroas(double rad) {
    return (rad * 648000000000.0) / GWHLMATHCONSTANTS_PI;
}


/**
 * @brief   radian to degree conversion
 * @return  degree
 */
double GWHL_radToDeg(double rad) {
    return (rad * 180.0) / GWHLMATHCONSTANTS_PI;
}
