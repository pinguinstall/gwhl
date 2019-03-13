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
#include "../include/GWHL_Constants.h"

/**
 * Rotation matrix rotating theta radian around the arbitrary vector u=(x,y,z)^T
 */
void GWHL_GeneralRotationMatrix(double theta, double *u, double *R){
  double x = u[0];
  double y = u[1];
  double z = u[2];

  double yxzpow = pow(x,2) + pow(y,2) + pow(z,2);
  double normxyz = sqrt(yxzpow);
  x = x / normxyz;
  y = y / normxyz;
  z = z / normxyz;
  
  double costheta = cos(theta);
  double sintheta = sin(theta);
  
  R[0] = costheta + x * x * (1-costheta);
  R[1] = x * y * (1-costheta) - z * sintheta;
  R[2] = x * z * (1-costheta) + y * sintheta;

  R[3] = y * x * (1-costheta) + z * sintheta;
  R[4] = costheta + y * y * (1-costheta);
  R[5] = y * z * (1-costheta) - x * sintheta;

  R[6] = z * x * (1-costheta) - y * sintheta;
  R[7] = z * y * (1-costheta) + x * sintheta;
  R[8] = costheta + z * z * (1-costheta);

}


/**
 * Rotation matrix around X
 */
void GWHL_RotationMatrixX(double theta, double *R){
  R[0] = 1.0;
  R[1] = 0.0;
  R[2] = 0.0;

  R[3] = 0.0;
  R[4] = cos(theta);
  R[5] = sin(theta);

  R[6] = 0.0;
  R[7] = -sin(theta);
  R[8] = cos(theta);
}


/**
 * Rotation matrix around Y
 */
void GWHL_RotationMatrixY(double theta, double *R){
  R[0] = cos(theta);
  R[1] = 0.0;
  R[2] = -sin(theta);

  R[3] = 0.0;
  R[4] = 1.0;
  R[5] = 0.0;

  R[6] = sin(theta);
  R[7] = 0.0;
  R[8] = cos(theta);
}


/**
 * Rotation matrix around Z
 */
void GWHL_RotationMatrixZ(double theta, double *R){
  R[0] = cos(theta);
  R[1] = sin(theta);
  R[2] = 0.0;

  R[3] = -sin(theta);
  R[4] = cos(theta);
  R[5] = 0.0;

  R[6] = 0.0;
  R[7] = 0.0;
  R[8] = 1.0;
}
