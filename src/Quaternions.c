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
#include "../include/GWHL_Quaternion.h"
#include "../include/GWHL_Constants.h"

#define EPS 0.000001
void GWHL_Quaternion_getSetQuaternion(double x, double y, double z, double w, GWHLQuaternion *quat){
    quat->x = x;
    quat->y = y;
    quat->z = z;
    quat->w = w;
}


void GWHL_Quaternion_getAsDoubleArray(GWHLQuaternion *quat, double *out){
    out[0] = quat->x;
    out[1] = quat->y;
    out[2] = quat->z;
    out[3] = quat->w;
}

/**
 * Constructs and initializes a Normalized Rotation Quaternion from
 * an axis. The vector modulus is giving the angle value. Its orientation
 * defines the vector direction.
 * 
 * @param vec - the axis of the rotation.
 * @param quat - output quaternion
 */
void GWHL_Quaternion_getNormalizedRotationFromAxis(double *vec, GWHLQuaternion *quat){
    double mag;
    double amag;
    double angle;

    double ax = vec[0];
    double ay = vec[1];
    double az = vec[2];
    
    // triad normalization
    amag = sqrt((ax * ax) + (ay * ay) + (az * az));
    angle = amag / 2.0; // Quat = cos(theta/2) + sin(theta/2)(rotation_axis)

    mag = sin(angle);
    
    quat->w = cos(angle);
    
    amag = 1.0 / amag;
    
    quat->x = ax * amag * mag;
    quat->y = ay * amag * mag;
    quat->z = az * amag * mag;
}


void GWHL_Quaternion_getNormalizedRotationFromAxisAngle(double *vec, double angle, GWHLQuaternion *quat) {
    double mag;
    double amag;

    // Quat = cos(theta/2) + sin(theta/2)(rotation_axis)
    double ax = vec[0];
    double ay = vec[1];
    double az = vec[2];
    // triad normalization
    amag = sqrt((ax * ax) + (ay * ay) + (az * az));

    //TODO better throw error?
    if (amag < EPS) {
	quat->w = 1.0;
	quat->x = 0.0;
	quat->y = 0.0;
	quat->z = 0.0;
    } else {
	amag = 1.0 / amag;
	mag = sin(angle / 2.0);
	quat->w = cos(angle / 2.0);
	quat->x = ax * amag * mag;
	quat->y = ay * amag * mag;
	quat->z = az * amag * mag;
    }
}


void GWHL_Quaternion_conjugate(GWHLQuaternion *quat){
    quat->x = -quat->x;
    quat->y = -quat->y;
    quat->z = -quat->z;
    // w is not negated during conjugation
}


void GWHL_Quaternion_negate(GWHLQuaternion *quat){
    quat->x = -quat->x;
    quat->y = -quat->y;
    quat->z = -quat->z;
    quat->w = -quat->w;
}


void GWHL_Quaternion_multiply(GWHLQuaternion *q1, GWHLQuaternion *q2, GWHLQuaternion *quat){
    double x;
    double y;
    double z;
    double w;

    w = (quat->w * quat->w) - (quat->x * quat->x) - (quat->y * quat->y) - (quat->z * quat->z);
    x = ((quat->w * quat->x) + (quat->w * quat->x) + (quat->y * quat->z)) - (quat->z * quat->y);
    y = ((quat->w * quat->y) + (quat->w * quat->y)) - (quat->x * quat->z) + (quat->z * quat->x);
    z = ((quat->w * quat->z) + (quat->w * quat->z) + (quat->x * quat->y)) - (quat->y * quat->x);

    quat->x = x;
    quat->y = y;
    quat->z = z;
    quat->w = w;
}


void GWHL_Quaternion_inverse(GWHLQuaternion *quat){
    double x;
    double y;
    double z;
    double w;
    
    double invNormSqu;

    invNormSqu = 1.0 / ((quat->x * quat->x) + (quat->y * quat->y)
	  + (quat->z * quat->z) + (quat->w * quat->w));
    w = invNormSqu * quat->w;
    x = -invNormSqu * quat->x;
    y = -invNormSqu * quat->y;
    z = -invNormSqu * quat->z;
    
    quat->x = x;
    quat->y = y;
    quat->z = z;
    quat->w = w;
}



void GWHL_Quaternion_normalize(GWHLQuaternion *quat){
    double x;
    double y;
    double z;
    double w;
    
    double norm;

    norm = ((quat->x * quat->x) + (quat->y * quat->y) + (quat->z * quat->z) + (quat->w * quat->w));

    if (norm > 0.0) {
	    norm = 1.0 / sqrt(norm);
	    x = norm * quat->x;
	    y = norm * quat->y;
	    z = norm * quat->z;
	    w = norm * quat->w;
    } else {
	    x = 0.0;
	    y = 0.0;
	    z = 0.0;
	    w = 0.0;
    }
    
    quat->x = x;
    quat->y = y;
    quat->z = z;
    quat->w = w;
}


void GWHL_Quaternion_add(GWHLQuaternion *q1, GWHLQuaternion *q2, GWHLQuaternion *quat){
    quat->x = q1->x + q2->x;
    quat->y = q1->y + q2->y;
    quat->z = q1->z + q2->z;
    quat->w = q1->w + q2->w;
}


void GWHL_Quaternion_scale(GWHLQuaternion *q1, double scale, GWHLQuaternion *quat){
    quat->x = q1->x * scale;
    quat->y = q1->y * scale;
    quat->z = q1->z * scale;
    quat->w = q1->w * scale;
}


void GWHL_Quaternion_getDirectionCosineMatrix(GWHLQuaternion *q1, double *dcm){
    dcm[0] = (q1->x * q1->x) - (q1->y * q1->y) - (q1->z * q1->z) + (q1->w * q1->w);
    dcm[1] = 2.0 * ((q1->x * q1->y) + (q1->z * q1->w));
    dcm[2] = 2.0 * ((q1->x * q1->z) - (q1->y * q1->w));
    // second row
    dcm[3] = 2.0 * ((q1->x * q1->y) - (q1->z * q1->w));
    dcm[4] = ((-q1->x * q1->x) + (q1->y * q1->y)) - (q1->z * q1->z) + (q1->w * q1->w);
    dcm[5] = 2.0 * ((q1->y * q1->z) + (q1->x * q1->w));
    // third row
    dcm[6] = 2.0 * ((q1->x * q1->z) + (q1->y * q1->w));
    dcm[7] = 2.0 * ((q1->y * q1->z) - (q1->x * q1->w));
    dcm[8] = (-q1->x * q1->x) - (q1->y * q1->y) + (q1->z * q1->z) + (q1->w * q1->w);
}


void GWHL_Quaternion_computeEulerAngles(GWHLQuaternion *q1, double *thetaPhiPsi){
    double dcm[9] = {0,0,0, 0,0,0, 0,0,0};
    
    GWHL_Quaternion_getDirectionCosineMatrix(q1, dcm);
    
    // get angles
    thetaPhiPsi[0] = asin(dcm[6]);
    thetaPhiPsi[1] = atan2(-dcm[7], dcm[8]);
    thetaPhiPsi[2] = atan2(-dcm[3], dcm[0]);
    
    //normalize
    if(thetaPhiPsi[1] < 0.0){
      thetaPhiPsi[1] = thetaPhiPsi[1] + 2.0 * GWHLMATHCONSTANTS_PI;
    }
    if(thetaPhiPsi[2] < 0.0){
      thetaPhiPsi[2] = thetaPhiPsi[2] + 2.0 * GWHLMATHCONSTANTS_PI;
    }

}
