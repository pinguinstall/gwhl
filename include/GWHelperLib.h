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

#ifndef ____GWHELPERLIB__H
#define ____GWHELPERLIB__H

#include <mpi.h>
#include <inttypes.h>
#include "GWHL_Constants.h"
#include "GWHL_Quaternion.h"
#include "GWHL_StatsContainer.h"

char*  GWHL_MPIChunkReadFile(MPI_File *fp, int overlap, MPI_Comm com);
void   GWHL_MPIgetFileChunkBorders(MPI_File *fp, int overlap, MPI_Comm com, long long int *start, long long int *end, long long int *size);
double GWHL_degToRad(double deg);
double GWHL_microasToRad(double mas);
double GWHL_radToMicroas(double rad);
double GWHL_radToDeg(double rad);
void   GWHL_getEHValpha(double a, double d, double *EHVa);
void   GWHL_getEHVdelta(double a, double d, double *EHVd);
void   GWHL_getEHVu(double a, double d, double *EHVu);
void   GWHL_getCelestCoordFromVec(const double *vec, double *a, double *d);
void   GWHL_getCelestCoordFromVec2(const double *vec, double *ad);
void   GWHL_getLogscaleList(double start, double end, long long int numPoints, double *result);
long long int GWHL_popcnt(int *vec, long long unsigned int n);
void   GWHL_getLinearList(double start, double end, long long int numPoints, double *result);
void   GWHL_getLinearList_int(long long int start, long long int end, long long int numPoints, long long int *result);
void   GWHL_getLinearCenterList(double start, double end, long long int numBuckets, double *result);
void   GWHL_matScale3d(double *matA, double a, double *matR);
void   GWHL_matMatAdd3d(double *matA, double *matB, double *matR);
void   GWHL_vecVecAddScaleN(double *aIn, double *bInOut, double s, uint64_t size);
void   GWHL_matVecMul3d(double *matA, double *vecB, double *vecR);
double GWHL_vecVecDot3d(double *vecA, double *vecB);
void   GWHL_matVecMulNd(double *matA, double *vecB, long long int nColsNVec, long long int nRows, double *vecR);
void   GWHL_matMatMul3d(double *matA, double *matB, double *matR);
void   GWHL_matTranspose3d(double *matA, double *matT);
void  GWHL_matTransposeNd(double *matA, long long int m, long long int n, double *matT);
void   GWHL_vecVecAdd3d(double *vecA, double *vecB, double *vecR);
void   GWHL_vecScale3d(double *vecA, double a, double *vecR);
void   GWHL_getEHVFromAD(double a, double d, double *EHVa, double *EHVd);
void   GWHL_getEHVFromAD2(const double *ad, double *EHVa, double *EHVd);
void   GWHL_localTriad(double alpha, double delta, double *p, double *q, double *r);
void   GWHL_localTriadPQOnly(double alpha, double delta, double *p, double *q);
void   GWHL_crossProductVec3d(double *vecA, double *vecB, double *cross);
void   GWHL_vecVecTMul(double *vecAIn, double *vecBIn, double *matOut, long long unsigned int size);
void   GWHL_vecVecTMulUpdate(double *vecAIn, double *vecBIn, double *matInOut, long long unsigned int size);
void   GWHL_vecScaleNdUpdate(double *vecA, double a, unsigned long long int length, double *vecRupdated);
void   GWHL_vecScaleNd(double *vecA, double a, unsigned long long int length, double *vecR);
void   GWHL_restoreShiftVecFromComponents(const double alpha, const double delta, const double dA, const double dD, double *deltaU);
double GWHL_vecVecDotNd(double *vecA, double *vecB, long long int len);
void   GWHL_vecVecAddN(double *aIn, double *bInOut, uint64_t size);

void   GWHL_normalizeVectorN(double *vec, int vl);
void   GWHL_zeroArray(double *array, uint64_t size);
void   GWHL_arrayAddN(double *aIn, double *bInOut, uint64_t size);
void   GWHL_addWithMask(double *vec, int vecl, double *result, int reslen, int *mask);
void   GWHL_setCompressedWithMask(double *vec, int vecl, double *result, int reslen, int *mask);
void   GWHL_setUncompressedWithMask(double *vec, int vecl, double *result, int reslen, int *mask);
void   GWHL_setArrayToXd(double *array, double x, uint64_t arraylen);
void   GWHL_setArrayToXi(int *array, int x, uint64_t arraylen);
void   GWHL_vecVecMulN(double *a, double *b, uint64_t len);
void   GWHL_vecVecDivN(double *a, double *b, uint64_t len);
void   GWHL_vecDivN(double *a, double b, uint64_t len);
void   GWHL_vecVecTMulScaleUpdate(double *vecAIn, double *vecBIn, double weight, double *matInOut, long long unsigned int size);

void   GWHL_Quaternion_getSetQuaternion(double x, double y, double z, double w, GWHLQuaternion *quat);
void   GWHL_Quaternion_getAsDoubleArray(GWHLQuaternion *quat, double *out);
void   GWHL_Quaternion_getNormalizedRotationFromAxis(double *vec, GWHLQuaternion *quat);
void   GWHL_Quaternion_getNormalizedRotationFromAxisAngle(double *vec, double angle, GWHLQuaternion *quat);
void   GWHL_Quaternion_conjugate(GWHLQuaternion *quat);
void   GWHL_Quaternion_negate(GWHLQuaternion *quat);
void   GWHL_Quaternion_multiply(GWHLQuaternion *q1, GWHLQuaternion *q2, GWHLQuaternion *quat);
void   GWHL_Quaternion_inverse(GWHLQuaternion *quat);
void   GWHL_Quaternion_normalize(GWHLQuaternion *quat);
void   GWHL_Quaternion_add(GWHLQuaternion *q1, GWHLQuaternion *q2, GWHLQuaternion *quat);
void   GWHL_Quaternion_scale(GWHLQuaternion *q1, double scale, GWHLQuaternion *quat);
void   GWHL_Quaternion_getDirectionCosineMatrix(GWHLQuaternion *q1, double *dcm);
void   GWHL_Quaternion_computeEulerAngles(GWHLQuaternion *q1, double *thetaPhiPsi);

int    GWHL_compareInt(const void *a, const void *b);
int    GWHL_compareLong(const void *a, const void *b);
int    GWHL_compareDouble(const void *a, const void *b);
int    GWHL_compareFloat(const void *a, const void *b);
int    GWHL_compareUInt64(const void *a, const void *b);

void   GWHL_GeneralRotationMatrix(double theta, double *u, double *R);
void   GWHL_RotationMatrixX(double theta, double *R);
void   GWHL_RotationMatrixY(double theta, double *R);
void   GWHL_RotationMatrixZ(double theta, double *R);

/**
 * @brief returns the angular distance (rad) of two points in celestial coordinates (rad)
 */
double GWHL_SphCo_getAngularDistanceCelest(double a1, double d1, double a2, double d2);

/**
 * @brief returns the angle (rad) between the two 3d vectors vec1 and vec2
 */
double GWHL_SphCo_getAngularDistanceVec(const double *vec1, const double *vec2);


double GWHL_Helpers_getMaxAbsFromData(const double *datain, long long int len);
double GWHL_Helpers_getMaxAbsFromUpperRT(const double *matIn, long long int N);


void GWHL_MPIJobDist_getMyWorkChunk(long long int numWorkUnitsIN, long long int *firstIdxOUT, long long int *lastIdxOUT);

/**
 * Returns the left-hand side of the VSH design equation,
 * the conventions here are the same as in GAIA-CA-TN-LO-SK-016-1
 * @param maxl - l_max in VSH expansion
 * @param a - right ascention alpha
 * @param d - declination delta
 * @param glmRI - the G part like 10 11R 11I 20 21R 21I .... (must be pre-allocated)
 * @param flmRI - the F part like 10 11R 11I 20 21R 21I .... (must be pre-allocated)
 */
void   GWHL_VSH_getGlmFlm(const int maxl, const double a, const double d, double *glmRI, double *flmRI);

/**
 * Frees and deallocate all internal buffers used by the computation of
 * G and F
 */
void   GWHL_VSH_freeBuffers();
unsigned long long int GWHL_VSH_getNumCoefficients(unsigned long long int maxl);

/**
 * Returns the normalized power per VSH degree l. The normalization used is the variance 
 * per coefficient. See Mignard, Klioner; A&A 2012 547 A59; Equation 86
 * 
 * This routine is made in the way that either S or T coefficients are handled.
 * Care must be taken that the right data is handed over  
 * 
 * @param maxl - l_max of the VSH expansion
 * @param vshCoefficients - unnormalized VSH coefficients
 * @param vshCoefficientsSigmas - sigmas for the coefficents
 * @param result - normalized powers, of length maxl (must be pre-allocated)
 */
void   GWHL_VSH_getPowerOfOrdersNormalizedWTilde(const int maxl, const double *vshCoefficients, const double *vshCoefficientsSigmas, double *result);
void   GWHL_VSH_printVSHCoefficients(FILE *fp, double *coeffs, int maxl);

/**
 * returns the x*sigma value of the power
 * @param Wl - power of the degree
 * @param l - degree l in VSH
 */
double GWHL_VSH_getZConfidence(double Wl, double l);

double GWHL_VSH_getZ(double Wl, double l);

/**
 * Returns the unnomalized powers as defined in Mignard, Klioner; A&A 2012 547 A59; Equation 76
 * @param maxl - l_max of the VSH expansion
 * @param vshCoefficients - unnormalized VSH coefficients
 * @param result - normalized powers, of length maxl (must be pre-allocated)
 */
void   GWHL_VSH_getPowerOfOrdersUnNormalized(const int maxl, const double *vshCoefficients, double *result);

/**
 * Returns the normalized power per VSH degree l. The normalization used is the variance 
 * for coefficient l0. See Mignard, Klioner; A&A 2012 547 A59; Equation 83
 * 
 * This routine is made in the way that either S or T coefficients are handled.
 * Care must be taken that the right data is handed over  
 * 
 * @param maxl - l_max of the VSH expansion
 * @param vshCoefficients - unnormalized VSH coefficients
 * @param vshCoefficientsSigmas - sigmas for the coefficents
 * @param result - normalized powers, of length maxl (must be pre-allocated)
 */
void   GWHL_VSH_getPowerOfOrdersNormalizedW(const int maxl, const double *vshCoefficients, const double *vshCoefficientsSigmas, double *result);
void GWHL_VSH_getPowerOfOrdersNormalizedWOver(const int maxl, const double *vshCoefficients, const double *vshCoefficientsSigmas, double *result);

/**
 * Returns the probability that the given power for the given degree l
 * is caused by pure white noise.
 * 
 * @param Wl - normalized power
 * @param l - degree at which power occured
 */
double GWHL_VSH_getSurvivalForPower(double Wl, double l);


double GWHL_VSH_getSurvivalForPowerSum(double Wl1, double Wl2, double l);

/**
 * gives the two points (in celestial coordinates) which define the symetry axis of the
 * VSH coefficents at l=2
 * @param maxl - maximum l for vshcofs
 * @param vshcoffs - array of all the VSH coefficients like (10R 11R 11I 20R 21R 21I 22R 22I ...)
 * @param ad - output coordinates like (alpha1 delta1 alpha2 delta2), must be pre-allocated
 */
void GWHL_VSH_getSymetryAxisPointsOfL2(const int maxl, double *vshcoffs, double *ad);

void GWHL_LegendreP_getLegendreP(double x, long long int lmax, double *P);
long long int GWHL_LegendreP_lm2idx(long long int l, long long int m, long long int N);
long long int GWHL_LegendreP_idx2m(long long int idx, long long int N);
long long int GWHL_LegendreP_idx2l(long long int idx, long long int N);
void          GWHL_LegendreP_idx2lm(long long int idx, long long int N, long long int *l, long long int *m);
long long int GWHL_LegendreP_getOutputArraySize(long long int maxl);
void          GWHL_LegendreP_getAlphas(long long int lmax, double *alphas);
void          GWHL_LegendreP_getCosSinMap(double lambda, long long int lmax, double *cosmap, double *sinmap);
void          GWHL_LegendreP_getBetas(long long int lmax, double *betas);
void          GWHL_LegendreP_getAB(double x, long long int lmax, double *A, double *B);


void GWHL_printMatrix(FILE *stream, const char *fmtString, const double *mat, const long long unsigned int sizeN, const long long unsigned int sizeM, int endWithNewline);
void GWHL_printVector(FILE *stream, const char *fmtString, const double *vec, const long long unsigned int length, int endWithNewline);


double GWHL_Cholesky_decompLL(const double *A, long long int n, double *R);
void GWHL_Cholesky_reduceRhsLL(const double *R, const double *b, long long int n, double *yRes);
void GWHL_Cholesky_solveLL(const double *R, const double *y, long long int n, double *xRes);
/**
 * Solves the system of equations Ax=b with the the method proposed by Lennart Lindegren in
 * Lindegren et al. A&A 538, A78 (2012) "The astrometric core solution for the Gaia mission"
 * 
 * 
 */
double GWHL_Cholesky_SolveSystemLL(double *A, double *b, long long int nSize, double *xRes);

GWHL_StatsContainer* GWHL_getStandardStatsFromDoubleData(double *datVec, unsigned long long int size, int inPlace);



#endif
