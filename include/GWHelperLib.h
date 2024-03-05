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

#include <inttypes.h>
#include <omp.h>
#include <stdio.h>
#include "GWHL_Constants.h"
#include "GWHL_Quaternion.h"
#include "GWHL_StatsContainer.h"
#include "GWHL_SparseTypes.h"
#include "GWHL_Bitmasks.h"

/**
 * @brief degrees to radians
 * 
 * @param deg : degrees
 * @return double radians
 */
double GWHL_degToRad(double deg);


/**
 * @brief uas (micro-arcseconds) to radians
 * 
 * @param mas : uas
 * @return double radians
 */
double GWHL_microasToRad(double mas);

/**
 * @brief uas (micro-arcseconds) to radians
 * 
 * @param mas : uas
 * @return double radians
 */
double GWHL_uasToRad(double mas);

/**
 * @brief convert rad to uas (micro-arcseconds)
 * 
 * @param rad : input in rad
 * @return double output in uas
 */
double GWHL_radToMicroas(double rad);


/**
 * @brief convert rad to degrees
 * 
 * @param rad : input in rad
 * @return double output in deg
 */
double GWHL_radToDeg(double rad);

/**
 * @brief convert rad to degrees for every entry in vector
 * 
 * @param vec : input vector in rad
 * @param len : length of vector
 */
void   GWHL_degToRadVec(double *vec, unsigned long long int len);


void   GWHL_microasToRadVec(double *vec, unsigned long long int len);
void   GWHL_radToMicroasVec(double *vec, unsigned long long int len);
void   GWHL_radToDegVec(double *vec, unsigned long long int len);
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
void   GWHL_matTransposeNd(double *matA, long long int m, long long int n, double *matT);
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
void   GWHL_vecSqrtElemWise(double *vInOut, long long unsigned int len);
double GWHL_vecVecDotNd(double *vecA, double *vecB, long long int len);
void   GWHL_vecVecAddN(double *aIn, double *bInOut, uint64_t size);
void   GWHL_vecVecAddNCopy(double *aIn, double *bIn, double *cOut, uint64_t size);
void   GWHL_setArrayToXllu(long long unsigned int *array, long long unsigned int x, uint64_t arraylen);
void   GWHL_setArrayToXll(long long int *array, long long int x, uint64_t arraylen);

void   GWHL_normalizeVectorN(double *vec, unsigned long long int vl);
void   GWHL_normalizeVectorNMaxNorm(double *vec, unsigned long long int vl);
void   GWHL_normalizeVector3d(double *vec);
void   GWHL_zeroArray(double *array, uint64_t size);
void   GWHL_arrayAddN(double *aIn, double *bInOut, uint64_t size);
void   GWHL_addWithMask(double *vec, long long int vecl, double *result, long long int reslen, int *mask);
void   GWHL_setCompressedWithMask(double *vec, long long int vecl, double *result, long long int reslen, int *mask);
void   GWHL_setUncompressedWithMask(double *vec, long long int vecl, double *result, long long int reslen, int *mask);
void   GWHL_setArrayToXd(double *array, double x, uint64_t arraylen);
void   GWHL_setArrayToXi(int *array, int x, uint64_t arraylen);
void   GWHL_vecVecMulN(double *aInOut, double *bIn, uint64_t len);
void   GWHL_vecVecMulN_copy(double *a, double *b, double *cOut, uint64_t len);
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

void   GWHL_vecVecTMulScaleUpdate_sparse_ompsafe(GWHL_SparseVector_double *a, GWHL_SparseVector_double *b, double w, double *MOut, omp_lock_t *locks);
void   GWHL_vecVecTMulScaleUpdate_sparse(GWHL_SparseVector_double *a, GWHL_SparseVector_double *b, double w, double *MOut);
void   GWHL_vecScaleNdUpdate_sparse(GWHL_SparseVector_double *a, double b, long long unsigned int len, double *Vout);

/**
 * @brief returns the angular distance (rad) of two points in celestial coordinates (rad)
 */
double GWHL_SphCo_getAngularDistanceCelest(double a1, double d1, double a2, double d2);

/**
 * @brief returns the angle (rad) between the two 3d vectors vec1 and vec2
 */
double GWHL_SphCo_getAngularDistanceVec(const double *vec1, const double *vec2);


double GWHL_Helpers_getMaxAbsFromData(const double *datain, long long int len);
double GWHL_Helpers_getMinAbsFromData(const double *datain, long long int len);
double GWHL_Helpers_getMaxFromData(const double *datain, long long int len);
double GWHL_Helpers_getMinFromData(const double *datain, long long int len);
double GWHL_Helpers_getMaxAbsFromUpperRT(const double *matIn, long long int N);
double GWHL_Helpers_getMaxFromUpperRT(const double *matIn, long long int N);
double GWHL_Helpers_getMinFromUpperRT(const double *matIn, long long int N);
void   GWHL_Helpers_print_2files(FILE *f1, FILE *f2, char const *fmtstr, ...);

void   GWHL_MPIJobDist_getMyWorkChunk(long long int numWorkUnitsIN, long long int *firstIdxOUT, long long int *lastIdxOUT);

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
double GWHL_VSH_getZConfidence_scaled(double Wl, double l, double N);

double GWHL_VSH_getZ(double Wl, double l);
double GWHL_VSH_getZ_scaled(double Wl, double l, double N);

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
void   GWHL_VSH_getPowerOfOrdersNormalizedWOver(const int maxl, const double *vshCoefficients, const double *vshCoefficientsSigmas, double *result);

/**
 * Returns the probability that the given power for the given degree l
 * is caused by pure white noise.
 * 
 * @param Wl - normalized power
 * @param l - degree at which power occured
 */
double GWHL_VSH_getSurvivalForPower(double Wl, double l);

/**
 * Returns the probability that the given power combination (e.g. from S and T)
 * for the given degree l, is caused by pure white noise.
 * 
 * @param Wl1 - normalized power
 * @param Wl2 - normalized power
 * @param l - degree at which power occured
 */
double GWHL_VSH_getSurvivalForPowerSum(double Wl1, double Wl2, double l);

/**
 * Returns the probability that the given power combination (e.g. from S and T)
 * for the given degree l, is caused by pure white noise.
 * Equivalent with GWHL_VSH_getSurvivalForPowerSum, if Wl is sum of two Wl's and n = 2
 * 
 * @param Wl - normalized power (sum)
 * @param l - degree at which power occured
 * @param n - number of summed up Wl's (S+T -> n = 2)
 */
double GWHL_VSH_getSurvivalForPower_scaled(double Wl, double l, double n);

/**
 * gives the two points (in celestial coordinates) which define the symetry axis of the
 * VSH coefficents at l=2
 * @param maxl - maximum l for vshcofs
 * @param vshcoffs - array of all the VSH coefficients like (10R 11R 11I 20R 21R 21I 22R 22I ...)
 * @param ad - output coordinates like (alpha1 delta1 alpha2 delta2), must be pre-allocated
 * @returns smallest eigenvalue of the internal matrix
 */
double        GWHL_VSH_getSymetryAxisPointsOfL2(const int maxl, double *vshcoffs, double *ad, double *evsOut);
double        GWHL_VSH_checkVSHGWConsistency(const int maxl, double *vshcoffs);

//void   GWHL_LegendreP_getLegendreP(double x, long long int lmax, double *P);
long long int GWHL_LegendreP_lm2idx(long long int l, long long int m, long long int N);
long long int GWHL_LegendreP_idx2m(long long int idx, long long int N);
long long int GWHL_LegendreP_idx2l(long long int idx, long long int N);
void          GWHL_LegendreP_idx2lm(long long int idx, long long int N, long long int *l, long long int *m);
long long int GWHL_LegendreP_getOutputArraySize(long long int maxl);
void          GWHL_LegendreP_getAlphas(long long int lmax, double *alphas);
void          GWHL_LegendreP_getCosSinMap(double lambda, long long int lmax, double *cosmap, double *sinmap);
void          GWHL_LegendreP_getBetas(long long int lmax, double *betas);
void          GWHL_LegendreP_getAB(double x, long long int lmax, double *A, double *B);


void   GWHL_printMatrix(FILE *stream, const char *fmtString, const double *mat, const long long unsigned int sizeN, const long long unsigned int sizeM, int endWithNewline);
void   GWHL_printVector(FILE *stream, const char *fmtString, const double *vec, const long long unsigned int length, int endWithNewline);


double GWHL_Cholesky_decompLL(const double *A, long long int n, double *R);
void   GWHL_Cholesky_reduceRhsLL(const double *R, const double *b, long long int n, double *yRes);
void   GWHL_Cholesky_solveLL(const double *R, const double *y, long long int n, double *xRes);

/**
 * @brief Solves the system of equations Ax=b with the the method proposed by Lennart Lindegren in Lindegren et al. A&A 538, A78 (2012) "The astrometric core solution for the Gaia mission"
 * 
 * @param A : matrix A (M x M)
 * @param b : vector b (M)
 * @param nSize : size of matrix and array (M)
 * @param xRes p_xRes: pointer to pre-allocated memory for result (dimension also M)
 * @return double 
 */
double GWHL_Cholesky_SolveSystemLL(double *A, double *b, long long int nSize, double *xRes);


double GWHL_LinAlg_getRowSumNorm(double *A, long long unsigned int M, long long unsigned int N);
double GWHL_LinAlg_getColSumNorm(double *A, long long unsigned int M, long long unsigned int N);

/**
 * @brief mirror the matrix, transfers upper to lower triangular
 * 
 * @param A p_A: M x N matrix
 * @param M p_M: dimension M
 * @param N p_N: dimension N
 */
void   GWHL_LinAlg_mirrorMatrixUpperToLower(double *A, long long unsigned int M, long long unsigned int N);

/**
 * @brief get the euclidean norm sqrt(x1^2 + x2^2 + ... + x^i) of the vector
 * 
 * @param vIn : vector in
 * @param len : vector length
 * @return double norm
 */
double GWHL_vecNormNd(double *vIn, long long int len);

/**
* @brief returns continous time in double, starting from some arbitrary time
* 
* @return double
*/
double GWHL_getTime();

/**
* @brief writes like fprintf but uses two separate outputs
* 
* @param stream1 p_stream1: first stream to write to
* @param stream2 p_stream2: second stream to write to
* @param format p_format: format string
* @return int - 1 or 2 if failing to stream 1 or 2 fails, 0 on success
*/
int GWHL_fprintf2(FILE* stream1, FILE* stream2, const char *format, ...);


/**
 * @brief Set all bits to zero
 * 
 * @param bitmask bitmask to manipulate
 */
void GWHL_BitMasks_setAllZero(uint64_t *bitmask);


/**
 * @brief Set all bits to one
 * 
 * @param bitmask bitmask to manipulate
 */
void WHL_BitMasks_setAllOne(uint64_t *bitmask);


/**
 * @brief Set n'th bit to zero (counting starts at 0)
 * 
 * @param bitmask bitmask to manipulate
 * @param nth position of bit
 */
void GWHL_BitMasks_setNthToZero(uint64_t *bitmask, int nth);

/**
 * @brief Set n'th bit to one (counting starts at 0)
 * 
 * @param bitmask bitmask to manipulate
 * @param nth position of bit
 */
void GWHL_BitMasks_setNthToOne(uint64_t *bitmask, int nth);

/**
 * @brief Flipp n'th bit (counting starts at 0)
 * 
 * @param bitmask bitmask to manipulate
 * @param nth position of bit
 */
void GWHL_BitMasks_flipNth(uint64_t *bitmask, int nth);

/**
 * @brief get value of n'th bit (counting starts at 0)
 * 
 * @param bitmask bitmask to query
 * @param nth position of bit
 * @return bit value as integer
 */
int GWHL_BitMasks_getNth(uint64_t *bitmask, int nth);


/**
 * @brief check if all bits are zero
 * 
 * @param bitmask input bitmask
 * @return 1 if all are zero, 0 if not
 */
int GWHL_BitMasks_checkAllZero(uint64_t *bitmask);

/**
 * @brief check if all bits are one
 * 
 * @param bitmask input bitmask
 * @return 1 if all are ones, 0 if not
 */
int GWHL_BitMasks_checkAllOne(uint64_t *bitmask);

/**
 * @brief compare two bitmasks
 * 
 * @param bm1 : bitmask 1
 * @param bm2 : bitmask 2
 * @return int 1 if equal, 0 if not
 */
int GWHL_BitMasks_compare(uint64_t *bm1, uint64_t *bm2);

/**
 * @brief relatively efficient popcnt (hamming weight) for simple 64bit bitmasks
 * 
 * @param bitmask : bitmask
 * @return int number of bits which are true (1)
 */
int GWHL_BitMasks_popcount(uint64_t *bitmask);

/**
 * @brief initialize and allocates arbitrary length bitmask
 * 
 * @param numBits : number of bits
 * @param bm : pointer to bitmask variable (unallocated)
 * @return int 0 if successful allocated
 */
int GWHL_BitMasks_large_map_init(unsigned long long int numBits, gwhl_bitmask_long *bm);

/**
 * @brief set all bits in arbitrary length bitmask to zero (false)
 * 
 * @param bm : bitmask
 */
void GWHL_BitMasks_large_map_setAllZero(gwhl_bitmask_long *bm);

/**
 * @brief set all bits in arbitrary length bitmask to one (true)
 * 
 * @param bm : bitmask
 */
void GWHL_BitMasks_large_map_setAllOne(gwhl_bitmask_long *bm);

/**
 * @brief set n'th bit in arbitrary length bitmask to zero (false)
 * 
 * @param bm : bitmask
 * @param nth : element number to set (count starts at 0)
 */
void GWHL_BitMasks_large_map_setNthToZero(gwhl_bitmask_long *bm, unsigned long long int nth);

/**
 * @brief set n'th bit in arbitrary length bitmask to one (true)
 * 
 * @param bm : bitmask
 * @param nth : element number to set (count starts at 0)
 */
void GWHL_BitMasks_large_map_setNthToOne(gwhl_bitmask_long *bm, unsigned long long int nth);

/**
 * @brief flip n'th bit in arbitrary length bitmask
 * 
 * @param bm : bitmask
 * @param nth : element number to flip (count starts at 0)
 */
void GWHL_BitMasks_large_map_flipNth(gwhl_bitmask_long *bm, unsigned long long int nth);

/**
 * @brief get n'th bit in arbitrary length bitmask
 * 
 * @param bm : bitmask
 * @param nth : element number to return (count starts at 0)
 * @return int bit value of nth element (0 = false, 1 = true)
 */
int GWHL_BitMasks_large_map_getNth(gwhl_bitmask_long *bm, unsigned long long int nth);

/**
 * @brief compare two arbitrary length bitmasks
 * 
 * @param bm1 : bitmask 1
 * @param bm2 : bitmaks 2
 * @return int 0 = not equal content, 1 = equal content
 */
int GWHL_BitMasks_large_map_compare(gwhl_bitmask_long *bm1, gwhl_bitmask_long *bm2);

/**
 * @brief hamming weight or popcount of arbitrary length bitmask
 * 
 * @param bm1 : bitmask
 * @return uint64_t hamming weight (number of 1s in bit-vector)
 */
uint64_t GWHL_BitMasks_large_map_popcnt(gwhl_bitmask_long *bm1);

/**
 * @brief print an arbitrary length bitmask
 * 
 * @param bm : bitmask
 * @param fp : file pointer
 * @param newline : number of new lines after print
 */
void GWHL_BitMasks_large_map_print(gwhl_bitmask_long *bm, FILE *fp, int newline);

/**
 * @brief parallel quicksort implementation
 * 
 * @param data : array of doubles to sort (this will happen in-place)
 * @param len : length of array
 * @todo make workable for arbitrary data types
 */
void GWHL_Sorting_parallel_quicksort(double *data, long long unsigned int len);

/**
 * @brief set all entries in vec to their absolute values (in-place)
 * 
 * @param vec : input vector
 * @param len : vector length
 */
void GWHL_Lists_setAbsoluteEntries(double *vec, long long unsigned int len);

/**
 * @brief set all entries to its squared value (in-place replacement)
 * 
 * @param vec : input vector
 * @param len : vector length
 */
void GWHL_Lists_setSquaredEntries(double *vec, long long unsigned int len);

/**
 * @brief computes the quantile of the given data in the vector
 * 
 * @param vec array of doubles
 * @param len length of array
 * @param quantile quantile
 * @param inPlace proces in-place (will sort data in the array), meaningless if isSortedData = 1 
 * @param isSortedData indicated whether data is already sorted (will not be checked internaly)
 * @return the quantile asked for
 */
double GWHL_Stat_quantile(double *vec, long long unsigned int len, double quantile, int inPlace, int isSortedData);

/**
 * @brief returns the RSE (robust scatter estimate) of the data in te array
 * 
 * @param vec array with the data
 * @param len length of data
 * @param inPlace proces in-place (will sort data in the array)
 * @return RSE
 */
double GWHL_Stat_RSE(double *vec, long long unsigned int len, int inPlace, int isSortedData);

/**
 * @brief returns the RMS of data, weighted if necessary
 * 
 * @param vec data input
 * @param weights equally long array of weights (must be set to NULL if unweighted)
 * @param len length of arrays
 * @return RMS or weighted RMS value
 */
double GWHL_Stat_RMS(double *vec, double *weights, long long unsigned int len);

/**
 * @brief computes the weighted or un-weighted R2 value (Coefficient of determination)
 * 
 * @param dat data (observed)
 * @param mod prediction (values for data-points given by model)
 * @param weights weights (un-weighted if this is NULL)
 * @param len length of arrays
 * @return R2
 */
double GWHL_Stat_R2(double *dat, double *mod, double *weights, long long unsigned int len);

/**
 * @brief computes the weighted or un-weighted R2 value (Coefficient of determination) from residuals
 * 
 * @param dat data (observed)
 * @param residuals residuals from the fit
 * @param weights weights of observations
 * @param len length of arrays
 * @return double
 */
double GWHL_Stat_R2_Residuals(double *dat, double *residuals, double *weights, long long unsigned int len);
#endif
