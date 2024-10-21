#pragma once
/**
 * @file This header file contains definitions and prototypes for the functions 
 * contained in the LST&E C program library.
 * @author Robert Freepartner, JaDa Systems/Raytheon/JPL
 * @copyright (c) Copyright 2015, Jet Propulsion Laboratories, Pasadena, CA
 */

#include "bool.h"
#include "hdf.h"
#include <math.h>
#include <stdlib.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif
// Natural logarithm base
#define E 2.71828182845904523536028747135266249775724709369995
// Value for floating point equality test
#define epsilon 1E-9
// Declare local not-a-number value as nan, but avoids redefining the builtin 
// function nan() by giving it the actual name _NaN_.
#define nan _NaN_
#define _NaN_ -32768
#define REAL_NAN (0.0/0.0)

#ifndef ABS
#define ABS(a) ((a) < 0 ? -(a) : (a))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

/// @constant Hard limits for various types of profile data:
/// Height Profile (m)
#define hard_g_min 0.0
#define hard_g_max 65000.0
/// Surface Elevation (m)
#define hard_h_min (-400.0)
#define hard_h_max 8840.0
#define hard_h_min_km (hard_h_min/1000.0)
#define hard_h_max_km (hard_h_max/1000.0)
/// Ozone
#define hard_o_min 0.0
#define hard_o_max 1000.0
/// Pressure
#define hard_p_min 400.0
#define hard_p_max 1100.0
/// Temperature (Coldest recorded temperature was ~180.)
#define hard_t_min 150.0
#define hard_t_max 350.0
/// Water Vapor
#define hard_v_min 0.0001
#define hard_v_max 599999.0
//#define hard_v_max 600000.0

/// Zenith: 84.0 is v9 max limit; v8 and before = 75.0
#define hard_z_min 0.0
#define hard_z_max 84.0

/**
 * @brief Compute the Correlation Coefficient for two vectors of real values.
 * 
 * corrcoef = Covariance(v1, v2)/(stddev(v1)*stddev(v2))
 * @param   v1  a double[n] array
 * @param   v2  another double[n] array
 * @param   n   the size of v1 and v2
 * @return      Correlation Coefficient
 */
double corrcoef(double *v1, double *v2, int n);

/**
 * @brief Compute the Covariance for two vectors of real values.
 *
 * The covariance is the mean of the products of the deviations of each pair 
 * of values in v1 and v2. 
 * @param   v1  a double[n] array
 * @param   v2  another double[n] array
 * @param   n   the size of v1 and v2
 * @return      Covariance
 */
double covvar(double *v1, double *v2, int n);

/**
 * @brief Uses epsilon to determine if two double values are equal.
 *
 * @param   v1 	first floating point value
 * @param   v2 	second floating point value
 * @return  true if abd difference is < epsilon
 */
bool dequal(double v1, double v2);

/**
 * @brief Recursive function to find v in d[] between imin and imax
 * where the values in d are monotonically decreasing.
 *
 * @param   d    Vector of monotonically decreasing values.
 * @param   v    Target value
 * @param   imin Starting index of search
 * @param   imax Ending index for search
 * @return  Index of closest value >= v.
 */
int dec_binsearch(double *d, double v, int imin, int imax);

/**
 * @brief Recursive function to find v in d[] between imin and imax
 * where the values in d are monotonically increasing.
 *
 * @param   d    Vector of monotonically increasing values.
 * @param   v    Target value
 * @param   imin Starting index of search
 * @param   imax Ending index for search
 * @return  Index of closest value <= v.
 */
int inc_binsearch(double *d, double v, int imin, int imax);

/**
 * @brief Return the distance in meters between two points.
 *
 * @param   lat1    Latitude of first point
 * @param   lon1    Longitude of first point
 * @param   lat2    Latitude of second point
 * @param   lon2    Longitude of second point
 * @return  distance between points 1 and 2 in meters
 */
double distance(double lat1, double lon1, double lat2, double lon2);

/**
 * @brief Get the nearest index to out[ix,iy] when out[] is smaller
 *
 * When resizing from a dataset that is nx_in by ny_in to 
 * a smaller dataset of the size nx_out by ny_out, finds the index
 * in the "in" system that is closest to (ix,iy) in the 
 * "out" system.
 *
 * @param   nx_in   pixels per line in original size
 * @param   ny_in   lines in original size
 * @param   nx_out  pixels per line in new size
 * @param   ny_out  lines in new size
 * @param   ix      pixel in new size
 * @param   iy      line in new size
 * @return          index in original size
 */
int getnearest(int nx_in, int ny_in, int nx_out, int ny_out, int ix, int iy);

/**
 * @brief Get the nearest index to out[ix,iy] when in[] is smaller
 *
 * When resizing from a dataset that is nx_in by ny_in to 
 * a larger dataset of the size nx_out by ny_out, finds the index
 * in the "in" system that is closest to (ix,iy) in the 
 * "out" system.
 *
 * @param   nx_in   pixels per line in original size
 * @param   ny_in   lines in original size
 * @param   nx_out  pixels per line in new size
 * @param   ny_out  lines in new size
 * @param   ix      pixel in new size
 * @param   iy      line in new size
 * @return          index in original size
 */
int getnearestInv(int nx_in, int ny_in, int nx_out, int ny_out, int ix, int iy);

/**
 * @brief Count NaNs
 *
 * Counts the number of NaNs contained in a dataset.
 * @param  in    Pointer to a dataset of double values.
 * @param  nx    Number of pixels per row
 * @param  ny    Number of rows
 * @return count of values considered "Not a Number"
 */
int get_nnans(double *in, int nx, int ny);

/**
 * @brief Test for any NaN in data
 *
 * Returns true when the first NaN is found in a dataset.
 * This is useful when the only information needed is if there is
 * at least one NaN. 
 * @param  in    Pointer to a dataset of double values.
 * @param  nx    Number of pixels per row
 * @param  ny    Number of rows
 * @return true if there is at least one "Not a Number" in the dataset.
 */
bool has_nan(double *in, int nx, int ny);

// inc_binsearch -- see dec_binsearch above

/**
 * @brief Determine if a number should be treated as not-a-number
 *
 * @param   val     Value to be tested for nan
 * @return  true    Value is not a number
 * @note If NEG_IS_NAN is defined, anhy value less than zero is considered 
 * as not a number.
 */
bool is_nan(double val);

/**
 * @brief Resizes a dataset.
 *
 * @param      in      Array containing uint8 values to be resized to a larger array.
 * @param      nx_in   Pixels per line in the in array
 * @param      ny_in   Number of lines in the input array
 * @param[out] out     Array to contain the output of the resize
 * @param      nx_out  Pixels per line in the out array
 * @param      ny_out  Number of lines in the out array
 */
void imresize_u8(uint8* in, int nx_in, int ny_in, uint8* out, int nx_out, int ny_out);
void imresize(double *in, int nx_in, int ny_in, double *out, int nx_out, int ny_out);

/**
 * @brief Linear interpolate the y value for given x between two of 
 * three x,y points.
 *
 * If xx < x[1], interpolates between x[0],y[0] and x[1],y[1].
 * Otherwise interpolates between x[1],y[1] and x[2],y[2].
 * @param   x   double[3] set of x coordinates
 * @param   y   double[3] set of y coordinates
 * @param   xx  given x value
 * @return      interpolated y value
 */
double interp1d(double *x, double *y, double xx);

/**
 * @brief Linear interpolate the y value for given x between two of
 * npts x,y points.
 *
 * Similar to interp1d, but any number of points may be given. There 
 * must at least two points. 
 * @param   x   double[npts] set of x coordinates
 * @param   y   double[npts] set of y coordinates
 * @param   xx  given x value
 * @param   npts number of points
 * @return      interpolated y value
 * @depricated use interp1 from interps.h
 */
double interp1d_npts(double *x, double *y, double xx, int npts);

/**
 * @brief Replace NaN values using Inverse Distance Weighting.
 *
 * @param[in,out]  in   nx by ny array
 * @param          nx   size of the row
 * @param          ny   numer of rows
 * @param          r2   radius (pixels) of maximum distance e.g., 5.0
 * @param          e    distance exponent, e.g., -2.0
 */
void interp_id(double *in, int nx, int ny, double r2, double e);

/**
 * @brief Replace NaN values using Inverse Distance Weighting where e = -2.0.
 *
 * This version is like interp_id, but optimizes performance for e == -2.0.
 * Optimizations include not calling sqrt and not calling pow in a
 * loop O(nNan * nGood).
 * @param[in,out]  in   nx by ny array
 * @param          nx   size of the row
 * @param          ny   numer of rows
 * @param          r2   radius (pixels) of maximum distance e.g., 5.0
 */
void interp_id_eneg2(double *in, int nx, int ny, double r2);

/**
 * @brief Perform interp_id_eneg2 on multiple datasets.
 *
 * Use this when multiple datasets need to be interpolated, and each dataset
 * contains nans in the same locations.
 *
 * The first dataset is the master, used to determine what points are to be
 * included in each interpolation. Once determined, the same indexes are applied 
 * to interpolate the additional datasets. 
 *
 * @param[in,out]  ds   list of datasets
 * @param          nds  number of datasets
 * @param          nx   size of the row
 * @param          ny   number of rows
 * @param          r2   radius (pixels) of maximum distance e.g., 5.0
 */
void multiset_interp_id(double *ds[], int nds, int nx, int ny, double r2);

/**
 * @brief Compute Mean
 * 
 * Sums n values from v and divides by n.
 * @param   v   Vector of values
 * @param   n   Number of values in v
 * @return      sum(v[0..n])/n
 */
double mean(double *v, int n);

/**
 * @brief Find the minimum and maximum values in a 2D dataset.
 *
 * @param[in]   d    dataset
 * @param[in]   nx   number of pixels per lin
 * @param[in]   ny   number of lines
 * @param[in]   fill fill value to be ignored
 * @param[out]  min  pointer to minimum value
 * @param[out]  max  pointer to maximum value
 */
 void minmax2d(double *d, int nx, int ny, double fill, double *min, double *max);

/**
 * @brief Rounds x to closest multiple of y. Halfway between
 * is rounded to the value farthest from zero. (see round).
 *
 * @param   x   value to be rounded
 * @param   y   granularity 
 * @return      closes multiple of Y to x
 */
double round2(double x, double y);

/**
 * @brief Smooth the values in a 2D array
 *
 * Performs smoothing over values within +/- dx,dy
 * For each in[i,j] uses averages of the values for the 
 * row within +/- dx and the column within +/- dy to smooth
 * the values in the in array.
 * @param[in,out]   in  nx by ny array of values
 * @param           nx  number of pixels per line
 * @param           ny  number of lines
 * @param           dx  smooth radius in number of pixels
 * @param           dy  smooth radius in number of lines
 */
void smooth2a(double *in, int nx, int ny, int dx, int dy);
void smooth2a_i16(int16 *in, int nx, int ny, int dx, int dy);

/**
 * @brief Compute Standard Deviation
 *
 * Computes the standard deviation from a vector of n values.
 * @param   v   vector of values
 * @param   n   number of values in v
 * @return      standard deviation of the values in v
 */
double stddev(double *v, int n);
