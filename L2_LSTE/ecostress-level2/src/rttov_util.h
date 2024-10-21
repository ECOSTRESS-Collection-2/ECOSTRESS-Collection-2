#pragma once
/// @file Utility functions for RTTOV processing.
//
//  References: 
/// @author Simon Latyshev, Raytheon, [SL]
/// @author Robert Freepartner, JPL, [BF]
//
/// @copyright (c) 2016 JPL, All rights reserved

#include "matrix.h"

/// @brief Copy data from matrix to double[].
void mat2d_copy_mat(double *mat, Matrix *mat_src, int nrows, int ncols);

/// @brief Copy data from double[] to double[].
void mat2d_copy_rttov(double *mat, double *mat_src, int nrows, int ncols);

/// @brief Return the mean value of values that are greater than the given value.
double mat2d_get_mean_more_val_rttov(double *mat, double val, int nrows, int ncols);

/// @brief Return the mean value of values that are greater than zero.
double mat2d_get_mean_more_zero_rttov(double *mat, int nrows, int ncols);

/// @brief Multiplies each non-nan value in the matrix by the given value.
void mat2d_scale_rttov(double *mat, double val, int nrows, int ncols);

/// @brief Replaces entries containing eps_val with val.
void mat2d_set_eps_to_val_rttov(
        double *mat, double eps_val, double val, int nrows, int ncols);

/// @brief Sets column [j] to the given value.
void mat2d_set_j_column_to_val_rttov(double *mat, double val, int j, int nrows, int ncols);

/// @brief Sets entries less than minval to val.
void mat2d_set_less_minval_to_val_rttov(double *mat, double minval, double val, int nrows, int ncols);

/// @brief Sets entries less than zero to val.
void mat2d_set_less_zero_to_val_rttov(double *mat, double val, int nrows, int ncols);

/// @brief Sets entries containing zero to val.
void mat2d_set_zero_to_val_rttov(double *mat, double val, int nrows, int ncols);

/// @brief transpose a k x j matrix to a j x k matrix in place.
void mat2d_transpose_rttov(double *mat, int kmax, int jmax);

/// @brief Deallocates 3d array.
void mat3d_clear_rttov(double *mat[], int imax);

/// @brief mat[i][k][j] = vect[i] for 0 <= i < imax.
void mat3d_copy_from_vector_kj(double vect[], int k, int j, Mat3d *mat);

/// @brief Copy data from a Mat3d matrix to a double*[] 3d representation.
void mat3d_copy_mat(double *mat[], Mat3d *mat_src, int imax, int kmax, int jmax);

/// @brief vect[i] = mat[i][k][j] for 0 <= i < imax.
void mat3d_copy_to_vector_kj_rttov(double vect[], int k, int j, Mat3d *mat, int imax, int kmax, int jmax);

/// @brief Returns the mean of values that are greater than the given val.
double mat3d_get_mean_more_eq_val_rttov(double *mat[], double val, int imax, int kmax, int jmax);

/// @brief Use calloc to allocate kmax x jmax doubles for each of imax entries in mat.
void mat3d_init_zeros_rttov(double *mat[], int imax, int kmax, int jmax);

/// @brief Multiplies each non-nan value in the matrix by the given value.
void mat3d_scale_rttov(double *mat[], double val, int imax, int kmax, int jmax);

/// @brief Set mat[i][k][j] = val for all i levels of a 3D matrix. 
void mat3d_set_kj_to_val(Mat3d *mat, double val, int k, int j);

/// @brief Replavce zeros with the given value.
void mat3d_set_zero_to_val(Mat3d *mat, double val);

// @brief Transpose (out-of-place), convert Mat[i][k][j] into Mat[i][j][k];
void mat3d_transpose_rttov(double *mat, double *mat_src[], int imax, int kmax, int jmax);

/// @brief Finds the value > -.- nearest the end of the vector.
double vect_get_last_positive_rttov(double *vect, int imax);

/// @brief Returns true if no value is >= 0.0.
bool vect_is_all_negative_rttov(double *vect, int imax);

/// @brief Returns true is there is any non-nan value less than zero.
bool vect_is_negative_found_rttov(double *vect, int imax);

/// @brief Sets non-nan values to fmod with val. 
/// Note: This is used to convert Longitudes from -180.0 to 180.0 to 0.0 to 360.0.
void vect_mod_rttov(double *vect, double val, int imax);

/// @brief Sets vect[i] = vect[i-1] when vect[i] < 0.0.
void vect_reset_negative_rttov(double *vect, int imax);

/// @brief Multiplies each non-nan value in the vector by the given value.
void vect_scale_rttov(double *vect, double val, int imax);
