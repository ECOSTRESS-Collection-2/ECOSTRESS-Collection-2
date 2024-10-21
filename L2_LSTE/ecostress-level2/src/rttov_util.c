// See rttov_util.h for credits and documentation.
// Copyright (c) 2016, JPL, All rights reserved

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "lste_lib.h"
#include "rttov_util.h"

void mat2d_copy_mat(double *mat, Matrix *mat_src, int nrows, int ncols)
{
    assert(mat_src != NULL);
    mat2d_copy_rttov(mat, mat_src->vals, nrows, ncols);
}

void mat2d_copy_rttov(double *mat, double *mat_src, int nrows, int ncols)
{
    assert(mat != NULL);
    assert(mat_src != NULL);
    size_t nbytes = nrows * ncols * sizeof(double);
    memcpy(mat, mat_src, nbytes);
}

// Mean(mat(mat > val)), get mean (average) of selected elements
double mat2d_get_mean_more_val_rttov(
        double *mat, double val, int nrows, int ncols)
{
    int k = 0;
    int j = 0;
    double dres = 0.0;	// default result
    double dsum = 0.0;	// sum of selected elements
    int count = 0;		// number of selected elements

    assert (mat != NULL);

    for(k = 0; k < nrows; k++) {
        for(j = 0; j < ncols; j++) {
            if(mat[k*ncols+j] > val) 
            {
                dsum += mat[k*ncols+j];
                count++;
            }
        }
    }

    // Compute arithmetic average (mean)
    if (count > 0) 
        dres = dsum /((double) count);	

    return dres; 
}

// Mean(mat(mat > 0)), get mean (average) of selected elements
double mat2d_get_mean_more_zero_rttov(double *mat, int nrows, int ncols)
{
    return mat2d_get_mean_more_val_rttov(mat, 0.0, nrows, ncols);
}

// Mat *= val, scale by scalar value
void mat2d_scale_rttov(double *mat, double val, int nrows, int ncols)
{
    int k = 0;
    int j = 0;

    assert(mat != NULL);
    for(k = 0; k < nrows; k++) {
        for(j = 0; j < ncols; j++) {
            if (!is_nan(mat[k*ncols+j]))
                mat[k*ncols+j] *= val;
        }
    }
}

// Mat(mat == eps_val) = val, replace zero values with val
void mat2d_set_eps_to_val_rttov(double *mat, double eps_val, double val, int nrows, int ncols)
{
    int k = 0;
    int j = 0;

    assert(mat != NULL);
    for(k = 0; k < nrows; k++) {
        for(j = 0; j < ncols; j++) {
            if(!is_nan(mat[k*ncols+j]) && fabs(mat[k*ncols+j] - eps_val) < epsilon)
                mat[k*ncols+j] = val;
        }
    }
}

// Mat[][j] = val, sel j-th column elements to val
void mat2d_set_j_column_to_val_rttov(double *mat, double val, int j, int nrows, int ncols)
{
    int k = 0;

    assert(mat != NULL);

    for(k = 0; k < nrows; k++) 
    {
        mat[k*ncols+j] = val;
    }
}

void mat2d_set_less_minval_to_val_rttov(double *mat, double minval, double val, int nrows, int ncols)
{
    int k = 0;
    int j = 0;
    assert(mat != NULL);
    for(k = 0; k < nrows; k++) {
        for(j = 0; j < ncols; j++) {
            if(!is_nan(mat[k*ncols+j]) && mat[k*ncols+j] < minval) 
                mat[k*ncols+j] = val;
        }
    }
}

void mat2d_set_less_zero_to_val_rttov(double *mat, double val, int nrows, int ncols)
{
    mat2d_set_less_minval_to_val_rttov(mat, 0.0, val, nrows, ncols);
}

void mat2d_set_zero_to_val_rttov(double *mat, double val, int nrows, int ncols)
{
    mat2d_set_eps_to_val_rttov(mat, 0.0, val, nrows, ncols);
}

// Mat = T(mat), IN-PLACE Transpose mat(NxM) to mat(MxN)
// 1. Original matrix dimension sizes (kmax, jmax)
// 2. Resulting dimension sizes (jmax,  kmax)

void mat2d_transpose_rttov(double *mat, int kmax, int jmax)
{
    int k = 0;
    int j = 0;
    double *temp = NULL;	// temp array

    assert(mat != NULL);

    // Allocate temp matrix
    temp = (double *) calloc(kmax * jmax, sizeof(double));
    assert(temp != NULL);

    // Copy matrix values to temp

    mat2d_copy_rttov(temp, mat, kmax, jmax); 

    // Copy element values in transposed order

    for(k = 0; k < kmax; k++) {
        for(j = 0; j < jmax; j++) {
            mat[j*kmax+k] = temp[k*jmax+j];
        }
    }

    free(temp);
}

// Clear, deallocate, 2d arrays for 3d matrix
// 1. Assuming 3d matrix is defined as double *mat[N], N == imax

void mat3d_clear_rttov(double *mat[], int imax)
{
    int i = 0;
    assert(mat != NULL);
    for(i = 0; i < imax; i++) 
    {
        if(mat[i] != NULL) free(mat[i]);
    }
}

// Mat[][k][j] = Vector[], copy values from vector

void mat3d_copy_from_vector_kj(double vect[], int k, int j, Mat3d *mat)
{
    int i;
    for(i = 0; i < mat->size1; i++) 
    {
        mat3d_set(mat, i, k, j, vect[i]);
    }
}

// Mat = mat_src, copy values
// 1. Assuming that both matrices allocated already
// 2. Assuming that mat[i][j][k] sizes match
// 3. UPDATE: use memcopy operation instead of for loops

void mat3d_copy_mat(double *mat[], Mat3d *mat_src, int imax, int kmax, int jmax)
{
    int i = 0;
    int k = 0;
    int j = 0;

    assert(mat != NULL);
    assert(mat_src != NULL);
    assert(mat_src->vals != NULL);

    for(i = 0; i < imax; i++) {
        for(k = 0; k < kmax; k++) {
            for(j = 0; j < jmax; j++) {
                mat[i][k*jmax+j] = mat3d_get(mat_src, i, k, j);
            }
        }
    }
}

// Vector[] = mat[][k][j], copy values, getting matrix slice as a vector
// 1. Assuming vector is defined as double vect[N]
// 2. Assuming matching mat[N][kmax][jmax] size, where N == imax

void mat3d_copy_to_vector_kj_rttov(double vect[], int k, int j, Mat3d *mat, int imax, int kmax, int jmax)
{
    int i = 0;
    assert(vect != NULL);
    assert(mat != NULL);
    assert(mat->vals != NULL);
    assert(k >= 0 && k < kmax);
    assert(j >= 0 && j < jmax);

    for(i = 0; i < imax; i++) 
    {
        vect[i] = mat3d_get(mat, i, k, j);
    }
}

// Mean(mat(mat > val)), get mean (average) of selected elements

double mat3d_get_mean_more_eq_val_rttov(double *mat[], double val, int imax, int kmax, int jmax)
{
    int i = 0;
    int k = 0;
    int j = 0;
    double dres = 0.0;	// default result
    double dsum = 0.0;	// sum of selected elements
    int count = 0;		// number of selected elements

    assert(mat != NULL);

    for(i = 0; i < imax; i++) {
        for(k = 0; k < kmax; k++) {
            for(j = 0; j < jmax; j++) {
                if(!is_nan(mat[i][k*jmax+j]) && mat[i][k*jmax+j] > val)
                {
                    dsum += mat[i][k*jmax+j];
                    count++;
                }
            }
        }
    }

    // Compute arithmetic average (mean)
    if(count > 0) 
        dres = dsum /((double) count);	

    return dres;
}

// Initialize, allocate, 2d arrays for 3d matrix
// 1. Assuming 3d matrix is defined as double *mat[N], N == imax
// 2. Setting values to zeros using calloc()

void mat3d_init_zeros_rttov(double *mat[], int imax, int kmax, int jmax)
{
    int i = 0;
    assert(mat != NULL);
    for(i = 0; i < imax; i++) 
    {
        mat[i] = (double *) calloc(kmax * jmax, sizeof(double));
        assert(mat[i] != NULL);
    }
}

// Mat *= val, scale by scalar value
// UPDATE: combine k and j loop, using nmax = kmax * jmax

void mat3d_scale_rttov(double *mat[], double val, int imax, int kmax, int jmax)
{
    int i = 0;
    int k = 0;
    int j = 0;

    assert(mat != NULL);

    for(i = 0; i < imax; i++) {
        for(k = 0; k < kmax; k++) {
            for(j = 0; j < jmax; j++) {
                if (!is_nan(mat[i][k*jmax+j]))
                    mat[i][k*jmax+j] *= val;
            }
        }
    }
}

// Mat[][k][j] = val, set matrix slice values to val

void mat3d_set_kj_to_val(Mat3d *mat, double val, int k, int j)
{
    int i = 0;
    for(i = 0; i < mat->size1; i++) 
    {
        mat3d_set(mat, i, k, j, val);
    }
}

// Mat(mat == 0) = val, replace zero values with val
void mat3d_set_zero_to_val(Mat3d *mat, double val)
{
    int i;
    for (i = 0; i < mat->size1 * mat->size2 * mat->size3; i++)
    {
        if (fabs(mat->vals[i]) < epsilon)
            mat->vals[i] = val;
    }
}

// Transpose (out-of-place), convert Mat[i][k][j] into Mat[i][j][k];

void mat3d_transpose_rttov(double *mat, double *mat_src[], int imax, int kmax, int jmax)
{
    int i = 0;
    int k = 0;
    int j = 0;

    assert(mat != NULL);
    assert(mat_src != NULL);

    for(i = 0; i < imax; i++) {
        for(k = 0; k < kmax; k++) {
            for(j = 0; j < jmax; j++) {
                mat[i*jmax*kmax + j*kmax+k] = mat_src[i][k*jmax+j];
            }
        }
    }
}

// Get, find and return, the last positive value
// 1. Returning 0.0 if not found // UPDATE: issue an error in this case

double vect_get_last_positive_rttov(double *vect, int imax)
{
    int i = 0;
    double val = 0.0; // default value to return

    assert(vect != NULL);

    for(i = imax-1; i >= 0; i--) // iterating in reverse order
    {
        if(vect[i] > 0.0) return vect[i];
    }
    return val;
}

// Check if all values are negative
// 1. vect[i] < 0 for all i=[0,imax]
// 2. Assuming vector is defined as vect[N], where N == imax

bool vect_is_all_negative_rttov(double *vect, int imax)
{
    int i = 0;

    assert(vect != NULL);
    for(i = 0; i < imax; i++) 
    {
        if(vect[i] >= 0.0) return false;
    }
    return true;
}

// Check if at least one negative value is found
// 1. vect[i] < 0 for least one element

bool vect_is_negative_found_rttov(double *vect, int imax)
{
    int i = 0;

    assert(vect != NULL);

    for(i = 0; i < imax; i++) 
    {
        if(!is_nan(vect[i]) && vect[i] < 0.0) return true;
    }
    return false;
}

// Mod(vector, val), reset array elements, to mod(x, val)
// 1. CHECK: Mod is defined as x = val - x*floor(x/val)
// 2. Don't update elements for val == 0.0
// [BF] The intention here is to convert Longitude in the
// range = -180.0 .. 180.0 to 0 .. 360.0
void vect_mod_rttov(double *vect, double val, int imax)
{
    int i = 0;
    assert(vect != NULL);
    for(i = 0; i < imax; i++) 
    {
        // [BF] This was erroneous: vect[i] = val - vect[i] * floor(vect[i] / val);
        if (!is_nan(vect[i]))
            vect[i] = fmod(vect[i] + val, val);
    }
}

// Reset negative values vect[i] = vect[i-1] iteratively
// CHECK: NOTE: if(vect[0] < 0) the value is reset to zero 

void vect_reset_negative_rttov(double *vect, int imax)
{
    int i = 0;
    double val = 0.0;        // the last non-negative value

    assert(vect != NULL);

    for(i = 0; i < imax; i++) 
    {
        if(vect[i] >= 0.0)
        {
            val = vect[i];
        }
        else
        {
            vect[i] = val;
        }
        // NOTE: if(vect[0] < 0) the value is reset to zero
    }
}

// Vect *= val, scale by scalar value

void vect_scale_rttov(double *vect, double val, int imax)
{
    int i = 0;

    assert(vect != NULL);

    for(i = 0; i < imax; i++) 
    {
        if (!is_nan(vect[i]))
            vect[i] *= val;
    }
}
