// matrix.c
// See matrix.h for credits and documentation.
//
// Copyright (c) JPL, All rights reserved

#include "matrix.h"
#include <assert.h>
#include <string.h>
#include "lste_lib.h"

// --------------------------------------------------------------
// Vector 1d float32: Initialize, allocate, setting values to zeros
// --------------------------------------------------------------

void vec_float32_init(VecFloat32 *vec)
{
    assert(vec != NULL);
    vec->size = 0;
    vec->vals = NULL;
}

void vec_float32_alloc(VecFloat32 *vec, size_t size)
{
    assert(vec != NULL);
    vec_float32_clear(vec);
    vec->size = size;
    vec->vals = (float32 *) calloc(size, sizeof(float32));
    assert(vec->vals != NULL);
}

float32 vec_float32_get(VecFloat32 *vec, size_t x)
{
    return vec->vals[x];
}

// Safe check version
float32 vec_float32_gets(VecFloat32 *vec, size_t x)
{
    if (x >= vec->size) return MATRIX_NAN;
    return vec->vals[x];
}

void vec_float32_set(VecFloat32 *vec, size_t x, float32 val)
{
    vec->vals[x] = val;
}

void vec_float32_sets(VecFloat32 *vec, size_t x, float32 val)
{
    if (x < vec->size)
        vec->vals[x] = val;
}
// --------------------------------------------------------------
// Vector 1d float32: Clear vector, deallocating elements array
// --------------------------------------------------------------

void vec_float32_clear(VecFloat32 *vec)
{
    assert(vec != NULL);
    if(vec->vals != NULL) 
    {
        free(vec->vals);    // deallocate
        vec->vals = NULL;    // reset to NULL
    }
    vec->size = 0;    // reset to zero 
}

// --------------------------------------------------------------
// Vector 1d double: Initialize, allocate, setting values to zeros
// --------------------------------------------------------------

void vec_init(Vector *vec)
{
    assert(vec != NULL);
    vec->size = 0;
    vec->vals = NULL;
}

void vec_alloc(Vector *vec, size_t size)
{
    assert(vec != NULL);
    vec_clear(vec);
    vec->size = size;
    vec->vals = (double *) calloc(size, sizeof(double));
    assert(vec->vals != NULL);
}

double vec_get(Vector *vec, size_t x)
{
    return vec->vals[x];
}

double vec_gets(Vector *vec, size_t x)
{
    if (x >= vec->size) return MATRIX_NAN;
    return vec->vals[x];
}

void vec_set(Vector *vec, size_t x, double val)
{
    vec->vals[x] = val;
}

void vec_sets(Vector *vec, size_t x, double val)
{
    if (x < vec->size)
        vec->vals[x] = val;
}

void vec_copy(Vector *dest, Vector *source)
{
    vec_alloc(dest, source->size);
    memcpy(dest->vals, source->vals, source->size*sizeof(double));
}
// --------------------------------------------------------------
// Vector 1d double: Clear vector, deallocating elements array
// --------------------------------------------------------------

void vec_clear(Vector *vec)
{
    assert(vec != NULL);
    if(vec->vals != NULL) 
    {
        free(vec->vals);    // deallocate array
        vec->vals = NULL;    // reset to NULL
    }
    vec->size = 0;    // reset to zero 
}

// --------------------------------------------------------------
// Vector 1d double: Copy values converting from float32 to double
// --------------------------------------------------------------

void vec_copy_float32(Vector *vec, const VecFloat32 *vec_src)
{
    assert(vec != NULL);
    assert(vec_src != NULL);

    // Copy values converting element types
    size_t i;
    for(i = 0; i < vec->size; i++)
    {
        vec->vals[i] = (double)vec_src->vals[i];
    }
}

// --------------------------------------------------------------
// Matrix 2d uint8: Initialize, allocate, setting values to zeros
// --------------------------------------------------------------

void mat_uint8_init(MatUint8 *mat)
{
    assert(mat != NULL);
    mat->size1 = 0;
    mat->size2 = 0;
    mat->vals = NULL;

}

void mat_uint8_alloc(MatUint8 *mat, size_t size1, size_t size2)
{
    assert(mat !=  NULL);
    mat_uint8_clear(mat);
    mat->size1 = size1;
    mat->size2 = size2;
    mat->vals = (uint8 *) calloc(size1 * size2, sizeof(uint8));
    assert(mat->vals != NULL);
}

// --------------------------------------------------------------
// Matrix 2d uint8: Clear matrix, deallocating elements array
// --------------------------------------------------------------

void mat_uint8_clear(MatUint8 *mat)
{
    assert(mat != NULL);
    if(mat->vals != NULL) 
    {
        free(mat->vals);    // deallocate
        mat->vals = NULL;    // reset to NULL
    }
    mat->size1 = 0;    // reset to zero 
    mat->size2 = 0;    // reset to zero
}

void mat_uint8_fill(MatUint8 *mat, uint8 value)
{
    int i;
    for (i = 0; i < mat->size1 * mat->size2; i++)
        mat->vals[i] = value;
}

uint8 mat_uint8_get(MatUint8 *mat, size_t y, size_t x)
{
    return mat->vals[y * mat->size2 + x];
}

uint8 mat_uint8_gets(MatUint8 *mat, size_t y, size_t x)
{
    if (x >= mat->size2) return MATRIX_UINT8_NAN;
    if (y >= mat->size1) return MATRIX_UINT8_NAN;
    return mat->vals[y * mat->size2 + x];
}

void mat_uint8_set(MatUint8 *mat, size_t y, size_t x, uint8 val)
{
    mat->vals[y * mat->size2 + x] = val;
}

void mat_uint8_sets(MatUint8 *mat, size_t y, size_t x, uint8 val)
{
    if (y < mat->size1 && x < mat->size2)
        mat->vals[y * mat->size2 + x] = val;
}

// --------------------------------------------------------------
// Matrix 2d int16: Initialize, allocate, setting values to zeros
// --------------------------------------------------------------

void mat_int16_init(MatInt16 *mat)
{
    assert(mat != NULL);
    mat->size1 = 0;
    mat->size2 = 0;
    mat->vals = NULL;
}

void mat_int16_alloc(MatInt16 *mat, size_t size1, size_t size2)
{
    assert(mat != NULL);
    mat_int16_clear(mat);
    mat->size1 = size1;
    mat->size2 = size2;
    mat->vals = (int16 *) calloc(size1 * size2, sizeof(int16));
    assert(mat->vals != NULL);
}

// --------------------------------------------------------------
// Matrix 2d int16: Clear matrix, deallocating elements array
// --------------------------------------------------------------

void mat_int16_clear(MatInt16 *mat)
{
    assert(mat != NULL);
    if(mat->vals != NULL) 
    {
        free(mat->vals);    // deallocate
        mat->vals = NULL;    // reset to NULL
    }
    mat->size1 = 0;    // reset to zero 
    mat->size2 = 0;    // reset to zero
}

int16 mat_int16_get(MatInt16 *mat, size_t y, size_t x)
{
    return mat->vals[y * mat->size2 + x];
}

int16 mat_int16_gets(MatInt16 *mat, size_t y, size_t x)
{
    if (x >= mat->size2) return MATRIX_NAN;
    if (y >= mat->size1) return MATRIX_NAN;
    return mat->vals[y * mat->size2 + x];
}

void mat_int16_set(MatInt16 *mat, size_t y, size_t x, int16 val)
{
    mat->vals[y * mat->size2 + x] = val;
}

void mat_int16_sets(MatInt16 *mat, size_t y, size_t x, int16 val)
{
    if (x < mat->size2 && y < mat->size1)
        mat->vals[y * mat->size2 + x] = val;
}

// --------------------------------------------------------------
// Matrix 2d int32: Initialize, allocate, setting values to zeros
// --------------------------------------------------------------

void mat_int32_init(MatInt32 *mat)
{
    assert(mat != NULL);
    mat->size1 = 0;
    mat->size2 = 0;
    mat->vals = NULL;
}

void mat_int32_alloc(MatInt32 *mat, size_t size1, size_t size2)
{
    assert(mat != NULL);
    mat_int32_clear(mat);
    mat->size1 = size1;
    mat->size2 = size2;
    mat->vals = (int32 *) calloc(size1 * size2, sizeof(int32));
    assert(mat->vals != NULL);
}

// --------------------------------------------------------------
// Matrix 2d int32: Clear matrix, deallocating elements array
// --------------------------------------------------------------

void mat_int32_clear(MatInt32 *mat)
{
    assert(mat != NULL);
    if(mat->vals != NULL) 
    {
        free(mat->vals);    // deallocate
        mat->vals = NULL;    // reset to NULL
    }
    mat->size1 = 0;    // reset to zero 
    mat->size2 = 0;    // reset to zero
}

int32 mat_int32_get(MatInt32 *mat, size_t y, size_t x)
{
    return mat->vals[y * mat->size2 + x];
}

int32 mat_int32_gets(MatInt32 *mat, size_t y, size_t x)
{
    if (x >= mat->size2) return MATRIX_NAN;
    if (y >= mat->size1) return MATRIX_NAN;
    return mat->vals[y * mat->size2 + x];
}

void mat_int32_set(MatInt32 *mat, size_t y, size_t x, int32 val)
{
    mat->vals[y * mat->size2 + x] = val;
}

void mat_int32_sets(MatInt32 *mat, size_t y, size_t x, int32 val)
{
    if (x < mat->size2 && y < mat->size1)
        mat->vals[y * mat->size2 + x] = val;
}

// --------------------------------------------------------------
// Unsigned Int16
// --------------------------------------------------------------

void mat_uint16_init(MatUint16 *mat)
{
    assert(mat != NULL);
    mat->size1 = 0;
    mat->size2 = 0;
    mat->vals = NULL;
}

void mat_uint16_alloc(MatUint16 *mat, size_t size1, size_t size2)
{
    assert(mat != NULL);
    mat_uint16_clear(mat);
    mat->size1 = size1;
    mat->size2 = size2;
    mat->vals = (uint16 *) calloc(size1 * size2, sizeof(uint16));
    assert(mat->vals != NULL);
}

void mat_uint16_clear(MatUint16 *mat)
{
    assert(mat != NULL);
    if(mat->vals != NULL) 
    {
        free(mat->vals);    // deallocate
        mat->vals = NULL;    // reset to NULL
    }
    mat->size1 = 0;    // reset to zero 
    mat->size2 = 0;    // reset to zero
}

uint16 mat_uint16_get(MatUint16 *mat, size_t y, size_t x)
{
    return mat->vals[y * mat->size2 + x];
}

uint16 mat_uint16_gets(MatUint16 *mat, size_t y, size_t x)
{
    if (x >= mat->size2) return MATRIX_NAN;
    if (y >= mat->size1) return MATRIX_NAN;
    return mat->vals[y * mat->size2 + x];
}

void mat_uint16_set(MatUint16 *mat, size_t y, size_t x, uint16 val)
{
    mat->vals[y * mat->size2 + x] = val;
}

void mat_uint16_sets(MatUint16 *mat, size_t y, size_t x, uint16 val)
{
    if (x < mat->size2 && y < mat->size1)
        mat->vals[y * mat->size2 + x] = val;
}

// --------------------------------------------------------------
// Matrix 2d float32: Initialize, allocate, setting values to zeros
// --------------------------------------------------------------

void mat_float32_init(MatFloat32 *mat)
{
    assert(mat != NULL);
    mat->size1 = 0;
    mat->size2 = 0;
    mat->vals = NULL;
}

void mat_float32_alloc(MatFloat32 *mat, size_t size1, size_t size2)
{
    assert(mat != NULL);
    mat_float32_clear(mat);
    mat->size1 = size1;
    mat->size2 = size2;
    mat->vals = (float32 *) calloc(size1 * size2, sizeof(float32));
    assert(mat != NULL);
}

// --------------------------------------------------------------
// Matrix 2d float32: Clear matrix, deallocating elements array
// --------------------------------------------------------------

void mat_float32_clear(MatFloat32 *mat)
{
    assert(mat != NULL);
    if(mat->vals != NULL) 
    {
        free(mat->vals);    // deallocate
        mat->vals = NULL;    // reset to NULL
    }
    mat->size1 = 0;    // reset to zero 
    mat->size2 = 0;    // reset to zero
}

float32 mat_float32_get(MatFloat32 *mat, size_t y, size_t x)
{
    return mat->vals[y * mat->size2 + x];
}

float32 mat_float32_gets(MatFloat32 *mat, size_t y, size_t x)
{
    if (x >= mat->size2) return MATRIX_NAN;
    if (y >= mat->size1) return MATRIX_NAN;
    return mat->vals[y * mat->size2 + x];
}

void mat_float32_set(MatFloat32 *mat, size_t y, size_t x, float32 val)
{
    mat->vals[y * mat->size2 + x] = val;
}

void mat_float32_sets(MatFloat32 *mat, size_t y, size_t x, float32 val)
{
    if (x < mat->size2 && y < mat->size1)
        mat->vals[y * mat->size2 + x] = val;
}

// --------------------------------------------------------------
// Matrix 2d double: Initialize, allocate, setting values to zeros
// --------------------------------------------------------------

void mat_init(Matrix *mat)
{
    assert(mat != NULL);
    mat->size1 = 0;
    mat->size2 = 0;
    mat->vals =  NULL;
}

void mat_alloc(Matrix *mat, size_t size1, size_t size2)
{
    assert(mat != NULL);
    mat_clear(mat);
    mat->size1 = size1;    // set sizes
    mat->size2 = size2;
    mat->vals = (double *) calloc(size1 * size2, sizeof(double)); // allocate elements
    assert(mat->vals != NULL);
}

// --------------------------------------------------------------
// Matrix 2d double: Clear matrix, deallocating elements array
// --------------------------------------------------------------

void mat_clear(Matrix *mat)
{
    assert(mat != NULL);
    if(mat->vals != NULL) 
    {
        free(mat->vals);    // deallocate array
        mat->vals = NULL;    // reset to NULL
    }
    mat->size1 = 0;    // reset to zero 
    mat->size2 = 0;    // reset to zero
}

double mat_get(Matrix *mat, size_t y, size_t x)
{
    return mat->vals[y * mat->size2 + x];
}

double mat_gets(Matrix *mat, size_t y, size_t x)
{
    if (x >= mat->size2) return MATRIX_NAN;
    if (y >= mat->size1) return MATRIX_NAN;
    return mat->vals[y * mat->size2 + x];
}

void mat_set(Matrix *mat, size_t y, size_t x, double val)
{
    mat->vals[y * mat->size2 + x] = val;
}

void mat_sets(Matrix *mat, size_t y, size_t x, double val)
{
    if (x < mat->size2 && y < mat->size1)
        mat->vals[y * mat->size2 + x] = val;
}

void mat_copy(Matrix *dest, Matrix *source)
{
    assert(source != NULL);
    assert(dest != NULL);
    assert(source->size1 != 0);
    assert(source->size2 != 0);
    mat_alloc(dest, source->size1, source->size2);
    memcpy(dest->vals, source->vals, 
        source->size1 * source->size2 * sizeof(double));
}

void mat_transpose(Matrix *mat_prime, Matrix *mat)
{
    mat_alloc(mat_prime, mat->size2, mat->size1);
    int i, j;
    for (i = 0; i < mat->size1; i++)
    {
        for (j = 0; j < mat->size2; j++)
        {
            mat_set(mat_prime, j, i, mat_get(mat, i, j));
        }
    }
}

// --------------------------------------------------------------
// Matrix 2d double: Copy values converting from int16 to double
// --------------------------------------------------------------

void mat_copy_int16(Matrix *mat, const int16 *mat_src)
{
    assert(mat != NULL);
    assert(mat_src != NULL);
    // Copy values converting element types
    
    size_t i = 0;
    size_t imax = mat->size1 * mat->size2;    // total size
    
    for(i = 0; i < imax; i++)
    {
        mat->vals[i] = (double) mat_src[i];
    }
}

void mat_copy_uint16(Matrix *mat, const uint16 *mat_src)
{
    assert(mat != NULL);
    assert(mat_src != NULL);
    // Copy values converting element types
    
    size_t i = 0;
    size_t imax = mat->size1 * mat->size2;    // total size
    
    for(i = 0; i < imax; i++)
    {
        mat->vals[i] = (double) mat_src[i];
    }
}

void mat_copy_uint_8_to_16(MatUint16 *mat, const MatUint8 *mat_src)
{
    assert(mat != NULL);
    assert(mat_src != NULL);
    // Copy values converting element types
    mat_uint16_alloc(mat, mat_src->size1, mat_src->size2);
    size_t i = 0;
    size_t imax = mat->size1 * mat->size2;    // total size
    
    for(i = 0; i < imax; i++)
    {
        mat->vals[i] = (uint16)mat_src->vals[i];
    }
}

void mat_fill(Matrix *mat, double value)
{
    assert(mat = NULL);
    assert(mat->vals != NULL);
    size_t i;
    for (i = 0; i < mat->size1 * mat->size2; i++)
    {
        mat->vals[i] = value;
    }
}

// --------------------------------------------------------------
// Matrix 2d double: Copy values converting from float32 to double
// --------------------------------------------------------------

void mat_copy_float32(Matrix *mat, const MatFloat32 *mat_src)
{
    assert(mat != NULL);
    assert(mat_src != NULL);
    // Copy values converting element types
    size_t i = 0;
    size_t imax = mat->size1 * mat->size2;    // total size
    
    for(i = 0; i < imax; i++)
    {
        mat->vals[i] = (double) mat_src->vals[i]; // float32 to double value
    }
}

// --------------------------------------------------------------
// Matrix 2d double: Set values corresponding to uint8 mask == 1 to NaN
// --------------------------------------------------------------

void mat_set_mask1_to_nan(Matrix *mat, const MatUint8 *mask)
{
    assert(mat != NULL);
    assert(mask != NULL);
    // Reset values to NaN
    size_t i = 0;
    size_t imax = mat->size1 * mat->size2;    // total size
    
    for(i = 0; i < imax; i++)
    {
        if(mask->vals[i] == 1) mat->vals[i] = MATRIX_NAN;
    }
}

// --------------------------------------------------------------
// Matrix 3d float32: Initialize, allocate, setting values to zeros
// --------------------------------------------------------------

void mat3d_float32_init(Mat3dFloat32 *mat)
{
    assert(mat != NULL);
    mat->size1 = 0;
    mat->size2 = 0;
    mat->size3 = 0;
    mat->vals =  NULL;
}

void mat3d_float32_alloc(Mat3dFloat32 *mat, size_t size1, size_t size2, size_t size3)
{
    assert(mat != NULL);
    mat->size1 = size1;    // set sizes
    mat->size2 = size2;
    mat->size3 = size3;
    mat->vals = (float32 *) calloc(size1 * size2 * size3, sizeof(float32)); // allocate elements
    assert(mat->vals !=  NULL);
}

// --------------------------------------------------------------
// Matrix 3d float32: Clear matrix, deallocating elements array
// --------------------------------------------------------------

void mat3d_float32_clear(Mat3dFloat32 *mat)
{
    assert(mat != NULL);
    if(mat->vals != NULL) 
    {
        free(mat->vals);    // deallocate array
        mat->vals = NULL;    // reset to NULL
    }
    
    mat->size1 = 0;    // reset sizes to zero 
    mat->size2 = 0;
    mat->size3 = 0;
}

float32 mat3d_float32_get(Mat3dFloat32 *mat, size_t z, size_t y, size_t x)
{
    return mat->vals[(z * mat->size2 + y) * mat->size3 + x];
}

float32 mat3d_float32_gets(Mat3dFloat32 *mat, size_t z, size_t y, size_t x)
{
    if (x >= mat->size3) return MATRIX_NAN;
    if (y >= mat->size2) return MATRIX_NAN;
    if (z >= mat->size1) return MATRIX_NAN;
    return mat->vals[(z * mat->size2 + y) * mat->size3 + x];
}

void mat3d_float32_set(Mat3dFloat32 *mat, size_t z, size_t y, size_t x, float32 val)
{
    mat->vals[(z * mat->size2 + y) * mat->size3 + x] = val;
}

void mat3d_float32_sets(Mat3dFloat32 *mat, size_t z, size_t y, size_t x, float32 val)
{
    if (x < mat->size3 && y < mat->size2 && z < mat->size1)
        mat->vals[(z * mat->size2 + y) * mat->size3 + x] = val;
}

// --------------------------------------------------------------
// Matrix 3d int16: Initialize, allocate, setting values to zeros
// --------------------------------------------------------------

void mat3d_int16_init(Mat3dInt16 *mat)
{
    assert(mat != NULL);
    mat->size1 = 0;
    mat->size2 = 0;
    mat->size3 = 0;
    mat->vals =  NULL;
}

void mat3d_int16_alloc(Mat3dInt16 *mat, size_t size1, size_t size2, size_t size3)
{
    assert(mat != NULL);
    mat3d_int16_clear(mat);
    mat->size1 = size1;    // set sizes
    mat->size2 = size2;
    mat->size3 = size3;
    mat->vals = (int16 *) calloc(size1 * size2 * size3, sizeof(int16));
    assert(mat->vals != NULL);
}

// --------------------------------------------------------------
// Matrix 3d int16: Clear matrix, deallocating elements array
// --------------------------------------------------------------

void mat3d_int16_clear(Mat3dInt16 *mat)
{
    assert(mat != NULL);
    if(mat->vals != NULL) 
    {
        free(mat->vals);    // deallocate array
        mat->vals = NULL;    // reset to NULL
    }
    mat->size1 = 0;    // reset sizes to zero 
    mat->size2 = 0;
    mat->size3 = 0;
}

int16 mat3d_int16_get(Mat3dInt16 *mat, size_t z, size_t y, size_t x)
{
    return mat->vals[(z * mat->size2 + y) * mat->size3 + x];
}

int16 mat3d_int16_gets(Mat3dInt16 *mat, size_t z, size_t y, size_t x)
{
    if (x >= mat->size3) return MATRIX_NAN;
    if (y >= mat->size2) return MATRIX_NAN;
    if (z >= mat->size1) return MATRIX_NAN;
    return mat->vals[(z * mat->size2 + y) * mat->size3 + x];
}

void mat3d_int16_set(Mat3dInt16 *mat, size_t z, size_t y, size_t x, double val)
{
    mat->vals[(z * mat->size2 + y) * mat->size3 + x] = val;
}

void mat3d_int16_sets(Mat3dInt16 *mat, size_t z, size_t y, size_t x, double val)
{
    if (x < mat->size3 && y < mat->size2 && z < mat->size1)
        mat->vals[(z * mat->size2 + y) * mat->size3 + x] = val;
}

// --------------------------------------------------------------
// Matrix 3d double: Initialize, allocate, setting values to zeros
// --------------------------------------------------------------

void mat3d_init(Mat3d *mat)
{
    assert(mat != NULL);
    mat->size1 = 0;
    mat->size2 = 0;
    mat->size3 = 0;
    mat->vals =  NULL;
}

void mat3d_alloc(Mat3d *mat, size_t size1, size_t size2, size_t size3)
{
    assert(mat != NULL);
    mat3d_clear(mat);
    mat->size1 = size1;    // set sizes
    mat->size2 = size2;
    mat->size3 = size3;
    mat->vals = (double *) calloc(size1 * size2 * size3, sizeof(double));
    assert(mat->vals != NULL);
}

// --------------------------------------------------------------
// Matrix 3d double: Clear matrix, deallocating elements array
// --------------------------------------------------------------

void mat3d_clear(Mat3d *mat)
{
    assert(mat != NULL);
    if(mat->vals != NULL) 
    {
        free(mat->vals);    // deallocate array
        mat->vals = NULL;    // reset to NULL
    }
    mat->size1 = 0;    // reset sizes to zero 
    mat->size2 = 0;
    mat->size3 = 0;
}

// Copy one major index of a 3D array to a 2D array.
void mat3d_copy_to_2d(Mat3d *mat_src, int index, Matrix *mat_dest)
{
    mat_alloc(mat_dest, mat_src->size2, mat_src->size3);
    size_t levelsize = mat_src->size2 * mat_src->size3 * sizeof(double);
    double *levelp = &mat_src->vals[index * mat_src->size2 * mat_src->size3];
    memcpy(mat_dest->vals, levelp, levelsize);
}

void mat3d_fill(Mat3d *mat, double value)
{
    int i;
    for (i = 0; i < mat->size1 * mat->size2 * mat->size3; i++)
    {
        mat->vals[i] = value;
    }
}

double mat3d_get(Mat3d *mat, size_t z, size_t y, size_t x)
{
    return mat->vals[(z * mat->size2 + y) * mat->size3 + x];
}

double mat3d_gets(Mat3d *mat, size_t z, size_t y, size_t x)
{
    if (x >= mat->size3) return MATRIX_NAN;
    if (y >= mat->size2) return MATRIX_NAN;
    if (z >= mat->size1) return MATRIX_NAN;
    return mat->vals[(z * mat->size2 + y) * mat->size3 + x];
}

void mat3d_set(Mat3d *mat, size_t z, size_t y, size_t x, double val)
{
    mat->vals[(z * mat->size2 + y) * mat->size3 + x] = val;
}

void mat3d_sets(Mat3d *mat, size_t z, size_t y, size_t x, double val)
{
    if (x < mat->size3 && y < mat->size2 && z < mat->size1)
        mat->vals[(z * mat->size2 + y) * mat->size3 + x] = val;
}

void mat3d_transpose(Mat3d *mat_prime, Mat3d *mat)
{
    mat3d_alloc(mat_prime, mat->size1, mat->size3, mat->size2);
    int i, j, k;
    for (i = 0; i < mat->size1; i++)
    {
        for (j = 0; j < mat->size2; j++)
        {
            for (k = 0; k < mat->size3; k++)
            {
                mat3d_set(mat_prime, i, k, j, mat3d_get(mat, i, j, k));
            }
        }
    }
}

// --------------------------------------------------------------
// Matrix 3d double: Copy values converting from float32 to double
// --------------------------------------------------------------

void mat3d_copy_float32(Mat3d *mat, const Mat3dFloat32 *mat_src)
{
    assert (mat != NULL);
    assert (mat_src != NULL);
    // Copy values converting element types
    size_t i = 0;
    size_t imax = mat->size1 * mat->size2 * mat->size3;    // total size
    
    for(i = 0; i < imax; i++)
    {
        mat->vals[i] = (double) mat_src->vals[i]; // float32 to double value
    }
}

// --------------------------------------------------------------
// Matrix 4d float32: Initialize, allocate, setting values to zeros
// --------------------------------------------------------------

void mat4d_float32_init(Mat4dFloat32 *mat)
{
    assert(mat != NULL);
    mat->size1 = 0;
    mat->size2 = 0;
    mat->size3 = 0;
    mat->size4 = 0;
    mat->vals =  NULL;
}

void mat4d_float32_alloc(Mat4dFloat32 *mat, size_t size1, size_t size2, size_t size3, size_t size4)
{
    assert(mat != NULL);
    mat4d_float32_clear(mat);
    mat->size1 = size1;    // set sizes
    mat->size2 = size2;
    mat->size3 = size3;
    mat->size4 = size4;
    mat->vals = (float32 *) calloc(size1 * size2 * size3 * size4, sizeof(float32)); // allocate elements
    assert(mat->vals != NULL);
}

// --------------------------------------------------------------
// Matrix 4d float32: Clear matrix, deallocating elements array
// --------------------------------------------------------------

void mat4d_float32_clear(Mat4dFloat32 *mat)
{
    assert(mat != NULL);
    if(mat->vals != NULL) 
    {
        free(mat->vals);    // deallocate array
        mat->vals = NULL;    // reset to NULL
    }
    mat->size1 = 0;    // reset sizes to zero 
    mat->size2 = 0;
    mat->size3 = 0;
    mat->size4 = 0;
}

float32 mat4D_float32_get(Mat4dFloat32 *mat, 
        size_t w, size_t z, size_t y, size_t x)
{
    return mat->vals[((w * mat->size2 + z) * mat->size3 + y) * mat->size4 + x];
}

float32 mat4D_float32_gets(Mat4dFloat32 *mat, 
        size_t w, size_t z, size_t y, size_t x)
{
    if (x >= mat->size4) return MATRIX_NAN;
    if (y >= mat->size3) return MATRIX_NAN;
    if (z >= mat->size2) return MATRIX_NAN;
    if (w >= mat->size1) return MATRIX_NAN;
    return mat->vals[((w * mat->size2 + z) * mat->size3 + y) * mat->size4 + x];
}

void mat4D_float32_set(Mat4dFloat32 *mat, 
        size_t w, size_t z, size_t y, size_t x, float32 val)
{
    mat->vals[((w * mat->size2 + z) * mat->size3 + y) * mat->size4 + x] = val;
}

void mat4D_float32_sets(Mat4dFloat32 *mat, 
        size_t w, size_t z, size_t y, size_t x, float32 val)
{
    if (x < mat->size4 && y < mat->size3 && z < mat->size2 && w < mat->size1)
        mat->vals[((w * mat->size2 + z) * mat->size3 + y) * mat->size4 + x] = val;
}

// --------------------------------------------------------------
// Matrix 4d double: Initialize, allocate, setting values to zeros
// --------------------------------------------------------------
void mat4d_init(Mat4d *mat)
{
    assert(mat != NULL);
    mat->size1 = 0;
    mat->size2 = 0;
    mat->size3 = 0;
    mat->size4 = 0;
    mat->vals =  NULL;
}

void mat4d_alloc(Mat4d *mat, size_t size1, size_t size2, size_t size3, size_t size4)
{
    assert(mat != NULL);
    mat4d_clear(mat);
    mat->size1 = size1;    // set sizes
    mat->size2 = size2;
    mat->size3 = size3;
    mat->size4 = size4;
    mat->vals = (double *) calloc(size1 * size2 * size3 * size4, sizeof(double));
    assert(mat->vals != NULL);
}

// --------------------------------------------------------------
// Matrix 4d double: Clear matrix, deallocating elements array
// --------------------------------------------------------------

void mat4d_clear(Mat4d *mat)
{
    assert(mat != NULL);
    if(mat->vals != NULL) 
    {
        free(mat->vals);    // deallocate array
        mat->vals = NULL;    // reset to NULL
    }
    mat->size1 = 0;    // reset sizes to zero 
    mat->size2 = 0;
    mat->size3 = 0;
    mat->size4 = 0;
}

double mat4d_get(Mat4d *mat, size_t w, size_t z, size_t y, size_t x)
{
    return mat->vals[((w * mat->size2 + z) * mat->size3 + y) * mat->size4 + x];
}

double mat4d_gets(Mat4d *mat, size_t w, size_t z, size_t y, size_t x)
{
    if (x >= mat->size4) return MATRIX_NAN;
    if (y >= mat->size3) return MATRIX_NAN;
    if (z >= mat->size2) return MATRIX_NAN;
    if (w >= mat->size1) return MATRIX_NAN;
    return mat->vals[((w * mat->size2 + z) * mat->size3 + y) * mat->size4 + x];
}

void mat4d_set(Mat4d *mat, size_t w, size_t z, size_t y, size_t x, double val)
{
    mat->vals[((w * mat->size2 + z) * mat->size3 + y) * mat->size4 + x] = val;
}

void mat4d_sets(Mat4d *mat, size_t w, size_t z, size_t y, size_t x, double val)
{
    if (x < mat->size4 && y < mat->size3 && z < mat->size2 && w < mat->size1)
        mat->vals[((w * mat->size2 + z) * mat->size3 + y) * mat->size4 + x] = val;
}

// --------------------------------------------------------------
// Matrix 4d double: Copy values converting from float32 to double
// --------------------------------------------------------------

void mat4d_copy_float32(Mat4d *mat, const Mat4dFloat32 *mat_src)
{
    assert(mat != NULL);
    assert(mat_src != NULL);
    // Copy values converting element types
    size_t i = 0;
    size_t imax = mat->size1 * mat->size2 * mat->size3 * mat->size4;    // total size
    
    for(i = 0; i < imax; i++)
    {
        mat->vals[i] = (double) mat_src->vals[i]; // float32 to double value
    }
}

void clamp1d(Vector *v, double minv, double maxv)
{
    int i;
    for (i = 0; i < v->size; i++)
    {
        if (v->vals[i] < minv)
            v->vals[i] = minv;
        else if (v->vals[i] > maxv)
            v->vals[i] = maxv;
    }
}

void clamp2d(Matrix *v, double minv, double maxv)
{
    int i;
    for (i = 0; i < v->size1 * v->size2; i++)
    {
        if (v->vals[i] < minv)
            v->vals[i] = minv;
        else if (v->vals[i] > maxv)
            v->vals[i] = maxv;
    }
}

void clamp3d(Mat3d *v, double minv, double maxv)
{
    int i;
    for (i = 0; i < v->size1 * v->size2 * v->size3; i++)
    {
        if (v->vals[i] < minv)
            v->vals[i] = minv;
        else if (v->vals[i] > maxv)
            v->vals[i] = maxv;
    }
}


void clamp4d(Mat4d *v, double minv, double maxv)
{
    int i;
    for (i = 0; i < v->size1 * v->size2 * v->size3 * v->size4; i++)
    {
        if (v->vals[i] < minv)
            v->vals[i] = minv;
        else if (v->vals[i] > maxv)
            v->vals[i] = maxv;
    }
}

void mat_flipdim(Matrix *mat)
{
    size_t rowsize = mat->size2 * sizeof(double);
    double *temp = (double*)malloc(rowsize);
    int i, j;
    for (i = 0; i < mat->size1 / 2; i++)
    {
        j = mat->size1 - 1 - i;
        memcpy(temp, &mat->vals[i * mat->size2], rowsize);
        memcpy(&mat->vals[i * mat->size2], &mat->vals[j  * mat->size2], rowsize);
        memcpy(&mat->vals[j * mat->size2], temp, rowsize);
    }
    free(temp);
}

void mat_uint8_flipdim(MatUint8 *mat)
{
    size_t rowsize = mat->size2 * sizeof(uint8);
    uint8 *temp = (uint8*)malloc(rowsize);
    int i, j;
    for (i = 0; i < mat->size1 / 2; i++)
    {
        j = mat->size1 - 1 - i;
        memcpy(temp, &mat->vals[i * mat->size2], rowsize);
        memcpy(&mat->vals[i * mat->size2], &mat->vals[j  * mat->size2], rowsize);
        memcpy(&mat->vals[j * mat->size2], temp, rowsize);
    }
    free(temp);
}

void mat3d_flipdim(Mat3d *mat)
{
    size_t levelix = mat->size2 * mat->size3;
    size_t levelsize = levelix * sizeof(double);
    double *temp = (double*)malloc(levelsize);
    int i, j;
    for (i = 0; i < mat->size1 / 2; i++)
    {
        j = mat->size1 - 1 - i;
        memmove(temp, &mat->vals[i * levelix], levelsize);
        memmove(&mat->vals[i * levelix], &mat->vals[j * levelix], levelsize);
        memmove(&mat->vals[j * levelix], temp, levelsize);
    }
    free(temp);
}

double mat_min(Matrix *mat)
{
    int i;
    double minval = mat->vals[0];
    for (i = 1; i < mat->size1*mat->size2; i++)
    {
        if (mat->vals[i] < minval || is_nan(minval))
            minval = mat->vals[i];
    }
    return minval;
}

double mat_max(Matrix *mat)
{
    int i;
    double maxval = mat->vals[0];
    for (i = 1; i < mat->size1*mat->size2; i++)
    {
        if (mat->vals[i] > maxval || is_nan(maxval))
            maxval = mat->vals[i];
    }
    return maxval;
}

void mat_set_fillvalues_to_nans(Matrix *mat, double fillvalue)
{
    int i;
    for (i = 0; i < mat->size1 * mat->size2; i++)
    {
        if (dequal(mat->vals[i], fillvalue))
            mat->vals[i] = REAL_NAN;
    }
}
