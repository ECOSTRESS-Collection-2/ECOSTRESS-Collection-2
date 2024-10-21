#pragma once
/** 
 * @file Matrix and Vector types and functions used for Modis, VIIRs, 
 * and other general purpose HDF and NetCDF data handling.
 * @author Robert Freepartner, JPL, January 2016
 * @author based on TES matrix by Simon Latyshev, Raytheon
 *
 * Copyright (c) JPL, All rights reserved
 *
 * List of types:
 *  Vector          One dimension array of double values
 *  VecFloat32      One dimension array of float32 values
 *  Matrix          2D matrix of double
 *  MatUint8        2D matrix of 8-bit unsigned
 *  MatInt16        2D matrix of 16-bit signed
 *  MatUint16       2D matrix of 16-bit unsigned
 *  MatFloat32      2D matrix of float32
 *  Mat3d           3D matrix of double
 *  Mat3dInt16      3D matrix of int16
 *  Mat3dFloat32    3D matrix of float32
 *  Mat4d           4D matrix of double
 *  Mat4dFloat32    4D matrix of float32
 */

#include <hdfi.h>
#include <unistd.h>

#define MATRIX_NAN -32768
#define MATRIX_UINT8_NAN ((uint8)0x80)

// --------------------------------------------------------------
// Vector 1d, type float32 (added for HDF arrays)
// --------------------------------------------------------------

typedef struct 
{
    size_t size;    // vector size
    float32 *vals;    // values array
} VecFloat32;

// Initialize, allocate, setting values to zeros
void vec_float32_init(VecFloat32 *vec);
void vec_float32_alloc(VecFloat32 *vec, size_t size);

// Clear, deallocating values array
void vec_float32_clear(VecFloat32 *vec);

// Access
float32 vec_float32_get(VecFloat32 *vec, size_t x); 
void vec_float32_set(VecFloat32 *vec, size_t x, float32 val); 

// --------------------------------------------------------------
// Vector 1d, type double
// --------------------------------------------------------------

typedef struct 
{
    size_t size;    // vector size
    double *vals;    // values array
} Vector;

// Initialize, allocate, setting values to zeros
void vec_init(Vector *vec);
void vec_alloc(Vector *vec, size_t size);

// Clear, deallocating values array
void vec_clear(Vector *vec);

// Copy values
void vec_copy(Vector *dest, Vector *source);

// Copy values converting from float32 to double
void vec_copy_float32(Vector *vec, const VecFloat32 *vec_src);

// Access
double vec_get(Vector *vec, size_t x);
void vec_set(Vector *vec, size_t x, double val);

// --------------------------------------------------------------
// Matrix 2d, type uint8
// --------------------------------------------------------------

typedef struct 
{
    size_t size1;    // size in the 1st dimension
    size_t size2;    // size in the 2nd dimension
    uint8 *vals;    // values array
} MatUint8;

// Initialize, allocate values array
void mat_uint8_init(MatUint8 *mat);
void mat_uint8_alloc(MatUint8 *mat, size_t size1, size_t size2);

// Clear, deallocate values array
void mat_uint8_clear(MatUint8 *mat);

// Fill
void mat_uint8_fill(MatUint8 *mat, uint8 value);

// Access 
uint8 mat_uint8_get(MatUint8 *mat, size_t y, size_t x);
void mat_uint8_set(MatUint8 *mat, size_t y, size_t x, uint8 val);

// Invert the rows of the matrix
void mat_uint8_flipdim(MatUint8 *mat);

// --------------------------------------------------------------
// Matrix 2d, type int16 (added for HDF arrays)
// --------------------------------------------------------------

typedef struct 
{
    size_t size1;    // size in the 1st dimension
    size_t size2;    // size in the 2nd dimension
    int16  *vals;    // values array
} MatInt16;

// Initialize, allocate values array
void mat_int16_init(MatInt16 *mat);
void mat_int16_alloc(MatInt16 *mat, size_t size1, size_t size2);

// Clear, deallocate values array
void mat_int16_clear(MatInt16 *mat);

// Access 
int16 mat_int16_get(MatInt16 *mat, size_t y, size_t x);
void mat_int16_set(MatInt16 *mat, size_t y, size_t x, int16 val);

// Int32
typedef struct 
{
    size_t size1;    // size in the 1st dimension
    size_t size2;    // size in the 2nd dimension
    int32  *vals;    // values array
} MatInt32;

// Initialize, allocate values array
void mat_int32_init(MatInt32 *mat);
void mat_int32_alloc(MatInt32 *mat, size_t size1, size_t size2);

// Clear, deallocate values array
void mat_int32_clear(MatInt32 *mat);

// Access 
int32 mat_int32_get(MatInt32 *mat, size_t y, size_t x);
void mat_int32_set(MatInt32 *mat, size_t y, size_t x, int32 val);

// --------------------------------------------------------------
// Matrix 2d, type uint16
// --------------------------------------------------------------

typedef struct 
{
    size_t size1;    // size in the 1st dimension
    size_t size2;    // size in the 2nd dimension
    uint16  *vals;   // values array
} MatUint16;

// Initialize, allocate values array
void mat_uint16_init(MatUint16 *mat);
void mat_uint16_alloc(MatUint16 *mat, size_t size1, size_t size2);

// Clear, deallocate values array
void mat_uint16_clear(MatUint16 *mat);

// Access 
uint16 mat_uint16_get(MatUint16 *mat, size_t y, size_t x);
void mat_uint16_set(MatUint16 *mat, size_t y, size_t x, uint16 val);

// --------------------------------------------------------------
// Matrix 2d, type float32 (added for HDF arrays)
// --------------------------------------------------------------

typedef struct 
{
    size_t size1;    // size in the 1st dimension
    size_t size2;    // size in the 2nd dimension
    float32 *vals;    // values arra
} MatFloat32;

// Initialize, allocate values array
void mat_float32_init(MatFloat32 *mat);
void mat_float32_alloc(MatFloat32 *mat, size_t size1, size_t size2);

// Clear, deallocate values array
void mat_float32_clear(MatFloat32 *mat);

// Access 
float32 mat_float32_get(MatFloat32 *mat, size_t y, size_t x);
void mat_float32_set(MatFloat32 *mat, size_t y, size_t x, float32 val);

// --------------------------------------------------------------
// Matrix 2d, type double
// --------------------------------------------------------------

typedef struct 
{
    size_t size1;    // size in the 1st dimension
    size_t size2;    // size in the 2nd dimension
    double *vals;    // values, elements array
} Matrix;

// Initialize, allocate, setting values to zeros
void mat_init(Matrix *mat);
void mat_alloc(Matrix *mat, size_t size1, size_t size2);

// Clear matrix, deallocating elements array
void mat_clear(Matrix *mat);

// Copy values
void mat_copy(Matrix *dest, Matrix *source);

// Copy values converting from int16 to double
void mat_copy_int16(Matrix *mat, const int16 *mat_src);

// Copy values converting from uint16 to double
void mat_copy_uint16(Matrix *mat, const uint16 *mat_src);

// Copy values converting from uint8 to uint16
void mat_copy_uint_8_to_16(MatUint16 *mat, const MatUint8 *mat_src);

// Copy values converting from float32 to double
void mat_copy_float32(Matrix *mat, const MatFloat32 *mat_src);

// Fill a matrix with a constant value.
void mat_fill(Matrix *mat, double value);

// Set values corresponding to uint8 mask == 1 to NaN
void mat_set_mask1_to_nan(Matrix *mat, const MatUint8 *mask);

// Access 
double mat_get(Matrix *mat, size_t y, size_t x);
void mat_set(Matrix *mat, size_t y, size_t x, double val);

// Invert the rows of the matrix
void mat_flipdim(Matrix *mat);

// Transpose = create a N x M array from an M x N array and flip the indexes.
void mat_transpose(Matrix *mat_prime, Matrix *mat);

// --------------------------------------------------------------
// Matrix 3d, type float32 (added for HDF arrays)
// --------------------------------------------------------------

typedef struct 
{
    size_t size1;    // size in the 1st dimension
    size_t size2;    // size in the 2nd dimension
    size_t size3;    // size in the 3rd dimension
    float32 *vals;    // values array
} Mat3dFloat32;

// Initialize, allocate, setting values to zeros
void mat3d_float32_init(Mat3dFloat32 *mat);
void mat3d_float32_alloc(Mat3dFloat32 *mat, size_t size1, size_t size2, size_t size3);

// Clear matrix, deallocating elements array
void mat3d_float32_clear(Mat3dFloat32 *mat);

// Access 
float32 mat3d_float32_get(Mat3dFloat32 *mat, size_t z, size_t y, size_t x);
void mat3d_float32_set(Mat3dFloat32 *mat, size_t z, size_t y, size_t x, float32 val);

// --------------------------------------------------------------
// Matrix 3d, type int16
// --------------------------------------------------------------

typedef struct 
{
    size_t size1;    // size in the 1st dimension
    size_t size2;    // size in the 2nd dimension
    size_t size3;    // size in the 3rd dimension
    int16  *vals;    // values array
} Mat3dInt16;

// Initialize, allocate, setting values to zeros
void mat3d_int16_init(Mat3dInt16 *mat);
void mat3d_int16_alloc(Mat3dInt16 *mat, size_t size1, size_t size2, size_t size3);

// Clear matrix, deallocating elements array
void mat3d_int16_clear(Mat3dInt16 *mat);

// Access 
int16 mat3d_int16_get(Mat3dInt16 *mat, size_t z, size_t y, size_t x);
void mat3d_int16_set(Mat3dInt16 *mat, size_t z, size_t y, size_t x, double val);

// --------------------------------------------------------------
// Matrix 3d, type double
// --------------------------------------------------------------

typedef struct 
{
    size_t size1;    // size in the 1st dimension
    size_t size2;    // size in the 2nd dimension
    size_t size3;    // size in the 3rd dimension
    double *vals;    // values array
} Mat3d;

// Initialize, allocate, setting values to zeros
void mat3d_init(Mat3d *mat);
void mat3d_alloc(Mat3d *mat, size_t size1, size_t size2, size_t size3);

// Clear matrix, deallocating elements array
void mat3d_clear(Mat3d *mat);

// Copy values converting from float32 to double
void mat3d_copy_float32(Mat3d *mat, const Mat3dFloat32 *mat_src);

// Copy one major index of a 3D array to a 2D array.
void mat3d_copy_to_2d(Mat3d *mat_src, int index, Matrix *mat_dest);

// Fill the 3D array with a value
void mat3d_fill(Mat3d *mat, double value);

// Access 
double mat3d_get(Mat3d *mat, size_t z, size_t y, size_t x);
void mat3d_set(Mat3d *mat, size_t z, size_t y, size_t x, double val);

// Flip Dim - flips the levels of a 3D array.
void mat3d_flipdim(Mat3d *mat);

// Transpose = create a K x N x M array from an K x M x N array and flip the indexes.
void mat3d_transpose(Mat3d *mat_prime, Mat3d *mat);

// --------------------------------------------------------------
// Matrix 4d, type float32 (added for HDF arrays)
// --------------------------------------------------------------

typedef struct 
{
    size_t size1;    // size in the 1st dimension
    size_t size2;    // size in the 2nd dimension
    size_t size3;    // size in the 3rd dimension
    size_t size4;    // size in the 4th dimension
    float32 *vals;    // values array
} Mat4dFloat32;

// Initialize, allocate, setting values to zeros
void mat4d_float32_init(Mat4dFloat32 *mat);
void mat4d_float32_alloc(Mat4dFloat32 *mat, size_t size1, size_t size2, size_t size3, size_t size4);

// Clear matrix, deallocating elements array
void mat4d_float32_clear(Mat4dFloat32 *mat);

// Access 
float32 mat4D_float32_get(Mat4dFloat32 *mat, size_t w, size_t z, size_t y, size_t x);
void mat4D_float32_set(Mat4dFloat32 *mat, size_t w, size_t z, size_t y, size_t x, float32 val);

// --------------------------------------------------------------
// Matrix 4d, type double
// --------------------------------------------------------------

typedef struct 
{
    size_t size1;    // size in the 1st dimension
    size_t size2;    // size in the 2nd dimension
    size_t size3;    // size in the 3rd dimension
    size_t size4;    // size in the 4th dimension
    double *vals;    // values array
} Mat4d;

// Initialize, allocate, setting values to zeros
void mat4d_init(Mat4d *mat);
void mat4d_alloc(Mat4d *mat, size_t size1, size_t size2, size_t size3, size_t size4);

// Clear matrix, deallocating elements array
void mat4d_clear(Mat4d *mat);

// Copy values converting from float32 to double
void mat4d_copy_float32(Mat4d *mat, const Mat4dFloat32 *mat_src);

// Access
double mat4d_get(Mat4d *mat, size_t w, size_t z, size_t y, size_t x);
void mat4d_set(Mat4d *mat, size_t w, size_t z, size_t y, size_t x, double val);

// --------------------------------------------------------------
// Clamping functions
// --------------------------------------------------------------

void clamp1d(Vector *v, double minv, double maxv);
void clamp2d(Matrix *v, double minv, double maxv);
void clamp3d(Mat3d *v, double minv, double maxv);
void clamp4d(Mat4d *v, double minv, double maxv);

// --------------------------------------------------------------
// Utility functions
// --------------------------------------------------------------

double mat_min(Matrix *mat);
double mat_max(Matrix *mat);
void mat_set_fillvalues_to_nans(Matrix *mat, double fillvalue);
