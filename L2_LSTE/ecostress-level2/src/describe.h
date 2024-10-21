#pragma once
/**
 * @file This header file contains definitions and prototypes for the functions 
 * that describe data.
 *
 * Each version of a describe function outputs the name and dimensions,
 * then calls describe_data to output the number of nans, zeroes, +inf, -inf, 
 * posirtive, and negative values. Positive and negative values include 
 * the mean values for each and the combined mean, excluding zeroes.
 *
 * @author Robert Freepartner, JaDa Systems/Raytheon/JPL
 *
 * @copyright (c) Copyright 2016, Jet Propulsion Laboratories, Pasadena, CA
 */
#include "matrix.h"

/**
 * @brief output a description of the data to stdout
 * @param   pd  Pointer to a vector of double values
 * @param   nd  Number of values in the data to be described
 */
void describe_data(double *pd, size_t nd);

/**
 * @brief Describes data that is arranged in nrows of consecutive
 * column values.
 * @param   name    A label to output for the data
 * @param   pd      Pointer to the data
 * @param   nrows   Number of rows in the data
 * @param   ncols   Nummber of columns in each row
 */
void describe_data_2d(const char *name, double *pd, size_t nrows, size_t ncols);
void describe_matrix(const char *name, Matrix *mat);
void describe_mat_u8(const char *name, MatUint8 *mat);
void describe_mat_u16_bitset(const char *name, MatUint16 *mat);
void describe_mat_u16(const char *name, MatUint16 *mat);
void describe_mat3d(const char *name, Mat3d *mat);
void describe_data_3d(const char *name, double *pd[], 
    size_t nlevels, size_t nrows, size_t ncols);
void describe_vector(const char *name, Vector *vec);

