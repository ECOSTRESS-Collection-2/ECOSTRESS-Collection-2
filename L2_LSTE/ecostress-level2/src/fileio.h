/// @file File I/O interface for HDF4, HDF5, and NetCDF files.
//
//  References: 
/// @author Robert Freepartner, JPL, [BF]
//  Based on tes_fileIO C programming by:
/// @author Simon Latyshev, Raytheon, [SL]
//
// Copyright (c) JPL, All rights reserved
 
#pragma once

#include <hdf.h>
#include <hdf5.h>
#include <mfhdf.h>
#include "matrix.h"

#define FILEIO_MAX_FILENAME 256
#define FILEIO_MAX_ATTRNAME 256
#define FILEIO_STR_SIZE 1024
#define FILEIO_NAN MATRIX_NAN

#ifndef NC_NOWRITE
#define NC_NOWRITE H5F_ACC_RDONLY
#endif

#ifndef NC_WRITE
#define NC_WRITE H5F_ACC_RDWR
#endif

typedef enum
{
    unknown,
    hdf4,
    hdf5,
    netcdf
} FileType;

/// @brief Returns type of file based on filename.
///
/// @param fname    file name string
/// @return         file type
FileType get_file_type(const char *fname);

// --------------------------------------------------------------
// HDF File I/O 
// --------------------------------------------------------------
typedef int32 hdf4id_t;

/// @brief open an hdf4 file
/// 
/// @param fname    File name string
/// @param  access_mode: DFACC_READ, DFACC_WRITE, or DFACC_CREATE
/// @return Failure: FAIL, else hdf file id
hdf4id_t open_hdf(const char* fname, int32 access_mode);

/// @brief creatge an hdf4 file
/// 
/// @param fname    File name string
/// @return Failure: FAIL, else hdf file id
hdf4id_t create_hdf(const char* fname);

/// @brief close an hdf4 file
//
/// @param hdfid     hdf file id
void close_hdf(hdf4id_t hdfid);

/// @brief Read HDF4 1D array into a vector of doubles
/// @param sd_fid   HDF open file handle from SDstart
/// @param dataset  name of the dataset to read
/// @param vec      Initialized vector
/// @return 0=success -1=error
int hdfread_vec(hdf4id_t sd_fid, const char *dataset, Vector *vec);

/// @brief Read HDF4 2D dataset into matrix of doubles
/// @param sd_fid   HDF open file handle from SDstart
/// @param dataset  name of the dataset to read
/// @param mat      Initialized matrix
/// @return 0=success -1=error
int hdfread_mat(hdf4id_t sd_fid, const char *dataset, Matrix *mat);

/// @brief Read HDF4 2D dataset into matrix of uint8
/// @param sd_fid   HDF open file handle from SDstart
/// @param dataset  name of the dataset to read
/// @param mat      Initialized matrix
/// @return 0=success -1=error
int hdfread_mat_uint8(hdf4id_t sd_fid, const char *dataset, MatUint8 *mat);

/// @brief Read HDF4 3D Float32 dataset into a 3D matrix of doubles
/// @param sd_fid   HDF open file handle from SDstart
/// @param dataset  name of the dataset to read
/// @param mat      Initialized 3D matrix
/// @return 0=success -1=error
int hdfread_mat3d(hdf4id_t sd_fid, const char *dataset, Mat3d *mat);

/// @brief Read HDF4 4D Float32 dataset into a 3D matrix of doubles
/// @param sd_fid   HDF open file handle from SDstart
/// @param dataset  name of the dataset to read
/// @param mat      Initialized 4D matrix
/// @return 0=success -1=error
int hdfread_mat4d(hdf4id_t sd_fid, const char *dataset, Mat4d *mat);

/// @brief Read scale and offset values from a 2D dataset and apply them.
/// @param sd_fid   HDF open file handle from SDstart
/// @param dataset  name of the dataset to be scaled
/// @param mat[in,out]    Initialized matrix
/// @return 0=success -1=error
int scale_and_offset_mat(hdf4id_t sd_fid, const char *dataset, Matrix *mat);

/// @brief Write HDF4 1D array from a vector of doubles
/// @param sd_fid   HDF open file handle from SDstart
/// @param dataset  name of the dataset to write
/// @param vec      Initialized vector
/// @return 0=success -1=error
int hdfwrite_vec(hdf4id_t sd_fid, const char *dataset, Vector *vec);

/// @brief Write HDF4 2D dataset from a matrix of doubles
/// @param sd_fid   HDF open file handle from SDstart
/// @param dataset  name of the dataset to write
/// @param mat      Initialized matrix
/// @return 0=success -1=error
int hdfwrite_mat(hdf4id_t sd_fid, const char *dataset, Matrix *mat);

/// @brief Write HDF4 2D dataset from a matrix of uint8
/// @param sd_fid   HDF open file handle from SDstart
/// @param dataset  name of the dataset to write
/// @param mat      Initialized matrix
/// @return 0=success -1=error
int hdfwrite_mat_uint8(hdf4id_t sd_fid, const char *dataset, MatUint8 *mat);

/// @brief Write HDF4 2D dataset from a matrix of uint16
/// @param sd_fid   HDF open file handle from SDstart
/// @param dataset  name of the dataset to write
/// @param mat      Initialized matrix
/// @return 0=success -1=error
int hdfwrite_mat_uint16(hdf4id_t sd_fid, const char *dataset, MatUint16 *mat);

/// @brief Write HDF4 3D dataset from a 3D matrix of doubles
/// @param sd_fid   HDF open file handle from SDstart
/// @param dataset  name of the dataset to write
/// @param mat      Initialized 3D matrix
/// @return 0=success -1=error
int hdfwrite_mat3d(hdf4id_t sd_fid, const char *dataset, Mat3d *mat);

/// @brief Write HDF4 4D dataset from a 3D matrix of doubles
/// @param sd_fid   HDF open file handle from SDstart
/// @param dataset  name of the dataset to write
/// @param mat      Initialized 4D matrix
/// @return 0=success -1=error
int hdfwrite_mat4d(hdf4id_t sd_fid, const char *dataset, Mat4d *mat);

// --------------------------------------------------------------
// HDF5 File I/O 
// --------------------------------------------------------------

/// @brief Specifies chunking/compression to be used for 2D dataset creation.
/// @param dim1         chunk dimension 1
/// @param dim2         chunk dimension 2
/// @param deflate      0 (no compression) to 9 (max compression)
void set_hdf5_2d_compression_on(hsize_t dim1, hsize_t dim2, uint deflate);

/// @brief Remove chunking/compression for 2D output.
void set_hdf5_2d_compression_off();

/// @brief Open an HDF5 file
/// @param fpname       file pathname
/// @param nc_attrs     HDF5 attributes for file:
///                     H5F_ACC_RDONLY, or H5F_ACC_RDWR
/// @return h5fid file handle, -1 = failure
hid_t open_hdf5(const char *fpname, unsigned nc_attrs);

/// @brief Create an HDF5 file
/// @param fpname       file pathname
/// @return h5fid file handle, -1 = failure
hid_t create_hdf5(const char *fpname);

// @brief Close an HDF5 file
// @param h5fid  File ID from open_hdf5
void close_hdf5(hid_t h5fid);

/// @brief Read HDF5 1D array into a vector of doubles
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to read
/// @param vec      Initialized vector 
int hdf5read_vec(hid_t h5fid, const char *dataset, Vector *vec);

/// @brief Read HDF5 2D array into a matrix of doubles
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to read
/// @param mat      Initialized matrix 
int hdf5read_mat(hid_t h5fid, const char *dataset, Matrix *mat);

/// @brief Read HDF5 2D Uint8 array into a 2D matrix of int16
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to read
/// @param mat      Initialized 2D matrix 
int hdf5read_mat_uint8(int hid_t, const char *dataset, MatUint8 *mat);

/// @brief Read HDF5 2D Int16 array into a 2D matrix of int16
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to read
/// @param mat      Initialized 2D matrix 
int hdf5read_mat_int16(int hid_t, const char *dataset, MatInt16 *mat);

/// @brief Read HDF5 2D Int32 array into a 2D matrix of int32
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to read
/// @param mat      Initialized 2D matrix 
int hdf5read_mat_int32(int hid_t, const char *dataset, MatInt32 *mat);

/// @brief Read HDF5 2D Uint16 array into a 2D matrix of uint16
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to read
/// @param mat      Initialized 2D matrix 
int hdf5read_mat_uint16(int hid_t, const char *dataset, MatUint16 *mat);

/// @brief Read HDF5 3D Float32 array into a 3D matrix of doubles
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to read
/// @param mat      Initialized 3D matrix 
int hdf5read_mat3d(hid_t h5fid, const char *dataset, Mat3d *mat);

/// @brief Read HDF5 3D Uint16 array into a 3D matrix of int16
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to read
/// @param mat      Initialized 3D matrix 
int hdf5read_mat3d_int16(int hid_t, const char *dataset, Mat3dInt16 *mat);

/// @brief Read HDF5 4D Float32 array into a 4D matrix of doubles
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to read
/// @param mat      Initialized 4D matrix 
int hdf5read_mat4d(hid_t h5fid, const char *dataset, Mat4d *mat);

/// @brief write HDF5 1D from a vector of doubles
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to write
/// @param vec      Dataset vector 
int hdf5write_vec(hid_t h5fid, const char *dataset, Vector *vec);

/// @brief write HDF5 2D from a matrix of doubles
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to write
/// @param mat      Dataset matrix 
int hdf5write_mat(hid_t h5fid, const char *dataset, Matrix *mat);

/// @brief write HDF5 2D from a matrix of unsigned int8
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to write
/// @param mat      Dataset matrix 
int hdf5write_mat_uint8(hid_t h5fid, const char *dataset, MatUint8 *mat);

/// @brief write HDF5 2D from a matrix of int16
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to write
/// @param mat      Dataset matrix 
int hdf5write_mat_int16(hid_t h5fid, const char *dataset, MatInt16 *mat);

/// @brief write HDF5 2D from a matrix of unsigned int16
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to write
/// @param mat      Dataset matrix 
int hdf5write_mat_uint16(hid_t h5fid, const char *dataset, MatUint16 *mat);

/// @brief write HDF5 2D as float32 from a matrix of double
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to write
/// @param mat      Dataset matrix 
int hdf5write_mat_asfloat32(hid_t h5fid, const char *dataset, Matrix *mat);

/// @brief write HDF5 3D from a 3D matrix of doubles
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to write
/// @param mat      Dataset 3D matrix 
int hdf5write_mat3d(hid_t h5fid, const char *dataset, Mat3d *mat);

/// @brief write HDF5 4D from a 4D matrix of doubles
/// @param h5fid    H5F open file handle from open_hdf5
/// @param dataset  name of the dataset to write
/// @param mat      Dataset 4D matrix 
int hdf5write_mat4d(hid_t h5fid, const char *dataset, Mat4d *mat);

/// @brief Read scale and offset values from a 2D dataset using hdf5.
/// @param h5fid        HDF5 open file handle
/// @param dataset      name of the dataset to be scaled
/// @param mat[in,out]  Initialized matrix
/// @return 0=success -1=error
int hdf5_scale_and_offset_mat(hid_t h5fid, const char *dataset, Matrix *mat);

// --------------------------------------------------------------
// NetCDF File I/O 
// --------------------------------------------------------------

/// @brief Open a NetCDF file using HDF5
/// @param fpname       file pathname
/// @param nc_attrs     HDF5 attributes for file:
///                     H5F_ACC_RDONLY, or H5F_ACC_RDWR
/// @return ncid file handle, -1 = failure
hid_t open_netcdf(const char *fpname, unsigned nc_attrs);

// @brief Close a NetCDF file
// @param ncid  File ID from nc_open or open_nettcdf
void close_netcdf(hid_t ncid);

/// @brief Read NetCDF 1D array into a vector of doubles
/// @param ncid     NetCDF open file handle from nc_open or open_netcdf
/// @param dataset  name of the dataset to read
/// @param vec      Initialized vector 
int netcdfread_vec(hid_t ncid, const char *dataset, Vector *vec);

/// @brief Read NetCDF 2D array into a matrix of doubles
/// @param ncid     NetCDF open file handle from nc_open or open_netcdf
/// @param dataset  name of the dataset to read
/// @param mat      Initialized matrix 
int netcdfread_mat(hid_t ncid, const char *dataset, Matrix *mat);

/// @brief Read NetCDF 3D Float32 array into a 3D matrix of doubles
/// @param ncid     NetCDF open file handle from nc_open or open_netcdf
/// @param dataset  name of the dataset to read
/// @param mat      Initialized 3D matrix 
int netcdfread_mat3d(hid_t ncid, const char *dataset, Mat3d *mat);

/// @brief Read NetCDF 4D Float32 array into a 4D matrix of doubles
/// @param ncid     NetCDF open file handle from nc_open or open_netcdf
/// @param dataset  name of the dataset to read
/// @param mat      Initialized 4D matrix 
int netcdfread_mat4d(hid_t ncid, const char *dataset, Mat4d *mat);

/// @brief write NetCDF 1D array into a vector of doubles
/// @param ncid     NetCDF open file handle from nc_open or open_netcdf
/// @param dataset  name of the dataset to write
/// @param vec      Initialized vector 
int netcdfwrite_vec(hid_t ncid, const char *dataset, Vector *vec);

/// @brief write NetCDF 2D array into a matrix of doubles
/// @param ncid     NetCDF open file handle from nc_open or open_netcdf
/// @param dataset  name of the dataset to write
/// @param mat      Initialized matrix 
int netcdfwrite_mat(hid_t ncid, const char *dataset, Matrix *mat);

/// @brief write NetCDF 3D Float32 array into a 3D matrix of doubles
/// @param ncid     NetCDF open file handle from nc_open or open_netcdf
/// @param dataset  name of the dataset to write
/// @param mat      Initialized 3D matrix 
int netcdfwrite_mat3d(hid_t ncid, const char *dataset, Mat3d *mat);

/// @brief write NetCDF 4D Float32 array into a 4D matrix of doubles
/// @param ncid     NetCDF open file handle from nc_open or open_netcdf
/// @param dataset  name of the dataset to write
/// @param mat      Initialized 4D matrix 
int netcdfwrite_mat4d(hid_t ncid, const char *dataset, Mat4d *mat);
