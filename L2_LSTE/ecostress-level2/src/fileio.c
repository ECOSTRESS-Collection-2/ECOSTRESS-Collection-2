// File I/O interface for HDF4, HDF5, and NetCDF files.
//
// See fileio.h for credits and documentation
//
// Copyright (c) JPL, All rights reserved

//#define DEBUG
//#define COMMENT

// Depending on which version of HDF5 is installed, the open call
// may need the second parameter or not. If there is a compiler error, 
// such as error: too many arguments to function 'H5Dopen1',
// then change the define for OLD_HDF5:
//#define OLD_HD5 
#ifdef OLD_HD5
#define H5D_OPEN
#else
#define H5D_OPEN ,H5P_DEFAULT
#endif

#include <assert.h>
#include "fileio.h"
#include "lste_lib.h"

#define DEFAULT_FILL_VALUE 9.9999999E14

FileType get_file_type(const char *fname)
{
    const char *find_dot = strrchr(fname, '.');
    if (find_dot == NULL) return unknown;
    if (strcasecmp(find_dot, ".hdf") == 0) return hdf4;
    if (strcasecmp(find_dot, ".h5") == 0) return hdf5;
    if (strcasecmp(find_dot, ".nc") == 0) return netcdf;
    // Less common file extensions
    if (strcasecmp(find_dot, ".nc4") == 0) return netcdf;
    if (strcasecmp(find_dot, ".cdf") == 0) return netcdf;
    if (strcasecmp(find_dot, ".netcdf") == 0) return netcdf;
    if (strcasecmp(find_dot, ".hdf5") == 0) return hdf5;
    if (strcasecmp(find_dot, ".he5") == 0) return hdf5;
    if (strcasecmp(find_dot, ".h4") == 0) return hdf4;
    if (strcasecmp(find_dot, ".hdf4") == 0) return hdf4;
    if (strcasecmp(find_dot, ".he2") == 0) return hdf4;
    return unknown;
}

hdf4id_t open_hdf(const char* fname, int32 access_mode)
{
   return SDstart(fname, access_mode);
}

hdf4id_t create_hdf(const char* fname)
{
    return open_hdf(fname, DFACC_CREATE);
}

void close_hdf(hdf4id_t hdfid)
{
    SDend(hdfid);
}

// --------------------------------------------------------------
// HDF File I/O: Read HDF dataset containing 1d array
// --------------------------------------------------------------

// Read HDF 1d array to double vector
// 1. Assuming that file is open already
// 2. Assuming that the vector is initialized by an _init call
// 3. Allocating and initializing double vector

int hdfread_vec(hdf4id_t sd_fid, const char *dataset, Vector *vec)
{
    // Open dataset
    hdf4id_t sds_index = SDnametoindex(sd_fid, dataset);
    if (sds_index == FAIL)
    {
        fprintf(stderr, 
                "error -- SDnametoindex failed on dataset %s in hdfread_vec\n", 
                dataset);
        return -1;
    }

    hdf4id_t sds_id = SDselect(sd_fid, sds_index);
    if (sds_id == FAIL) 
    {
        fprintf(stderr, 
                "error -- SDselect failed on dataset %s in hdfread_vec\n", 
                dataset);
        return -1;
    }
        
    // Read data
    int32 rank = 0;                    // number of array dimensions
    int32 dim_sizes[2] = {0, 0};    // dimension sizes are specified for 2d array
                                    // suitable for representing (1xN) or (Nx1) matrix
    int32 data_type = 0;            // HDF data type
    int32 n_attrs = 0;                // number of attributes
    char sds_name[60];                // string for data info
    
    int32 start[2] = {0, 0};        // data sizes for HDF read
    int32 edge[2] = {0, 0};
    
    int32 status = SDgetinfo(sds_id, sds_name, &rank, dim_sizes, &data_type, &n_attrs);
    if (status == FAIL) 
    {
        fprintf(stderr, 
                "error -- SDgetinfo failed on dataset %s in hdfread_vec\n", 
                dataset);
        return -1;
    }
    else if (rank != 1) 
    {
        fprintf(stderr, "error -- dataset %s has %d dimensions in hdfread_vec\n", 
                dataset, rank);
        return -1;
    }
    else if (dim_sizes[0] <= 0)
    {
        fprintf(stderr, "error -- dataset %s has a  size of zero in hdfread_vec\n", 
                dataset);
        return -1;
    }
    else if (data_type == DFNT_FLOAT64) 
    {
        // Read double data directly
        
        // Note: start[2] values remain {0, 0}
        // CHECK: switch edge[0], edge[1]
        edge[0] = dim_sizes[0]; 
        edge[1] = 1;
        
        // Vector size
        size_t size = (size_t) dim_sizes[0];
        
        // Initialize vector
        vec_alloc(vec, size);
        
        // Read data
        status = SDreaddata(sds_id, start, NULL, edge, (VOIDP)(vec->vals));
        if (status == FAIL) 
        {
            fprintf(stderr, 
                    "error -- SDreaddata (float64) failed on dataset %s in hdfread_vec\n", 
                    dataset);
            return -1;
        }
    }
    else if (data_type == DFNT_FLOAT32) 
    {
        // Read float data first, then convert values to double
        
        // Note: start[2] values remain {0, 0}
        // CHECK: switch edge[0], edge[1]
        edge[0] = dim_sizes[0]; 
        edge[1] = 1;
        
        // Vector size
        size_t size = (size_t) dim_sizes[0]; // CHECK: (size_t) dim_sizes[1];
        
        VecFloat32 temp; // temporary float32 vector
        vec_float32_init(&temp);
        vec_float32_alloc(&temp, size); // allocate float vector
        
        status = SDreaddata(sds_id, start, NULL, edge, (VOIDP)(temp.vals)); // read data
        if (status == FAIL) 
        {
            fprintf(stderr, 
                    "error -- SDreaddata (float32) failed on dataset %s in hdfread_vec\n", 
                    dataset);
            return -1;
        }
        else
        {
            vec_alloc(vec, size);            // initialize, allocate double vector
            vec_copy_float32(vec, &temp);   // copy values converting from float to double 
        }
        
        vec_float32_clear(&temp);    // clear, deallocate, temp vector
    }
    else
    {
        fprintf(stderr, 
                "error -- unrecognized data type on dataset %s in hdfread_vec\n", 
                dataset);
        return -1;
    }

    // Close dataset
    SDendaccess(sds_id);    // UPDATE: check for errors

    return 0;
}

// --------------------------------------------------------------
// HDF File I/O: Read HDF dataset containing 2d array
// --------------------------------------------------------------

// Read HDF 2d array to 2d double matrix
// 1. Assuming that file is open already
// 2. Assuming that the matrix is not initialized yet
// 3. Allocating and initializing double 2d matrix

int hdfread_mat(hdf4id_t sd_fid, const char *dataset, Matrix *mat)
{
    // Open dataset
    // Open dataset
    hdf4id_t sds_index = SDnametoindex(sd_fid, dataset);
    if (sds_index == FAIL)
    {
        fprintf(stderr, 
                "error -- SDnametoindex failed on dataset %s in hdfread_mat\n", 
                dataset);
        return -1;
    }

    hdf4id_t sds_id = SDselect(sd_fid, sds_index);
    if (sds_id == FAIL) 
    {
        fprintf(stderr, 
                "error -- SDselect failed on dataset %s in hdfread_mat\n", 
                dataset);
        return -1;
    }
        
    // Read data
    int32 rank = 0;                    // number of array dimensions
    int32 dim_sizes[2] = {0, 0};    // dimension sizes are specified for 2d array
                                    // suitable for representing (1xN) or (Nx1) matrix
    int32 data_type = 0;            // HDF data type
    int32 n_attrs = 0;                // number of attributes
    char sds_name[60];                // string for data info
    
    int32 start[2] = {0, 0};        // data sizes for HDF read
    int32 edge[2] = {0, 0};
    
    int32 status = SDgetinfo(sds_id, sds_name, &rank, dim_sizes, &data_type, &n_attrs);
    if (status == FAIL) 
    {
        fprintf(stderr, 
                "error -- SDgetinfo failed on dataset %s in hdfread_mat\n", 
                dataset);
        return -1;
    }
    else if (rank != 2) 
    {
        fprintf(stderr, "error -- dataset %s has %d dimensions in hdfread_mat\n", 
                dataset, rank);
        return -1;
    }
    else if (dim_sizes[0] <= 0 || dim_sizes[1] <= 0)
    {
        fprintf(stderr, "error -- dataset %s has a  size of zero in hdfread_mat\n", 
                dataset);
        return -1;
    }
    else if (data_type == DFNT_FLOAT64) 
    {
        // Read double data directly
        //
        // Note: start[2] values remain {0, 0}
        // CHECK: edge values = dimension sizes for reading all of the elements
        
        edge[0] = dim_sizes[0]; 
        edge[1] = dim_sizes[1];
        
        // Array sizes
        size_t size1 = (size_t) dim_sizes[0];
        size_t size2 = (size_t) dim_sizes[1];
        
        // Initialize matrix
        mat_alloc(mat, size1, size2);
        
        // Read data
        status = SDreaddata(sds_id, start, NULL, edge, (VOIDP)(mat->vals));
        if (status == FAIL) 
        {
            fprintf(stderr, 
                    "error -- SDreaddata (float64) failed on dataset %s in hdfread_mat\n", 
                    dataset);
        }
    }
    else if (data_type == DFNT_FLOAT32) 
    {
        // Read float data first, then convert values to double
        
        edge[0] = dim_sizes[0];     // Note: start[2] values remain {0, 0}
        edge[1] = dim_sizes[1];
        
        size_t size1 = (size_t) dim_sizes[0]; // array sizes
        size_t size2 = (size_t) dim_sizes[1];
        
        MatFloat32 temp; // temporary float32 matrix
        mat_float32_init(&temp);
        mat_float32_alloc(&temp, size1, size2);    // allocate float matrix
        
        status = SDreaddata(sds_id, start, NULL, edge, (VOIDP)(temp.vals)); // read data
        if (status == FAIL) 
        {
            fprintf(stderr, 
                    "error -- SDreaddata (float32) failed on dataset %s in hdfread_mat\n", 
                    dataset);
            return -1;
        }
        else
        {
            mat_alloc(mat, size1, size2);    // initialize, allocate double matrix
            mat_copy_float32(mat, &temp);   // copy values converting from float to double 
        }
        
        mat_float32_clear(&temp);    // clear, deallocate, temp matrix
    }
    else if (data_type == DFNT_UINT16 || data_type == DFNT_INT16) 
    {
        // Read int16 data first, then convert values to double
        
        edge[0] = dim_sizes[0];     // Note: start[2] values remain {0, 0}
        edge[1] = dim_sizes[1];
        
        size_t size1 = (size_t) dim_sizes[0]; // array sizes
        size_t size2 = (size_t) dim_sizes[1];
        
        MatInt16 temp;
        mat_int16_init(&temp);
        mat_int16_alloc(&temp, size1, size2);    // allocate int16 matrix

        status = SDreaddata(sds_id, start, NULL, edge, (VOIDP)(temp.vals)); // read data
        if (status == FAIL) 
        {
            fprintf(stderr, 
                    "error -- SDreaddata (int16) failed on dataset %s in hdfread_mat\n", 
                    dataset);
            mat_int16_clear(&temp);    // clear, deallocate, temp vector
            return -1;
        }
        else
        {
            mat_alloc(mat, size1, size2);  // initialize, allocate double matrix
            if (data_type == DFNT_UINT16)
                mat_copy_uint16(mat, (uint16*)temp.vals);
            else
                mat_copy_int16(mat, temp.vals);    // copy values converting from int16 to double 
        }
        
        mat_int16_clear(&temp);    // clear, deallocate, temp matrix
    }
    else
    {
        fprintf(stderr, 
                "error -- unrecognized data type on dataset %s in hdfread_mat\n", 
                dataset);
    }

    // Close dataset
    SDendaccess(sds_id);    // UPDATE: check for errors

    return 0;
}

// --------------------------------------------------------------
// HDF File I/O: Read HDF dataset containing 2d UINT8 array
// --------------------------------------------------------------

// Read HDF 2d array to 2d uint8 matrix

int hdfread_mat_uint8(hdf4id_t sd_fid, const char *dataset, MatUint8 *mat)
{
    // Open dataset
    hdf4id_t sds_index = SDnametoindex(sd_fid, dataset);
    if (sds_index == FAIL)
    {
        fprintf(stderr,
                "error -- SDnametoindex failed on dataset %s in hdfread_mat_uint8\n", 
                dataset);
        return -1;
    }

    hdf4id_t sds_id = SDselect(sd_fid, sds_index);
    if (sds_id == FAIL) 
    {
        fprintf(stderr, 
                "error -- SDselect failed on dataset %s in hdfread_mat_uint8\n", 
                dataset);
        return -1;
    }

    // Read data
    int32 rank = 0;                    // number of array dimensions
    int32 dim_sizes[2] = {0, 0};    // dimension sizes for 2d array
    int32 data_type = 0;            // HDF data type
    int32 n_attrs = 0;                // number of attributes
    char sds_name[60];                // string for data info
    int32 start[2] = {0, 0};        // data sizes for HDF read
    int32 edge[2] = {0, 0};
    
    int32 status = SDgetinfo(sds_id, sds_name, &rank, dim_sizes, &data_type, &n_attrs);
    if (status == FAIL) 
    {
        fprintf(stderr, 
                "error -- SDgetinfo failed on dataset %s in hdfread_mat_uint8\n", 
                dataset);
        return -1;
    }
    else if (rank != 2) 
    {
        fprintf(stderr, "error -- dataset %s has %d dimensions in hdfread_mat_uint8\n", 
                dataset, rank);
        return -1;
    }
    else if (dim_sizes[0] <= 0 || dim_sizes[1] <= 0)
    {
        fprintf(stderr, "error -- dataset %s has a  size of zero in hdfread_mat_uint8\n", 
                dataset);
        return -1;
    }
    else if (data_type == DFNT_UINT8 || data_type == DFNT_INT8) 
    {
        // Read uint8 array from the dataset
        // Note: start[2] values remain {0, 0}
        edge[0] = dim_sizes[0]; 
        edge[1] = dim_sizes[1];
        
        // Array sizes
        size_t size1 = (size_t) dim_sizes[0];
        size_t size2 = (size_t) dim_sizes[1];
        
        // Initialize matrix
        mat_uint8_alloc(mat, size1, size2);
        
        // Read data
        status = SDreaddata(sds_id, start, NULL, edge, (VOIDP)(mat->vals));
        if (status == FAIL) 
        {
            fprintf(stderr, 
                    "error -- SDreaddata (uint8) failed on dataset %s in hdfread_mat_uint8\n", 
                    dataset);
            return -1;
        }
    }
    else
    {
        fprintf(stderr, 
                "error -- data type not DFNT_UINT8 or DFNT_INT8 on dataset %s in hdfread_mat_uint8\n", 
                dataset);
        return -1;
    }

    // Close dataset
    SDendaccess(sds_id);    // UPDATE: check for errors

    return 0;
}

// --------------------------------------------------------------
// HDF File I/O: Read HDF FLOAT32 dataset setting 3d double matrix
// --------------------------------------------------------------

// Read HDF 3d array to 3d double matrix
// 1. Assuming that file is open already
// 2. Assuming that the matrix is not initialized yet
// 3. Allocating and initializing double 3d matrix

int hdfread_mat3d(hdf4id_t sd_fid, const char *dataset, Mat3d *mat)
{
    // Open dataset
    hdf4id_t sds_index = SDnametoindex(sd_fid, dataset);
    if (sds_index == FAIL)
    {
        fprintf(stderr,
                "error -- SDnametoindex failed on dataset %s in hdfread_mat3d\n", 
                dataset);
        return -1;
    }

    hdf4id_t sds_id = SDselect(sd_fid, sds_index);
    if (sds_id == FAIL) 
    {
        fprintf(stderr, 
                "error -- SDselect failed on dataset %s in hdfread_mat3d\n", 
                dataset);
        return -1;
    }
        
    // Read data
    int32 rank = 0;                 // number of array dimensions
    int32 dim_sizes[3] = {0, 0, 0}; // dimension sizes for 3d array
    int32 data_type = 0;            // HDF data type
    int32 n_attrs = 0;              // number of attributes
    char sds_name[60];              // string for data info
    int32 start[3] = {0, 0, 0};     // data sizes for HDF read
    int32 edge[3] = {0, 0, 0};
    
    int32 status = SDgetinfo(sds_id, sds_name, &rank, dim_sizes, &data_type, &n_attrs);
    if (status == FAIL) 
    {
        fprintf(stderr, 
                "error -- SDgetinfo failed on dataset %s in hdfread_mat3d\n", 
                dataset);
        return -1;
    }
    else if (rank != 3) 
    {
        fprintf(stderr, "error -- dataset %s has %d dimensions in hdfread_mat3d\n", 
                dataset, rank);
        return -1;
    }
    else if (dim_sizes[0] <= 0 || dim_sizes[1] <= 0 || dim_sizes[2] <= 0)
    {
        fprintf(stderr, "error -- dataset %s has a  size of zero in hdfread_mat3d\n", 
                dataset);
        return -1;
    }
    else if (data_type == DFNT_FLOAT32) 
    {
        // Read float data first, then convert values to double
        
        edge[0] = dim_sizes[0];     // Note: start values remain {0, 0, 0}
        edge[1] = dim_sizes[1];
        edge[2] = dim_sizes[2];
        
        size_t size1 = (size_t) dim_sizes[0]; // array sizes
        size_t size2 = (size_t) dim_sizes[1];
        size_t size3 = (size_t) dim_sizes[2];
        
        Mat3dFloat32 temp;
        mat3d_float32_init(&temp);
        mat3d_float32_alloc(&temp, size1, size2, size3);    // allocate float matrix
        
        status = SDreaddata(sds_id, start, NULL, edge, (VOIDP)(temp.vals)); // read data
        if (status == FAIL) 
        {
            fprintf(stderr, 
                    "error -- SDreaddata (float32) failed on dataset %s in hdfread_mat3d\n", 
                    dataset);
            return -1;
        }
        else
        {
            mat3d_alloc(mat, size1, size2, size3);    // initialize, allocate double matrix
            mat3d_copy_float32(mat, &temp);           // copy values converting from float to double 
        }
        
        mat3d_float32_clear(&temp);    // clear, deallocate, temp matrix
    }
    else if (data_type == DFNT_FLOAT64)
    {
        edge[0] = dim_sizes[0];     // Note: start values remain {0, 0, 0}
        edge[1] = dim_sizes[1];
        edge[2] = dim_sizes[2];

        size_t size1 = (size_t) dim_sizes[0]; // array sizes
        size_t size2 = (size_t) dim_sizes[1];
        size_t size3 = (size_t) dim_sizes[2];
        mat3d_alloc(mat, size1, size2, size3);
        status = SDreaddata(sds_id, start, NULL, edge, (VOIDP)(mat->vals)); // read data
        if (status == FAIL)
        {
            fprintf(stderr,
                    "error -- SDreaddata (float64) failed on dataset %s in hdfread_mat3d\n",
                    dataset);
            return -1;
        }
    }
    else
    {
        fprintf(stderr, 
                "error -- data type not DFNT_FLOAT32 ot DFNT_FLOAT64 on dataset %s in hdfread_mat3d\n", 
                dataset);
        return -1;
    }
    
    // Close dataset
    SDendaccess(sds_id);    // UPDATE: check for errors

    return 0;
}

// --------------------------------------------------------------
// HDF File I/O: Read HDF FLOAT32 dataset setting 4d double matrix
// --------------------------------------------------------------

// Read HDF 4d array to 4d double matrix
// 1. Assuming that file is open already
// 2. Assuming that the matrix is not initialized yet
// 3. Allocating and initializing double 4d matrix

int hdfread_mat4d(hdf4id_t sd_fid, const char *dataset, Mat4d *mat)
{
    // Open dataset
    hdf4id_t sds_index = SDnametoindex(sd_fid, dataset);
    if (sds_index == FAIL)
    {
        fprintf(stderr,
                "error -- SDnametoindex failed on dataset %s in hdfread_mat4d\n", 
                dataset);
        return -1;
    }

    hdf4id_t sds_id = SDselect(sd_fid, sds_index);
    if (sds_id == FAIL) 
    {
        fprintf(stderr, 
                "error -- SDselect failed on dataset %s in hdfread_mat4d\n", 
                dataset);
        return -1;
    }

    // Read data
    int32 rank = 0;                        // number of array dimensions
    int32 dim_sizes[4] = {0, 0, 0, 0};    // dimension sizes for 4d array
    int32 data_type = 0;                // HDF data type
    int32 n_attrs = 0;                    // number of attributes
    char sds_name[60];                    // string for data info
    
    int32 start[4] = {0, 0, 0, 0};        // data sizes for HDF read
    int32 edge[4] = {0, 0, 0, 0};
    
    int32 status = SDgetinfo(sds_id, sds_name, &rank, dim_sizes, &data_type, &n_attrs);
    if (status == FAIL) 
    {
        fprintf(stderr, 
                "error -- SDgetinfo failed on dataset %s in hdfread_mat4d\n", 
                dataset);
        return -1;
    }
    else if (rank != 4) 
    {
        fprintf(stderr, "error -- dataset %s has %d dimensions in hdfread_mat4d\n", 
                dataset, rank);
        return -1;
    }
    else if (dim_sizes[0] <= 0 || dim_sizes[1] <= 0 || dim_sizes[2] <= 0 || dim_sizes[3] <= 0)
    {
        fprintf(stderr, "error -- dataset %s has a  size of zero in hdfread_mat4d\n", 
                dataset);
        return -1;
    }
    else if (data_type == DFNT_FLOAT32) 
    {
        // Read float data first, then convert values to double
        
        edge[0] = dim_sizes[0];     // Note: start values remain {0, 0, 0, 0}
        edge[1] = dim_sizes[1];
        edge[2] = dim_sizes[2];
        edge[3] = dim_sizes[3];
        
        size_t size1 = (size_t) dim_sizes[0]; // array sizes
        size_t size2 = (size_t) dim_sizes[1];
        size_t size3 = (size_t) dim_sizes[2];
        size_t size4 = (size_t) dim_sizes[3];
        
        Mat4dFloat32 temp;
        mat4d_float32_init(&temp);
        mat4d_float32_alloc(&temp, size1, size2, size3, size4);    // allocate float matrix
        
        status = SDreaddata(sds_id, start, NULL, edge, (VOIDP)(temp.vals)); // read data
        
        if (status == FAIL) 
        {
            fprintf(stderr, 
                    "error -- SDreaddata (float32) failed on dataset %s in hdfread_mat4d\n", 
                    dataset);
            return -1;
        }
        else
        {
            mat4d_alloc(mat, size1, size2, size3, size4);    // initialize, allocate double matrix
            mat4d_copy_float32(mat, &temp);                   // copy values converting from float to double 
        }
        
        mat4d_float32_clear(&temp);    // clear, deallocate, temp matrix
    }
    else
    {
        fprintf(stderr, 
                "error -- data type not DFNT_FLOAT32 on dataset %s in hdfread_mat3d\n", 
                dataset);
        return -1;
    }
    
    // Close dataset
    SDendaccess(sds_id);    // UPDATE: check for errors

    return 0;
}

int scale_and_offset_mat(hdf4id_t sd_fid, const char *dataset, Matrix *mat)
{
    int status;
    hdf4id_t attr_index, attr_type, attr_len;
    char attr_name[FILEIO_MAX_ATTRNAME];
    float32 scale32, offset32, fillvalue32;
    double scale, offset;
    double fillvalue = REAL_NAN;

    // Get dataset index
    hdf4id_t sds_index = SDnametoindex(sd_fid, dataset);
    if (sds_index == FAIL) 
    {
        fprintf(stderr, "SDnametoindex sds_index failed for %s in scale_and_offset\n",
                dataset);
        return -1;
    }

    // Access the dataset
    hdf4id_t sds_id = SDselect(sd_fid, sds_index);
    if (sds_id == FAIL) 
    {
        fprintf(stderr, "SDselect sds_id failed for %s in scale_and_offset\n",
                dataset);
        return -1;
    }

    // Get attribute index for "Scale"
    attr_index = SDfindattr(sds_id, "Scale");
    if (attr_index == FAIL)
    {
        // Some are named "_Scale"
        attr_index = SDfindattr(sds_id, "_Scale");
    }
    if (attr_index == FAIL)
    {
        // Some are named "scale_factor"
        attr_index = SDfindattr(sds_id, "scale_factor");
    }

    // If still failed, assume scale = 1.0 and skip scaling
    if (attr_index == FAIL)
        scale = 1.0;
    else
    {
        // Get attribute4 info for "Scale"
        status = SDattrinfo(sds_id, attr_index, attr_name, &attr_type, &attr_len);
        if (status == FAIL)
        {
            fprintf(stderr, "SDattrinfo failed for \"Scale\" or "
                    "\"scale_factpr\" on %s in scale_and_offset\n",
                    dataset);
            return -1;
        }

        // Verify data type for scale
        if (attr_type == DFNT_FLOAT32)
        {
            status = SDreadattr(sds_id, attr_index, &scale32);
            scale = scale32;
        }
        else if (attr_type == DFNT_FLOAT64)
        {
            status = SDreadattr(sds_id, attr_index, &scale);
        }
        else
        {
            fprintf(stderr, "SDattrinfo unexpected datatype for \"Scale\" or "
                    "\"scale_factor\" on %s in scale_and_offset\n",
                    dataset);
            fprintf(stderr, "  expected DFNT_FLOAT64(%d) or DFNT_FLOAT32(%d) but got %d\n",
                    DFNT_FLOAT64, DFNT_FLOAT32, attr_type);
            return -1;
        }

        // Get the value of "Scale"
        if (status == FAIL)
        {
            fprintf(stderr, "SDreadattr failed for \"Scale\" or "
                    "\"scale_factor\" on %s in scale_and_offset\n",
                    dataset);
            return -1;
        }
    }

    // Get attribute index for "Offset"
    attr_index = SDfindattr(sds_id, "Offset");
    if (attr_index == FAIL)
    {
        offset = 0.0;
    }
    else
    {
        // Get attribute info for "Offset"
        status = SDattrinfo(sds_id, attr_index, attr_name, &attr_type, &attr_len);
        if (status == FAIL)
        {
            fprintf(stderr, "SDattrinfo failed for \"Offset\" on %s in scale_and_offset\n",
                    dataset);
            return -1;
        }

        // Verify data type for Offset
        if (attr_type == DFNT_FLOAT32)
        {
            status = SDreadattr(sds_id, attr_index, &offset32);
            offset = offset32;
        }
        else if (attr_type == DFNT_FLOAT64)
        {
            status = SDreadattr(sds_id, attr_index, &offset);
        }
        else
        {
            fprintf(stderr, "SDattrinfo unexpected datatype for \"Offset\" on %s in scale_and_offset\n",
                    dataset);
            fprintf(stderr, "  expected DFNT_FLOAT64(%d) or DFNT_FLOAT32(%d) but got %d\n",
                    DFNT_FLOAT64, DFNT_FLOAT32, attr_type);
            return -1;
        }

        // Get the value of "Offset"
        if (status == FAIL)
        {
            fprintf(stderr, "SDreadattr failed for \"Offset\" on %s in scale_and_offset\n",
                    dataset);
            return -1;
        }
    }

    // Get attribute index for "_FillValue"
    attr_index = SDfindattr(sds_id, "_FillValue");
    if (attr_index != FAIL)
    {
        // Get attribute info for "_FillValue"
        status = SDattrinfo(sds_id, attr_index, attr_name, &attr_type, &attr_len);
        if (status == FAIL)
        {
            fprintf(stderr, "SDattrinfo failed for \"_FillValue\" on %s "
                    "in scale_and_offset\n", dataset);
            return -1;
        }

        // Verify data type for _FillValue
        if (attr_type == DFNT_FLOAT32)
        {
            status = SDreadattr(sds_id, attr_index, &fillvalue32);
            fillvalue = fillvalue32;
        }
        else if (attr_type == DFNT_FLOAT64)
        {
            status = SDreadattr(sds_id, attr_index, &fillvalue);
        }
        else
        {
            fprintf(stderr, "SDattrinfo unexpected datatype for \"_FillValue\" on %s "
                    "in scale_and_offset\n", dataset);
            fprintf(stderr, "  expected DFNT_FLOAT64(%d) or DFNT_FLOAT32(%d) but got %d\n",
                    DFNT_FLOAT64, DFNT_FLOAT32, attr_type);
            return -1;
        }

        // Get the value of "_FillValue"
        if (status == FAIL)
        {
            fprintf(stderr, "SDreadattr failed for \"_FillValue\" on %s in scale_and_offset\n",
                    dataset);
            return -1;
        }
    }

    // Done with dataset
    SDendaccess(sds_id);

    // Change FillValues to NaN, then apply the Scale and Offset to the matrix.
    int i;
    for (i = 0; i < mat->size1 * mat->size2; i++)
    {
        if (dequal(mat->vals[i], fillvalue))
        {
            mat->vals[i] = REAL_NAN;
        }
        else
        {
        mat->vals[i] *= scale;
        mat->vals[i] += offset;
        }
    }
#ifdef COMMENT
    printf("scale_and_offset_mat: %s scaled by %.3f and offset by %.2f\n", 
           dataset, scale, offset);
#endif
    return 0;
}

// Write HDF 1D array from double vector
int hdfwrite_vec(hdf4id_t sd_fid, const char *dataset, Vector *vec)
{
    assert(sd_fid != -1);
    hdf4id_t dhdsid, dh_dim_sizes[2], dh_start[2], dh_edge[2];
    dh_dim_sizes[0] = vec->size;
    dh_dim_sizes[1] = 1;
    dhdsid = SDcreate(sd_fid, dataset, DFNT_FLOAT64, 1, dh_dim_sizes);
    if(dhdsid == FAIL) 
    {
        fprintf(stderr, "hdfwrite_vec: SDcreate failed for %s\n",
            dataset);
        return -1;
    }
    dh_start[0] = 0;
    dh_start[1] = 0;
    dh_edge[0] = vec->size;
    dh_edge[1] = 1;
    int stat;

    stat = SDwritedata(dhdsid, dh_start, NULL, dh_edge, vec->vals);
    if (stat == FAIL)
    {
        fprintf(stderr, "hdfwrite_vec: SDwritedata failed on %s\n",
                dataset);
        return -1;
    }

    SDendaccess(dhdsid);
    return 0;
}

// Write HDF 2D data from double matrix
int hdfwrite_mat(hdf4id_t sd_fid, const char *dataset, Matrix *mat)
{
    assert(sd_fid != -1);
    hdf4id_t dhdsid, dh_dim_sizes[2], dh_start[2], dh_edge[2];
    dh_dim_sizes[0] = mat->size1;
    dh_dim_sizes[1] = mat->size2;
    dhdsid = SDcreate(sd_fid, dataset, DFNT_FLOAT64, 2, dh_dim_sizes);
    if(dhdsid == FAIL) 
    {
        fprintf(stderr, "hdfwrite_mat: SDcreate failed for %s\n",
            dataset);
        return -1;
    }
    dh_start[0] = 0;
    dh_start[1] = 0;
    dh_edge[0] = mat->size1;
    dh_edge[1] = mat->size2;
    int stat;

    stat = SDwritedata(dhdsid, dh_start, NULL, dh_edge, mat->vals);
    if (stat == FAIL)
    {
        fprintf(stderr, "hdfwrite_mat: SDwritedata failed on %s\n",
                dataset);
        return -1;
    }

    SDendaccess(dhdsid);
    return 0;
}

// Write HDF 2D data from uint8 matrix
int hdfwrite_mat_uint8(hdf4id_t sd_fid, const char *dataset, MatUint8 *mat)
{
    assert(sd_fid != -1);
    hdf4id_t dhdsid, dh_dim_sizes[2], dh_start[2], dh_edge[2];
    dh_dim_sizes[0] = mat->size1;
    dh_dim_sizes[1] = mat->size2;
    dhdsid = SDcreate(sd_fid, dataset, DFNT_UINT8, 2, dh_dim_sizes);
    if(dhdsid == FAIL) 
    {
        fprintf(stderr, "hdfwrite_mat_uint8: SDcreate failed for %s\n",
            dataset);
        return -1;
    }
    dh_start[0] = 0;
    dh_start[1] = 0;
    dh_edge[0] = mat->size1;
    dh_edge[1] = mat->size2;
    int stat;

    stat = SDwritedata(dhdsid, dh_start, NULL, dh_edge, mat->vals);
    if (stat == FAIL)
    {
        fprintf(stderr, "hdfwrite_mat_uint8: SDwritedata failed on %s\n",
                dataset);
        return -1;
    }

    SDendaccess(dhdsid);
    return 0;
}

// Write HDF 2D data from uint16 matrix
int hdfwrite_mat_uint16(hdf4id_t sd_fid, const char *dataset, MatUint16 *mat)
{
    assert(sd_fid != -1);
    hdf4id_t dhdsid, dh_dim_sizes[2], dh_start[2], dh_edge[2];
    dh_dim_sizes[0] = mat->size1;
    dh_dim_sizes[1] = mat->size2;
    dhdsid = SDcreate(sd_fid, dataset, DFNT_UINT16, 2, dh_dim_sizes);
    if(dhdsid == FAIL) 
    {
        fprintf(stderr, "hdfwrite_mat_uint16: SDcreate failed for %s\n",
            dataset);
        return -1;
    }
    dh_start[0] = 0;
    dh_start[1] = 0;
    dh_edge[0] = mat->size1;
    dh_edge[1] = mat->size2;
    int stat;

    stat = SDwritedata(dhdsid, dh_start, NULL, dh_edge, mat->vals);
    if (stat == FAIL)
    {
        fprintf(stderr, "hdfwrite_mat_uint16: SDwritedata failed on %s\n",
                dataset);
        return -1;
    }

    SDendaccess(dhdsid);
    return 0;
}

// Write HDF 3D data from double matrix
int hdfwrite_mat3d(hdf4id_t sd_fid, const char *dataset, Mat3d *mat)
{
    assert(sd_fid != -1);
    hdf4id_t dhdsid, dh_dim_sizes[3], dh_start[3], dh_edge[3];
    dh_dim_sizes[0] = mat->size1;
    dh_dim_sizes[1] = mat->size2;
    dh_dim_sizes[2] = mat->size3;
    dhdsid = SDcreate(sd_fid, dataset, DFNT_FLOAT64, 3, dh_dim_sizes);
    if(dhdsid == FAIL) 
    {
        fprintf(stderr, "hdfwrite_mat3d: SDcreate failed for %s\n",
            dataset);
        return -1;
    }
    dh_start[0] = 0;
    dh_start[1] = 0;
    dh_start[2] = 0;
    dh_edge[0] = mat->size1;
    dh_edge[1] = mat->size2;
    dh_edge[2] = mat->size3;
    int stat;

    stat = SDwritedata(dhdsid, dh_start, NULL, dh_edge, mat->vals);
    if (stat == FAIL)
    {
        fprintf(stderr, "hdfwrite_mat3d: SDwritedata failed on %s\n",
                dataset);
        return -1;
    }

    SDendaccess(dhdsid);
    return 0;
}

// Write HDF 4D data from double matrix
int hdfwrite_mat4d(hdf4id_t sd_fid, const char *dataset, Mat4d *mat)
{
    assert(sd_fid != -1);
    hdf4id_t dhdsid, dh_dim_sizes[4], dh_start[4], dh_edge[4];
    dh_dim_sizes[0] = mat->size1;
    dh_dim_sizes[1] = mat->size2;
    dh_dim_sizes[2] = mat->size3;
    dh_dim_sizes[3] = mat->size4;
    dhdsid = SDcreate(sd_fid, dataset, DFNT_FLOAT64, 4, dh_dim_sizes);
    if(dhdsid == FAIL) 
    {
        fprintf(stderr, "hdfwrite_mat4d: SDcreate failed for %s\n",
            dataset);
        return -1;
    }
    dh_start[0] = 0;
    dh_start[1] = 0;
    dh_start[2] = 0;
    dh_start[3] = 0;
    dh_edge[0] = mat->size1;
    dh_edge[1] = mat->size2;
    dh_edge[2] = mat->size3;
    dh_edge[3] = mat->size4;
    int stat;

    stat = SDwritedata(dhdsid, dh_start, NULL, dh_edge, mat->vals);
    if (stat == FAIL)
    {
        fprintf(stderr, "hdfwrite_mat4d: SDwritedata failed on %s\n",
                dataset);
        return -1;
    }

    SDendaccess(dhdsid);
    return 0;
}

// --------------------------------------------------------------
// HDF5 File I/O
// --------------------------------------------------------------
static bool hdf5_2d_compress = false;
static hsize_t hdf5_2d_chunk[2] = {1, 1};
static uint hdf5_2d_deflate = 0;

void set_hdf5_2d_compression_on(hsize_t dim1, hsize_t dim2, uint deflate)
{
    hdf5_2d_compress = true;
    hdf5_2d_chunk[0] = dim1;
    hdf5_2d_chunk[1] = dim2;
    hdf5_2d_deflate = MIN(deflate, 9);
}

void set_hdf5_2d_compression_off()
{
    hdf5_2d_compress = false;
}

// --------------------------------------------------------------
// HDF5 File I/O: Read HDF5 dataset containing 1d array
// --------------------------------------------------------------
hid_t open_hdf5(const char *fpname, unsigned nc_attrs)
{
    hid_t fid;

#ifdef COMMENT
    printf("Opening hdf5 file %s\n", fpname);
    fflush(stdout);
#endif

    fid = H5Fopen(fpname, nc_attrs, H5P_DEFAULT);

#ifdef COMMENT
    printf("Open file id=%d\n", fid);
    fflush(stdout);
#endif

    if (fid < 0)
    {
        fprintf(stderr, "H5Fopen failed for %s in open_hdf5\n", 
            fpname);
        return -1;
    }
    return fid;
}

hid_t create_hdf5(const char *fpname)
{
    hid_t fid = H5Fcreate(fpname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (fid < 0)
    {
        fprintf(stderr, "H5Fcreate failed for %s in create_hdf5\n", 
            fpname);
        return -1;
    }
    return fid;
}

void close_hdf5(int fid)
{
    H5Fclose(fid);
}

int hdf5read_vec(hid_t h5fid, const char *dataset, Vector *vec)
{
    hid_t dsetid;
    herr_t stat;
    
    // Open dataset
    dsetid = H5Dopen(h5fid, dataset H5D_OPEN);
    if (dsetid < 0)
    {
        //fprintf(stderr, "H5Dopen failed for %s in id=%d in hdf5read_vec\n", 
        //    dataset, h5fid);
        return -1;
    }

    // Read data
    //int dimid[2];
    hsize_t dim_sizes[2] = {0, 0};      // dimension sizes are specified for 2d array
    //nc_type data_type = 0;              // NetCDF data type
    //float32 fill_value_f32;         // fill value for float
    double fill_value;
    int i;
    //int n_attrs = 0;                // number of attributes
    //char nc_var_name[TES_STR_SIZE]; // string for data info
    //int rank = 0;                   // number of array dimensions
    hsize_t ndims;

    hid_t dspcid = H5Dget_space(dsetid);
    assert(dspcid >= 0);
    ndims = H5Sget_simple_extent_dims(dspcid, dim_sizes, NULL);
    assert(ndims == 1);
    assert(dim_sizes[0] > 0);
    
    // Get datatype
    hid_t dtypeid = H5Dget_type(dsetid);
    assert(dtypeid >= 0);
    H5T_class_t dtype = H5Tget_class(dtypeid);
    size_t dsize = H5Tget_size(dtypeid);
    H5Tclose(dtypeid);

    if (dtype == H5T_FLOAT && dsize > 4) 
    {
        // Vector size
        size_t size = (size_t) dim_sizes[0];

        // Initialize vector
        vec_alloc(vec, size);
        
        // Read double data directly
        //stat = nc_get_var_double(h5fid, dsetid, vec->vals);
        stat = H5Dread(dsetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec->vals);
        if (stat < 0) 
        {
            fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s -- hdf5read_vec\n", 
                dataset, h5fid, "DOUBLE");
            return -1;
        }
    }
    else if (dtype == H5T_FLOAT) 
    {
        // Read float data first, then convert values to double
        
        // Vector size
        size_t size = (size_t) dim_sizes[0]; // CHECK: (size_t) dim_sizes[1];
        
        VecFloat32 temp;
        vec_float32_init(&temp);
        vec_float32_alloc(&temp, size);        // allocate float vector
        
        // TBD -- H5 can read float data into double array based on arg 2.
        stat = H5Dread(dsetid, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp.vals);
        if (stat < 0) 
        {
            fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s -- hdf5read_vec\n", 
                dataset, h5fid, "FLOAT");
            return -1;
        }
        else
        {
            vec_alloc(vec, size);            // initialize, allocate double vector
            vec_copy_float32(vec, &temp);   // copy values converting from float to double 
        }
        
        vec_float32_clear(&temp);    // clear, deallocate, temp vector
    }
    else
    {
        fprintf(stderr, "unrecognised HDF5 data type for %s in h5fid=%d -- hdf5read_vec",
            dataset, h5fid);
        return -1;
    }

    // Get _FillValue attribute
    fill_value = DEFAULT_FILL_VALUE;
    if (H5Aexists_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT))
    {
        hid_t atid_fv = H5Aopen_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT, H5P_DEFAULT);
        if (atid_fv >= 0)
        {
            stat = H5Aread(atid_fv, H5T_NATIVE_DOUBLE, &fill_value);
            assert(stat >= 0);
            H5Aclose(atid_fv);
        }
    }

    // Convert fills to NaN
#ifdef DEBUG
    int nnans = 0;
#endif
    for (i = 0; i < vec->size; i++)
    {
        if (vec->vals[i] == fill_value)
        {
            vec->vals[i] = FILEIO_NAN;
#ifdef DEBUG
            nnans++;
#endif
        }
    }

#ifdef DEBUG
    if (nnans > 0)
        printf("%s contains %d NaN values.\n", dataset, nnans);
#endif

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    return 0;
}

int hdf5read_mat(int h5fid, const char *dataset, Matrix *mat)
{
    herr_t stat;
    hid_t dsetid;

    // Open dataset
    dsetid = H5Dopen(h5fid, dataset H5D_OPEN);
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Dopen failed for %s in id=%d -- hdf5read_mat",
            dataset, h5fid);
        return -1;
    }

    // Read data
    hsize_t dim_sizes[2] = {0, 0};    // dimension sizes are specified for 2d array
    double fill_value;
    int i;
    hsize_t ndims = 0;                // number of array dimensions
    
    hid_t dspcid = H5Dget_space(dsetid);
    assert(dspcid >= 0);
    ndims = H5Sget_simple_extent_dims(dspcid, dim_sizes, NULL);
    assert(ndims == 2);
    assert(dim_sizes[0] > 0);
    assert(dim_sizes[1] > 0);

    // Get datatype
    hid_t dtypeid = H5Dget_type(dsetid);
    assert(dtypeid >= 0);
    H5T_class_t dtype = H5Tget_class(dtypeid);
    H5T_sign_t dsign = H5Tget_sign(dtypeid);
    size_t dsize = H5Tget_size(dtypeid);

    // TBD Should be able to do away with reading in float32 and converting
    // to double simply by using H5T_NATIVE_DOUBLE for the memory type
    // for either float32 or float64.
    if (dtype == H5T_FLOAT && dsize > 4) 
    {
        // Array sizes
        size_t size1 = (size_t) dim_sizes[0];
        size_t size2 = (size_t) dim_sizes[1];
        
        // Initialize matrix
        mat_alloc(mat, size1, size2);

        // Read double data directly
        stat = H5Dread(dsetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat->vals);
        if (stat < 0)
        {
            fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s -- hdf5read_mat\n", 
                dataset, h5fid, "DOUBLE");
            return -1;
        }
    }
    else if (dtype == H5T_FLOAT) 
    {
        // Read float data first, then convert values to double
        
        // Vector size
        // Array sizes
        size_t size1 = (size_t) dim_sizes[0];
        size_t size2 = (size_t) dim_sizes[1];
        
        MatFloat32 temp;
        mat_float32_init(&temp);
        mat_float32_alloc(&temp, size1, size2);    // allocate float vector
        
        // TBD -- H5 can read float data into double array based on arg 2.
        stat = H5Dread(dsetid, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp.vals);
        if (stat < 0) 
        {
            fprintf(stderr, "H5Dread failed for %s in h5fid=%d -- hdf5read_mat\n", 
                dataset, h5fid);
            return -1;
        }
        else
        {
            mat_alloc(mat, size1, size2);     // initialize, allocate double vector
            mat_copy_float32(mat, &temp);   // copy values converting from float to double 
        }
        
        mat_float32_clear(&temp);    // clear, deallocate, temp vector
    }
    else if (dtype == H5T_INTEGER) 
    {
        if (dsign == H5T_SGN_NONE)
        {
            // Read uint16 data first, then convert values to double
            size_t size1 = (size_t) dim_sizes[0]; // array sizes
            size_t size2 = (size_t) dim_sizes[1];
            
            MatUint16 temp;
            mat_uint16_init(&temp);
            mat_uint16_alloc(&temp, size1, size2);    // allocate int16 matrix
            
            stat = H5Dread(dsetid, H5T_NATIVE_UINT16, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp.vals);
            if (stat < 0) 
            {
                fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s -- hdf5read_mat\n", 
                    dataset, h5fid, "UINT16");
                mat_uint16_clear(&temp);    // clear, deallocate, temp vector
                return -1;
            }
            else
            {
                mat_alloc(mat, size1, size2); // initialize, allocate double matrix
                mat_copy_uint16(mat, temp.vals);
            }
            mat_uint16_clear(&temp);    // clear, deallocate, temp vector
        }
        else
        {
            // Read uint16 data first, then convert values to double
            size_t size1 = (size_t) dim_sizes[0]; // array sizes
            size_t size2 = (size_t) dim_sizes[1];
            
            MatInt16 temp;
            mat_int16_init(&temp);
            mat_int16_alloc(&temp, size1, size2);    // allocate int16 matrix
            
            stat = H5Dread(dsetid, H5T_NATIVE_INT16, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp.vals);
            if (stat < 0) 
            {
                fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s -- hdf5read_mat\n", 
                    dataset, h5fid, "INT16");
                mat_int16_clear(&temp);    // clear, deallocate, temp vector
                return -1;
            }
            else
            {
                mat_alloc(mat, size1, size2); // initialize, allocate double matrix
                mat_copy_int16(mat, temp.vals);
            }
            mat_int16_clear(&temp);    // clear, deallocate, temp vector
        }
    }
    else
    {
        fprintf(stderr, "unrecognised HDF5 data type %d for %s in h5fid=%d "
            "-- hdf5read_mat\n",
            dtype, dataset, h5fid);
        return -1;
    }

    // Get _FillValue attribute
    fill_value = DEFAULT_FILL_VALUE;
    if (H5Aexists_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT))
    {
        hid_t atid_fv = H5Aopen_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT, H5P_DEFAULT);
        if (atid_fv >= 0)
        {
            stat = H5Aread(atid_fv, H5T_NATIVE_DOUBLE, &fill_value);
            assert(stat >= 0);
            H5Aclose(atid_fv);
        }
    }

    // Convert fills to NaN
#ifdef DEBUG
    int nnans = 0;
#endif

    for (i = 0; i < mat->size1 * mat->size2; i++)
    {
        if (mat->vals[i] == fill_value)
        {
            mat->vals[i] = FILEIO_NAN;
#ifdef DEBUG
            nnans++;
#endif
        }
    }

#ifdef DEBUG
    if (nnans > 0)
        printf("%s contains %d NaN values.\n", dataset, nnans);
#endif

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    return 0;
}

int hdf5read_mat_uint8(int h5fid, const char *dataset, MatUint8 *mat)
{
    herr_t stat;
    hid_t dsetid;

    // Open dataset
    dsetid = H5Dopen(h5fid, dataset H5D_OPEN);
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Dopen failed for %s in id=%d -- hdf5read_mat_int16",
            dataset, h5fid);
        return -1;
    }

    // Read data
    hsize_t dim_sizes[2] = {0, 0};    // dimension sizes are specified for 2d array
    hsize_t ndims = 0;                // number of array dimensions

    hid_t dspcid = H5Dget_space(dsetid);
    assert(dspcid >= 0);
    ndims = H5Sget_simple_extent_dims(dspcid, dim_sizes, NULL);
    assert(ndims == 2);
    assert(dim_sizes[0] > 0);
    assert(dim_sizes[1] > 0);

    // Get datatype
    hid_t dtypeid = H5Dget_type(dsetid);
    assert(dtypeid >= 0);
    H5T_class_t dtype = H5Tget_class(dtypeid);
    //size_t dsize = H5Tget_size(dtypeid);

    if (dtype == H5T_INTEGER) 
    {
        hid_t memtype = H5T_NATIVE_UINT8;
        size_t size1 = (size_t) dim_sizes[0]; // array sizes
        size_t size2 = (size_t) dim_sizes[1];
        mat_uint8_alloc(mat, size1, size2); // initialize, allocate double matrix

        stat = H5Dread(dsetid, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                       mat->vals);
        if (stat < 0) 
        {
            fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s -- "
                "hdf5read_mat_uint8\n", dataset, h5fid, "INT16");
            return -1;
        }
    }
    else
    {
        fprintf(stderr, "unrecognised HDF5 data type %d for %s in h5fid=%d "
            "-- hdf5read_mat_uint8\n",
            dtype, dataset, h5fid);
        return -1;
    }

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    return 0;
}

int hdf5read_mat_int16(int h5fid, const char *dataset, MatInt16 *mat)
{
    herr_t stat;
    hid_t dsetid;

    // Open dataset
    dsetid = H5Dopen(h5fid, dataset H5D_OPEN);
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Dopen failed for %s in id=%d -- hdf5read_mat_int16",
            dataset, h5fid);
        return -1;
    }

    // Read data
    hsize_t dim_sizes[2] = {0, 0};    // dimension sizes are specified for 2d array
    double fill_value;
    int i;
    hsize_t ndims = 0;                // number of array dimensions

    hid_t dspcid = H5Dget_space(dsetid);
    assert(dspcid >= 0);
    ndims = H5Sget_simple_extent_dims(dspcid, dim_sizes, NULL);
    assert(ndims == 2);
    assert(dim_sizes[0] > 0);
    assert(dim_sizes[1] > 0);

    // Get datatype
    hid_t dtypeid = H5Dget_type(dsetid);
    assert(dtypeid >= 0);
    H5T_class_t dtype = H5Tget_class(dtypeid);
    //size_t dsize = H5Tget_size(dtypeid);

    if (dtype == H5T_INTEGER) 
    {
        size_t size1 = (size_t) dim_sizes[0]; // array sizes
        size_t size2 = (size_t) dim_sizes[1];
        mat_int16_alloc(mat, size1, size2); // initialize, allocate double matrix

        stat = H5Dread(dsetid, H5T_NATIVE_INT16, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                       mat->vals);
        if (stat < 0) 
        {
            fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s -- "
                "hdf5read_mat_int16\n", dataset, h5fid, "INT16");
            return -1;
        }
    }
    else
    {
        fprintf(stderr, "unrecognised HDF5 data type %d for %s in h5fid=%d "
            "-- hdf5read_mat_int16\n",
            dtype, dataset, h5fid);
        return -1;
    }

    // Get _FillValue attribute
    fill_value = DEFAULT_FILL_VALUE;
    if (H5Aexists_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT))
    {
        hid_t atid_fv = H5Aopen_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT, H5P_DEFAULT);
        if (atid_fv >= 0)
        {
            stat = H5Aread(atid_fv, H5T_NATIVE_DOUBLE, &fill_value);
            assert(stat >= 0);
            H5Aclose(atid_fv);
        }
    }

    // Convert fills to NaN
#ifdef DEBUG
    int nnans = 0;
#endif

    for (i = 0; i < mat->size1 * mat->size2; i++)
    {
        if (mat->vals[i] == fill_value)
        {
            mat->vals[i] = FILEIO_NAN;
#ifdef DEBUG
            nnans++;
#endif
        }
    }

#ifdef DEBUG
    if (nnans > 0)
        printf("%s contains %d NaN values.\n", dataset, nnans);
#endif

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    return 0;
}

int hdf5read_mat_int32(int h5fid, const char *dataset, MatInt32 *mat)
{
    herr_t stat;
    hid_t dsetid;

    // Open dataset
    dsetid = H5Dopen(h5fid, dataset H5D_OPEN);
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Dopen failed for %s in id=%d -- hdf5read_mat_int32\n",
            dataset, h5fid);
        return -1;
    }

    // Read data
    hsize_t dim_sizes[2] = {0, 0};    // dimension sizes are specified for 2d array
    double fill_value;
    int i;
    hsize_t ndims = 0;                // number of array dimensions

    hid_t dspcid = H5Dget_space(dsetid);
    assert(dspcid >= 0);
    ndims = H5Sget_simple_extent_dims(dspcid, dim_sizes, NULL);
    assert(ndims == 2);
    assert(dim_sizes[0] > 0);
    assert(dim_sizes[1] > 0);

    // Get datatype
    hid_t dtypeid = H5Dget_type(dsetid);
    assert(dtypeid >= 0);
    //H5T_class_t dtype = H5Tget_class(dtypeid);
    //size_t dsize = H5Tget_size(dtypeid);

    size_t size1 = (size_t) dim_sizes[0]; // array sizes
    size_t size2 = (size_t) dim_sizes[1];
    mat_int32_alloc(mat, size1, size2); // initialize, allocate double matrix

    stat = H5Dread(dsetid, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                   mat->vals);
    if (stat < 0) 
    {
        fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s -- "
            "hdf5read_mat_int32\n", dataset, h5fid, "INT16");
        return -1;
    }

    // Get _FillValue attribute
    fill_value = DEFAULT_FILL_VALUE;
    if (H5Aexists_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT))
    {
        hid_t atid_fv = H5Aopen_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT, H5P_DEFAULT);
        if (atid_fv >= 0)
        {
            stat = H5Aread(atid_fv, H5T_NATIVE_DOUBLE, &fill_value);
            assert(stat >= 0);
            H5Aclose(atid_fv);
        }
    }

    // Convert fills to NaN
#ifdef DEBUG
    int nnans = 0;
#endif

    for (i = 0; i < mat->size1 * mat->size2; i++)
    {
        if (mat->vals[i] == fill_value)
        {
            mat->vals[i] = FILEIO_NAN;
#ifdef DEBUG
            nnans++;
#endif
        }
    }

#ifdef DEBUG
    if (nnans > 0)
        printf("%s contains %d NaN values.\n", dataset, nnans);
#endif

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    return 0;
}

int hdf5read_mat_uint16(int h5fid, const char *dataset, MatUint16 *mat)
{
    herr_t stat;
    hid_t dsetid;
    double fill_value;

    // Open dataset
    dsetid = H5Dopen(h5fid, dataset H5D_OPEN);
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Dopen failed for %s in id=%d -- hdf5read_mat_int16",
            dataset, h5fid);
        return -1;
    }

    // Read data
    hsize_t dim_sizes[2] = {0, 0};    // dimension sizes are specified for 2d array
    int i;
    hsize_t ndims = 0;                // number of array dimensions

    hid_t dspcid = H5Dget_space(dsetid);
    assert(dspcid >= 0);
    ndims = H5Sget_simple_extent_dims(dspcid, dim_sizes, NULL);
    assert(ndims == 2);
    assert(dim_sizes[0] > 0);
    assert(dim_sizes[1] > 0);

    // Get datatype
    hid_t dtypeid = H5Dget_type(dsetid);
    assert(dtypeid >= 0);
    H5T_class_t dtype = H5Tget_class(dtypeid);
    //size_t dsize = H5Tget_size(dtypeid);

    if (dtype == H5T_INTEGER) 
    {
        size_t size1 = (size_t) dim_sizes[0]; // array sizes
        size_t size2 = (size_t) dim_sizes[1];
        mat_uint16_alloc(mat, size1, size2); // initialize, allocate double matrix

        stat = H5Dread(dsetid, H5T_NATIVE_UINT16, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                       mat->vals);
        if (stat < 0) 
        {
            fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s -- "
                "hdf5read_mat_uint16\n", dataset, h5fid, "INT16");
            return -1;
        }
    }
    else
    {
        fprintf(stderr, "unrecognised HDF5 data type %d for %s in h5fid=%d "
            "-- hdf5read_mat_uint16\n",
            dtype, dataset, h5fid);
        return -1;
    }

    // Get _FillValue attribute
    fill_value = DEFAULT_FILL_VALUE;
    if (H5Aexists_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT))
    {
        hid_t atid_fv = H5Aopen_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT, H5P_DEFAULT);
        if (atid_fv >= 0)
        {
            stat = H5Aread(atid_fv, H5T_NATIVE_DOUBLE, &fill_value);
            assert(stat >= 0);
            H5Aclose(atid_fv);
        }
    }

    // Convert fills to NaN
#ifdef DEBUG
    int nnans = 0;
#endif

    for (i = 0; i < mat->size1 * mat->size2; i++)
    {
        if (mat->vals[i] == fill_value)
        {
            mat->vals[i] = FILEIO_NAN;
#ifdef DEBUG
            nnans++;
#endif
        }
    }

#ifdef DEBUG
    if (nnans > 0)
        printf("%s contains %d NaN values.\n", dataset, nnans);
#endif

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    return 0;
}

int hdf5read_mat3d(int h5fid, const char *dataset, Mat3d *mat)
{
    herr_t stat;
    hid_t dsetid;

    // Open dataset
    dsetid = H5Dopen(h5fid, dataset H5D_OPEN);
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Dopen failed for %s in id=%d -- hdf5read_mat3d\n",
            dataset, h5fid);
        return -1;
    }

    // Read data
    hsize_t dim_sizes[3] = {0, 0, 0}; // dimension sizes are specified for 2d array
    double fill_value;
    int i;
    hsize_t ndims = 0;                // number of array dimensions
    
    hid_t dspcid = H5Dget_space(dsetid);
    assert(dspcid >= 0);
    ndims = H5Sget_simple_extent_dims(dspcid, dim_sizes, NULL);
    assert(ndims == 3);
    assert(dim_sizes[0] > 0);
    assert(dim_sizes[1] > 0);
    assert(dim_sizes[2] > 0);

    // Get datatype
    hid_t dtypeid = H5Dget_type(dsetid);
    assert(dtypeid >= 0);
    H5T_class_t dtype = H5Tget_class(dtypeid);
    //H5T_sign_t dsign = H5Tget_sign(dtypeid);
    size_t dsize = H5Tget_size(dtypeid);

    if (dtype == H5T_FLOAT && dsize > 4) 
    {
        // Array sizes
        size_t size1 = (size_t) dim_sizes[0];
        size_t size2 = (size_t) dim_sizes[1];
        size_t size3 = (size_t) dim_sizes[2];
        
        // Initialize matrix
        mat3d_alloc(mat, size1, size2, size3);

        // Read double data directly
        stat = H5Dread(dsetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat->vals);
        if (stat < 0) 
        {
            fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s "
                "-- hdf5read_mat3d\n", 
                dataset, h5fid, "DOUBLE");
            return -1;
        }
    }
    else if (dtype == H5T_FLOAT) 
    {
        // Read float data first, then convert values to double

        // Array sizes
        size_t size1 = (size_t) dim_sizes[0];
        size_t size2 = (size_t) dim_sizes[1];
        size_t size3 = (size_t) dim_sizes[2];
        
        Mat3dFloat32 temp; // temporary float32 vector
        mat3d_float32_init(&temp);
        mat3d_float32_alloc(&temp, size1, size2, size3); // allocate float vector
        
        stat = H5Dread(dsetid, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp.vals);
        if (stat < 0) 
        {
            fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s "
                "-- hdf5read_mat3d\n", 
                dataset, h5fid, "FLOAT");
            return -1;
        }
        else
        {
            mat3d_alloc(mat, size1, size2, size3); // initialize, allocate double vector
            mat3d_copy_float32(mat, &temp);   // copy values converting from float to double 
        }
        
        mat3d_float32_clear(&temp);    // clear, deallocate, temp vector
    }
    else
    {
        fprintf(stderr, "unrecognised HDF5 data type %d for %s in h5fid=%d"
            " -- hdf5read_mat3d\n",
            dtype, dataset, h5fid);
        return -1;
    }

    // Get _FillValue attribute
    fill_value = DEFAULT_FILL_VALUE;
    if (H5Aexists_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT))
    {
        hid_t atid_fv = H5Aopen_by_name(h5fid, dataset, "_FillValue", 
            H5P_DEFAULT, H5P_DEFAULT);
        if (atid_fv >= 0)
        {
            stat = H5Aread(atid_fv, H5T_NATIVE_DOUBLE, &fill_value);
            assert(stat >= 0);
            H5Aclose(atid_fv);
        }
    }
 
    // Convert fills to NaN
#ifdef DEBUG
    int nnans = 0;
#endif
    for (i = 0; i < mat->size1 * mat->size2 * mat->size3; i++)
    {
        if (mat->vals[i] == fill_value)
        {
            mat->vals[i] = FILEIO_NAN;
#ifdef DEBUG
            nnans++;
#endif
        }
    }

#ifdef DEBUG
    if (nnans > 0)
        printf("%s contains %d NaN values.\n", dataset, nnans);
#endif

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    return 0;
}

int hdf5read_mat3d_int16(int h5fid, const char *dataset, Mat3dInt16 *mat)
{
    herr_t stat;
    hid_t dsetid;

    // Open dataset
    dsetid = H5Dopen(h5fid, dataset H5D_OPEN);
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Dopen failed for %s in id=%d -- hdf5read_mat3d_int16\n",
            dataset, h5fid);
        return -1;
    }

    // Read data
    hsize_t dim_sizes[3] = {0, 0, 0}; // dimension sizes are specified for 2d array
    double fill_value;
    int i;
    hsize_t ndims = 0;                // number of array dimensions
    
    hid_t dspcid = H5Dget_space(dsetid);
    assert(dspcid >= 0);
    ndims = H5Sget_simple_extent_dims(dspcid, dim_sizes, NULL);
    assert(ndims == 3);
    assert(dim_sizes[0] > 0);
    assert(dim_sizes[1] > 0);
    assert(dim_sizes[2] > 0);

    // Get datatype
    hid_t dtypeid = H5Dget_type(dsetid);
    assert(dtypeid >= 0);
    H5T_class_t dtype = H5Tget_class(dtypeid);
    //size_t dsize = H5Tget_size(dtypeid);

    if (dtype == H5T_INTEGER) 
    {
        // Array sizes
        size_t size1 = (size_t) dim_sizes[0];
        size_t size2 = (size_t) dim_sizes[1];
        size_t size3 = (size_t) dim_sizes[2];
        
        // Initialize matrix
        mat3d_int16_alloc(mat, size1, size2, size3);

        // Read integer data directly
        stat = H5Dread(dsetid, H5T_NATIVE_INT16, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat->vals);
        if (stat < 0) 
        {
            fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s "
                "-- hdf5read_mat3d_int16\n", 
                dataset, h5fid, "DOUBLE");
            return -1;
        }
    }
    else
    {
        fprintf(stderr, "unrecognised HDF5 data type %d for %s in h5fid=%d"
            " -- hdf5read_mat3d_int16\n",
            dtype, dataset, h5fid);
        return -1;
    }

    // Get _FillValue attribute
    fill_value = DEFAULT_FILL_VALUE;
    if (H5Aexists_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT))
    {
        hid_t atid_fv = H5Aopen_by_name(h5fid, dataset, "_FillValue", 
            H5P_DEFAULT, H5P_DEFAULT);
        if (atid_fv >= 0)
        {
            stat = H5Aread(atid_fv, H5T_NATIVE_DOUBLE, &fill_value);
            assert(stat >= 0);
            H5Aclose(atid_fv);
        }
    }
 
    // Convert fills to NaN
#ifdef DEBUG
    int nnans = 0;
#endif
    for (i = 0; i < mat->size1 * mat->size2 * mat->size3; i++)
    {
        if (mat->vals[i] == fill_value)
        {
            mat->vals[i] = FILEIO_NAN;
#ifdef DEBUG
            nnans++;
#endif
        }
    }

#ifdef DEBUG
    if (nnans > 0)
        printf("%s contains %d NaN values.\n", dataset, nnans);
#endif

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    return 0;
}

int hdf5read_mat4d(int h5fid, const char *dataset, Mat4d *mat)
{
    herr_t stat;
    hid_t dsetid;

    // Open dataset
    dsetid = H5Dopen(h5fid, dataset H5D_OPEN);
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Dopen failed for %s in id=%d -- hdf5read_mat4d\n",
            dataset, h5fid);
        return -1;
    }

    // Read data
    hsize_t dim_sizes[4] = {0, 0, 0, 0}; // dimension sizes are specified for 2d array
    double fill_value;
    int i;
    hsize_t ndims = 0;                // number of array dimensions
    
    hid_t dspcid = H5Dget_space(dsetid);
    assert(dspcid >= 0);
    ndims = H5Sget_simple_extent_dims(dspcid, dim_sizes, NULL);
    assert(ndims == 4);
    assert(dim_sizes[0] > 0);
    assert(dim_sizes[1] > 0);
    assert(dim_sizes[2] > 0);
    assert(dim_sizes[3] > 0);

    // Get datatype
    hid_t dtypeid = H5Dget_type(dsetid);
    assert(dtypeid >= 0);
    H5T_class_t dtype = H5Tget_class(dtypeid);
    size_t dsize = H5Tget_size(dtypeid);

    if (dtype == H5T_FLOAT && dsize > 4) 
    {
        // Array sizes
        size_t size1 = (size_t) dim_sizes[0];
        size_t size2 = (size_t) dim_sizes[1];
        size_t size3 = (size_t) dim_sizes[2];
        size_t size4 = (size_t) dim_sizes[3];
        
        // Initialize matrix
        mat4d_alloc(mat, size1, size2, size3, size4);

        // Read double data directly
        stat = H5Dread(dsetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat->vals);
        if (stat < 0) 
        {
            fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s "
                "-- hdf5read_mat4d\n", 
                dataset, h5fid, "DOUBLE");
            return -1;
        }
    }
    else if (dtype == H5T_FLOAT) 
    {
        // Read float data first, then convert values to double
        
        // Array sizes
        size_t size1 = (size_t) dim_sizes[0];
        size_t size2 = (size_t) dim_sizes[1];
        size_t size3 = (size_t) dim_sizes[2];
        size_t size4 = (size_t) dim_sizes[3];
        
        Mat4dFloat32 temp;
        mat4d_float32_init(&temp);
        mat4d_float32_alloc(&temp, size1, size2, size3, size4); // allocate float vector
        
        stat = H5Dread(dsetid, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp.vals);
        if (stat < 0) 
        {
            fprintf(stderr, "H5Dread failed for %s in h5fid=%d type %s"
                " -- hdf5read_mat4d\n", 
                dataset, h5fid, "DOUBLE");
            return -1;
        }
        else
        {
            mat4d_alloc(mat, size1, size2, size3, size4); // initialize, allocate double vector
            mat4d_copy_float32(mat, &temp);   // copy values converting from float to double 
        }
        
        mat4d_float32_clear(&temp);    // clear, deallocate, temp vector
    }
    else
    {
        fprintf(stderr, "unrecognised HDF5 data type %d  for %s in h5fid=%d"
            " -- hdf5read_mat4d\n",
            dtype, dataset, h5fid);
        return -1;
    }

    // Get _FillValue attribute
    fill_value = DEFAULT_FILL_VALUE;
    if (H5Aexists_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT))
    {
        hid_t atid_fv = H5Aopen_by_name(h5fid, dataset, "_FillValue", H5P_DEFAULT, H5P_DEFAULT);
        if (atid_fv >= 0)
        {
            stat = H5Aread(atid_fv, H5T_NATIVE_DOUBLE, &fill_value);
            assert(stat >= 0);
            H5Aclose(atid_fv);
        }
    }

    // Convert fills to NaN
#ifdef DEBUG
    int nnans = 0;
#endif
    for (i = 0; i < mat->size1 * mat->size2 * mat->size3 * mat->size4; i++)
    {
        if (mat->vals[i] == fill_value)
        {
            mat->vals[i] = FILEIO_NAN;
#ifdef DEBUG
            nnans++;
#endif
        }
    }

#ifdef DEBUG
    if (nnans > 0)
        printf("%s contains %d NaN values.\n", dataset, nnans);
#endif

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    return 0;
}

int hdf5write_vec(hid_t h5fid, const char *dataset, Vector *vec)
{
    hid_t dsetid;
    herr_t stat;
    hsize_t dim_sizes[2] = {vec->size, 1};
    
    // Create dataset space
    hid_t dspcid = H5Screate_simple(1, dim_sizes, NULL);
    if (dspcid < 0)
    {
        fprintf(stderr, "H5Screate_simple failed for %s in id=%d in hdf5write_vec\n", 
            dataset, h5fid);
        return -1;
    }

#ifdef OLD_HDF5
    dsetid = H5Dcreate1(h5fid, dataset, H5T_NATIVE_DOUBLE, dspcid, H5P_DEFAULT);
#else
    dsetid = H5Dcreate2(h5fid, dataset, H5T_NATIVE_DOUBLE, dspcid, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Screate failed for %s in id=%d in hdf5write_vec\n", 
            dataset, h5fid);
        return -1;
    }

    stat = H5Dwrite(dsetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                    vec->vals);
    if (stat < 0)
    {
        fprintf(stderr, "H5Dwrite failed for %s in id=%d in hdf5write_vec\n", 
            dataset, h5fid);
    }

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    return stat;
}

int hdf5write_mat(hid_t h5fid, const char *dataset, Matrix *mat)
{
    hid_t dsetid;
    herr_t stat;
    hsize_t dim_sizes[2] = {mat->size1, mat->size2};
    hid_t dsprop_id = H5P_DEFAULT;
    if (hdf5_2d_compress)
    {
        // To perform compression, create a Dataset Property list.
        dsprop_id = H5Pcreate(H5P_DATASET_CREATE);
        // Set the chunk dimensions and deflate level on the property list.
        H5Pset_chunk(dsprop_id, 2, hdf5_2d_chunk);
        H5Pset_deflate(dsprop_id, hdf5_2d_deflate);
    }
    
    // Create dataset space
    hid_t dspcid = H5Screate_simple(2, dim_sizes, NULL);
    if (dspcid < 0)
    {
        fprintf(stderr, "H5Screate_simple failed for %s in id=%d in hdf5write_mat\n", 
            dataset, h5fid);
        return -1;
    }

#ifdef OLD_HDF5
    dsetid = H5Dcreate1(h5fid, dataset, H5T_NATIVE_DOUBLE, dspcid, H5P_DEFAULT);
#else
    dsetid = H5Dcreate2(h5fid, dataset, H5T_NATIVE_DOUBLE, dspcid, 
                          H5P_DEFAULT, dsprop_id, H5P_DEFAULT);
#endif
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Screate failed for %s in id=%d in hdf5write_mat\n", 
            dataset, h5fid);
        return -1;
    }

    stat = H5Dwrite(dsetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                    mat->vals);
    if (stat < 0)
    {
        fprintf(stderr, "H5Dwrite failed for %s in id=%d in hdf5write_mat\n", 
            dataset, h5fid);
    }

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    if (hdf5_2d_compress)
        H5Pclose(dsprop_id);
    return stat;
}

int hdf5write_mat_asfloat32(hid_t h5fid, const char *dataset, Matrix *mat)
{
    hid_t dsetid;
    herr_t stat;
    hsize_t dim_sizes[2] = {mat->size1, mat->size2};
    
    // Create dataset space
    hid_t dspcid = H5Screate_simple(2, dim_sizes, NULL);
    if (dspcid < 0)
    {
        fprintf(stderr, "H5Screate_simple failed for %s in id=%d in hdf5write_mat\n", 
            dataset, h5fid);
        return -1;
    }

    hid_t dsprop_id = H5P_DEFAULT;
    if (hdf5_2d_compress)
    {
        // To perform compression, create a Dataset Property list.
        dsprop_id = H5Pcreate(H5P_DATASET_CREATE);
        // Set the chunk dimensions and deflate level on the property list.
        H5Pset_chunk(dsprop_id, 2, hdf5_2d_chunk);
        H5Pset_deflate(dsprop_id, hdf5_2d_deflate);
    }
    
#ifdef OLD_HDF5
    dsetid = H5Dcreate1(h5fid, dataset, H5T_NATIVE_FLOAT, dspcid, H5P_DEFAULT);
#else
    dsetid = H5Dcreate2(h5fid, dataset, H5T_NATIVE_FLOAT, dspcid, 
                          H5P_DEFAULT, dsprop_id, H5P_DEFAULT);
#endif
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Screate failed for %s in id=%d in hdf5write_mat\n", 
            dataset, h5fid);
        return -1;
    }

    stat = H5Dwrite(dsetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                    mat->vals);
    if (stat < 0)
    {
        fprintf(stderr, "H5Dwrite failed for %s in id=%d in hdf5write_mat\n", 
            dataset, h5fid);
    }

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    if (hdf5_2d_compress)
        H5Pclose(dsprop_id);
    return stat;
}

int hdf5write_mat_uint8(hid_t h5fid, const char *dataset, MatUint8 *mat)
{
    hid_t dsetid;
    herr_t stat;
    hsize_t dim_sizes[2] = {mat->size1, mat->size2};
    
    // Create dataset space
    hid_t dspcid = H5Screate_simple(2, dim_sizes, NULL);
    if (dspcid < 0)
    {
        fprintf(stderr, "H5Screate_simple failed for %s in id=%d in hdf5write_mat\n", 
            dataset, h5fid);
        return -1;
    }

    hid_t dsprop_id = H5P_DEFAULT;
    if (hdf5_2d_compress)
    {
        // To perform compression, create a Dataset Property list.
        dsprop_id = H5Pcreate(H5P_DATASET_CREATE);
        // Set the chunk dimensions and deflate level on the property list.
        H5Pset_chunk(dsprop_id, 2, hdf5_2d_chunk);
        H5Pset_deflate(dsprop_id, hdf5_2d_deflate);
    }
    
#ifdef OLD_HDF5
    dsetid = H5Dcreate1(h5fid, dataset, H5T_NATIVE_UINT8, dspcid, H5P_DEFAULT);
#else
    dsetid = H5Dcreate2(h5fid, dataset, H5T_NATIVE_UINT8, dspcid, 
                          H5P_DEFAULT, dsprop_id, H5P_DEFAULT);
#endif
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Screate failed for %s in id=%d in hdf5write_mat\n", 
            dataset, h5fid);
        return -1;
    }

    stat = H5Dwrite(dsetid, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                    mat->vals);
    if (stat < 0)
    {
        fprintf(stderr, "H5Dwrite failed for %s in id=%d in hdf5write_mat\n", 
            dataset, h5fid);
    }

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    if (hdf5_2d_compress)
        H5Pclose(dsprop_id);
    return stat;
}

int hdf5write_mat_int16(hid_t h5fid, const char *dataset, MatInt16 *mat)
{
    hid_t dsetid;
    herr_t stat;
    hsize_t dim_sizes[2] = {mat->size1, mat->size2};
    
    // Create dataset space
    hid_t dspcid = H5Screate_simple(2, dim_sizes, NULL);
    if (dspcid < 0)
    {
        fprintf(stderr, "H5Screate_simple failed for %s in id=%d in hdf5write_mat_int16\n", 
            dataset, h5fid);
        return -1;
    }

    hid_t dsprop_id = H5P_DEFAULT;
    if (hdf5_2d_compress)
    {
        // To perform compression, create a Dataset Property list.
        dsprop_id = H5Pcreate(H5P_DATASET_CREATE);
        // Set the chunk dimensions and deflate level on the property list.
        H5Pset_chunk(dsprop_id, 2, hdf5_2d_chunk);
        H5Pset_deflate(dsprop_id, hdf5_2d_deflate);
    }
    
#ifdef OLD_HDF5
    dsetid = H5Dcreate1(h5fid, dataset, H5T_NATIVE_INT16, dspcid, H5P_DEFAULT);
#else
    dsetid = H5Dcreate2(h5fid, dataset, H5T_NATIVE_INT16, dspcid, 
                          H5P_DEFAULT, dsprop_id, H5P_DEFAULT);
#endif
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Screate failed for %s in id=%d in hdf5write_mat_int16\n", 
            dataset, h5fid);
        return -1;
    }

    stat = H5Dwrite(dsetid, H5T_NATIVE_INT16, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                    mat->vals);
    if (stat < 0)
    {
        fprintf(stderr, "H5Dwrite failed for %s in id=%d in hdf5write_mat_int16\n", 
            dataset, h5fid);
    }

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    if (hdf5_2d_compress)
        H5Pclose(dsprop_id);
    return stat;
}

int hdf5write_mat_uint16(hid_t h5fid, const char *dataset, MatUint16 *mat)
{
    hid_t dsetid;
    herr_t stat;
    hsize_t dim_sizes[2] = {mat->size1, mat->size2};
    
    // Create dataset space
    hid_t dspcid = H5Screate_simple(2, dim_sizes, NULL);
    if (dspcid < 0)
    {
        fprintf(stderr, "H5Screate_simple failed for %s in id=%d in hdf5write_mat_uint16\n", 
            dataset, h5fid);
        return -1;
    }

    hid_t dsprop_id = H5P_DEFAULT;
    if (hdf5_2d_compress)
    {
        // To perform compression, create a Dataset Property list.
        dsprop_id = H5Pcreate(H5P_DATASET_CREATE);
        // Set the chunk dimensions and deflate level on the property list.
        H5Pset_chunk(dsprop_id, 2, hdf5_2d_chunk);
        H5Pset_deflate(dsprop_id, hdf5_2d_deflate);
    }
    
#ifdef OLD_HDF5
    dsetid = H5Dcreate1(h5fid, dataset, H5T_NATIVE_UINT16, dspcid, H5P_DEFAULT);
#else
    dsetid = H5Dcreate2(h5fid, dataset, H5T_NATIVE_UINT16, dspcid, 
                          H5P_DEFAULT, dsprop_id, H5P_DEFAULT);
#endif
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Screate failed for %s in id=%d in hdf5write_mat_uint16\n", 
            dataset, h5fid);
        return -1;
    }

    stat = H5Dwrite(dsetid, H5T_NATIVE_UINT16, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                    mat->vals);
    if (stat < 0)
    {
        fprintf(stderr, "H5Dwrite failed for %s in id=%d in hdf5write_mat_uint16\n", 
            dataset, h5fid);
    }

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    if (hdf5_2d_compress)
        H5Pclose(dsprop_id);
    return stat;
}

int hdf5write_mat3d(hid_t h5fid, const char *dataset, Mat3d *mat)
{
    hid_t dsetid;
    herr_t stat;
    hsize_t dim_sizes[3] = {mat->size1, mat->size2, mat->size3};
    
    // Create dataset space
    hid_t dspcid = H5Screate_simple(3, dim_sizes, NULL);
    if (dspcid < 0)
    {
        fprintf(stderr, "H5Screate_simple failed for %s in id=%d in hdf5write_mat3d\n", 
            dataset, h5fid);
        return -1;
    }

#ifdef OLD_HDF5
    dsetid = H5Dcreate1(h5fid, dataset, H5T_NATIVE_DOUBLE, dspcid, H5P_DEFAULT);
#else
    dsetid = H5Dcreate2(h5fid, dataset, H5T_NATIVE_DOUBLE, dspcid, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Screate failed for %s in id=%d in hdf5write_mat3d\n", 
            dataset, h5fid);
        return -1;
    }

    stat = H5Dwrite(dsetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                    mat->vals);
    if (stat < 0)
    {
        fprintf(stderr, "H5Dwrite failed for %s in id=%d in hdf5write_mat3d\n", 
            dataset, h5fid);
    }

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    return stat;
}

int hdf5write_mat4d(hid_t h5fid, const char *dataset, Mat4d *mat)
{
    hid_t dsetid;
    herr_t stat;
    hsize_t dim_sizes[4] = {mat->size1, mat->size2, mat->size3, mat->size4};
    
    // Create dataset space
    hid_t dspcid = H5Screate_simple(4, dim_sizes, NULL);
    if (dspcid < 0)
    {
        fprintf(stderr, "H5Screate_simple failed for %s in id=%d in hdf5write_mat4d\n", 
            dataset, h5fid);
        return -1;
    }

#ifdef OLD_HDF5
    dsetid = H5Dcreate1(h5fid, dataset, H5T_NATIVE_DOUBLE, dspcid, H5P_DEFAULT);
#else
    dsetid = H5Dcreate2(h5fid, dataset, H5T_NATIVE_DOUBLE, dspcid, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Screate failed for %s in id=%d in hdf5write_mat4d\n", 
            dataset, h5fid);
        return -1;
    }

    stat = H5Dwrite(dsetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                    mat->vals);
    if (stat < 0)
    {
        fprintf(stderr, "H5Dwrite failed for %s in id=%d in hdf5write_mat4d\n", 
            dataset, h5fid);
    }

    H5Sclose(dspcid);
    H5Dclose(dsetid);
    return stat;
}

int hdf5_scale_and_offset_mat(hid_t h5fid, const char *dataset, Matrix *mat)
{
    herr_t stat;
    hid_t dsetid;

    // Open dataset
    dsetid = H5Dopen(h5fid, dataset H5D_OPEN);
    if (dsetid < 0)
    {
        fprintf(stderr, "H5Dopen failed for %s in id=%d -- "
            "hdf5_scale_and_offset_mat",
            dataset, h5fid);
        return -1;
    }

    double add_offset = 0.0;
    if (H5Aexists_by_name(h5fid, dataset, "add_offest", H5P_DEFAULT))
    {
        hid_t atid_fv = H5Aopen_by_name(h5fid, dataset, "add_offset",
                H5P_DEFAULT, H5P_DEFAULT);
        if (atid_fv >= 0)
        {
            stat = H5Aread(atid_fv, H5T_NATIVE_DOUBLE, &add_offset);
            assert(stat >= 0);
            H5Aclose(atid_fv);
        }
    }

    double scale_factor = 1.0;
    if (H5Aexists_by_name(h5fid, dataset, "scale_factor", H5P_DEFAULT))
    {
        hid_t atid_fv = H5Aopen_by_name(h5fid, dataset, "scale_factor", H5P_DEFAULT, H5P_DEFAULT);
        if (atid_fv >= 0)
        {
            stat = H5Aread(atid_fv, H5T_NATIVE_DOUBLE, &scale_factor);
            assert(stat >= 0);
            if (scale_factor == 0.0) scale_factor = 1.0;
            H5Aclose(atid_fv);
        }
    }
    H5Dclose(dsetid);

    // No need to process "multiply by 1 and add zero."
    if (scale_factor == 1.0 && add_offset == 0.0)
    {
        // No need to process
        return 0;
    }
    
    int i;
    for (i = 0; i < mat->size1 * mat->size2; i++)
    {
        mat->vals[i] *= scale_factor;
        mat->vals[i] += add_offset;
    }
    return 0;
}

// --------------------------------------------------------------
// NetCDF File I/O: Uses HDF5 to read NetCDF files
// --------------------------------------------------------------
hid_t open_netcdf(const char *fpname, unsigned nc_attrs)
{
    return open_hdf5(fpname, nc_attrs);
}

hid_t create_netcdf(const char *fpname)
{
    return create_hdf5(fpname);
}

void close_netcdf(hid_t ncid)
{
    H5Fclose(ncid);
}

int netcdfread_vec(hid_t ncid, const char *dataset, Vector *vec)
{
    return hdf5read_vec(ncid, dataset, vec);
}

int netcdfread_mat(int ncid, const char *dataset, Matrix *mat)
{
    return hdf5read_mat(ncid, dataset, mat);
}

int netcdfread_mat3d(int ncid, const char *dataset, Mat3d *mat)
{
    return hdf5read_mat3d(ncid, dataset, mat);
}

int netcdfread_mat4d(int ncid, const char *dataset, Mat4d *mat)
{
    return hdf5read_mat4d(ncid, dataset, mat);
}

int netcdfwrite_vec(hid_t ncid, const char *dataset, Vector *vec)
{
    return hdf5write_vec(ncid, dataset, vec);
}

int netcdfwrite_mat(int ncid, const char *dataset, Matrix *mat)
{
    return hdf5write_mat(ncid, dataset, mat);
}

int netcdfwrite_mat3d(int ncid, const char *dataset, Mat3d *mat)
{
    return hdf5write_mat3d(ncid, dataset, mat);
}

int netcdfwrite_mat4d(int ncid, const char *dataset, Mat4d *mat)
{
    return hdf5write_mat4d(ncid, dataset, mat);
}
