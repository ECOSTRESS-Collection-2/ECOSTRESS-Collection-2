#pragma once
#include <stdint.h>

/// @file Morphology Functions
///
// References: 
/// @author Robert Freepartner, JaDa Systems/Raytheon/JPL
///
/// @copyright (c) 2018 JPL, All rights reserved

/// @brief Type for points, equivalent to uint8.
typedef uint8_t mpoint_t;

/// @brief Offset of a Coordinate Point in a structuring element, i.e., kernel.
typedef struct
{
    int dy; ///< row, or "vertical" offset
    int dx; ///< col, or "horizontal" offset
} OffsetCoordinate;

// @brief A structuring element represented as a set of offset points.
typedef struct
{
    int noffs; ///< Number of offset points 
    OffsetCoordinate* offsets; ///< Offset points set
} Kernel;

/// @brief Sets initial NULL values.
void construct_Kernel(Kernel *kernel);

/// @brief Frees offsets and frees kernel.
void destruct_Kernel(Kernel *kernel);

/// @brief create a disk-shaped kernel.
///
/// @param  radius Radius of the disk in pixels
/// @return kernel structuring element
Kernel *disk_kernel(int radius);

/// @brief create a square-shaped kernel.
///
/// Example: A radius of 1 will include one pixel before and one pixel after,
/// both in rows and columns, resulting in a 3x3 kernel.
///
/// @param  radius Radius of the disk in pixels
/// @return kernel structuring element
Kernel *square_kernel(int radius);

/// @brief Perform morphological dilation on an image
///
/// @param points two-dimensional image stored in a one-dimensional array
///               arranged row after row.
/// @param nrows  number of rows in the image
/// @param ncols  number of columns in the image
/// @param kernel structuring element
/// @return new set of points after dilation
mpoint_t *dilation(const mpoint_t *points, int nrows, int ncols, const Kernel* kernel);

/// @brief Perform morphological erosion on an image
///
/// @param points two-dimensional image stored in a one-dimensional array
///               arranged row after row.
/// @param nrows  number of rows in the image
/// @param ncols  number of columns in the image
/// @param kernel structuring elemnt
/// @return new set of points after erosion
mpoint_t *erosion(const mpoint_t *points, int nrows, int ncols, const Kernel* kernel);

/// @brief Perform a morphological close operation on an image by performing dilation and erosion.
///
/// This is the same as imclose in MatLab.
///
/// @param points two-dimensional image stored in a one-dimensional array
///               arranged row after row.
/// @param nrows  number of rows in the image
/// @param ncols  number of columns in the image
/// @param kernel structuring elemnt
/// @return new set of points after dilation
mpoint_t *imclose(const mpoint_t *points, int nrows, int ncols, const Kernel* kernel);

