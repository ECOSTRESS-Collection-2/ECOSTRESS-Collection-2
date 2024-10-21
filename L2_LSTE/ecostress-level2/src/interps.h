#pragma once
#include "matrix.h"

/**
 * @file This header provides the declaration for the interp functions.
 * @author Robert Freepartner
 */

/**
 * @brief Performs 2D interpolation equivalent to interp2 in MatLab
 * @author Robert Freepartner
 * @date June 2016
 *
 * 2D datasets are organized as a vector of doubles sized
 * [nsampsx * nsampsy] organized so that there
 * are nsampsx rows of nsampsy columns.
 *
 * @param[in] X         Sample data X coordinates (rows)
 * @param[in] Y	        Sample data Y coordinates (columns per row)
 * @param[in] V	        Sample data
 * @param[in] nsampsx	Number of sample rows
 * @param[in] nsampsy   Number of sample columns
 * @param[in] Xq        Query point X coordinates
 * @param[in] Yq        Query point Y coordinates
 * @param[out] Vq       Interpolated values for each query point
 * @param[in] nqpoints  Number of query points
 */
void interp2(
    const double *X, const double *Y, double *V,
    int nsampsx, int nsampsy,
    const double *Xq, const double *Yq, double *Vq,
    int nqpoints);

/**
 * @brief Performs 2D interpolation on multiple query sets using
 * a single 2D sample set.
 * @author Robert Freepartner
 * @date November 2015
 *
 * 2D datasets are organized as a vector of doubles sized
 * [nsampsx * nsampsy] organized so that there
 * are nsampsx rows of nsampsy columns.
 *
 * @param[in] X         Sample data X coordinates (rows per layer)
 * @param[in] Y	        Sample data Y coordinates (columns per row)
 * @param[in] nsampsx	Number of sample rows
 * @param[in] nsampsy   Number of sample columns
 * @param[in] Xq        Query point X coordinates
 * @param[in] Yq        Query point Y coordinates
 * @param[in] nqpoints  Number of query points
 * @param[in] insets	List of sample value (V) datasets (nsampsx x nsampsy)
 * @param[in] outsets   List of output (Vq) datasets (nqpoints)
 * @param[in] nsets     Number of datasets
 */

void multi_interp2(
    const double *X, const double *Y, int nsampsx, int nsampsy,
    const double *Xq, const double *Yq, int nqpoints,
	double *insets[], double *outsets[], int nsets);

/**
 * @brief Performs 3D interpolation equivalent to interp3 in MatLab
 * @author Robert Freepartner
 * @date November 2015
 *
 * 3D datasets are organized as a vector of doubles sized
 * [nsampsx * nsampsy * nsampsz] organized so that there
 * are nsampsx rows of nsampsy columns nsampsz times.
 *
 * @param[in] X         Sample data X coordinates (rows per layer)
 * @param[in] Y	        Sample data Y coordinates (columns per row)
 * @param[in] Z         Sample data Z coordinates (layers)
 * @param[in] V	        Sample data
 * @param[in] nsampsx	Number of sample rows
 * @param[in] nsampsy   Number of sample columns
 * @param[in] nsampsz   Number of sample layers
 * @param[in] Xq        Query point X coordinates
 * @param[in] Yq        Query point Y coordinates
 * @param[in] Zq        Query point Z coordinates
 * @param[out] Vq       Interpolated values for each query point
 * @param[in] nqpoints  Number of query points
 */

void interp3(
    const double *X, const double *Y, const double *Z, const double *V,
    int nsampsx, int nsampsy, int nsampsz, 
    const double *Xq, const double *Yq, const double *Zq, double *Vq,
    int nqpoints);

/**
 * @brief Performs 3D interpolation on multiple query sets
 * @author Robert Freepartner
 * @date November 2015
 *
 * 3D datasets are organized as a vector of doubles sized
 * [nsampsx * nsampsy * nsampsz] organized so that there
 * are nsampsx rows of nsampsy columns nsampsz times.
 *
 * @param[in] X         Sample data X coordinates (rows per layer)
 * @param[in] Y	        Sample data Y coordinates (columns per row)
 * @param[in] Z         Sample data Z coordinates (layers)
 * @param[in] nsampsx	Number of sample rows
 * @param[in] nsampsy   Number of sample columns
 * @param[in] nsampsz   Number of sample layers
 * @param[in] Xq        Query point X coordinates
 * @param[in] Yq        Query point Y coordinates
 * @param[in] Zq        Query point Z coordinates
 * @param[in] nqpoints  Number of query points
 * @param[in] insets	List of sample value (V) datasets (nsampsx x nsampsy)
 * @param[in] outsets   List of output (Vq) datasets (nqpoints)
 * @param[in] nsets     Number of datasets
 */

void multi_interp3(
    const double *X, const double *Y, const double *Z, 
    int nsampsx, int nsampsy, int nsampsz, 
    const double *Xq, const double *Yq, const double *Zq, int nqpoints,
    const double *insets[], double *outsets[], int nsets);

