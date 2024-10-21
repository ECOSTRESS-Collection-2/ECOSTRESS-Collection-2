// maptogrid.h
#pragma once
#include "lste_lib.h"
#include "matrix.h"

/** 
 * @file The C language version of griddata maps 2D data stored in lat/lon form
 * onto the swath grid. The algorithm is based on map_ancillary_data developed 
 * by:
 *  @author Jim Ray
 *  Science Systems and Applications, Inc.
 *  MODIS Land Science Team
 *  NASA Goddard Space Flight Center
 *  Greenbelt, Maryland 20771
 *  jim@ltdri.org
 * 
 * The modifications from map_ancillary_data to maptogrid for MOD21 PGE116 were
 * developed by:
 *  Robert Freepartner
 *  JPL
 *  October, 2015
 *
 */

 /**
  * @struct GridMapData
  * @brief This struct contains all the parameters needed for maptogrid. 
  *
  * The dsetout element is an input parameter that provides pointers to 
  * the locations for interpolated values from dsetin.
  * @author  Robert Freepartner, JPL
  */
 typedef struct
 {
     double    *lat;       ///< Grid latitude query points, nlines x ncols
     double    *lon;       ///< Grid longitude query points, nlines x ncols
     int        nlines;    ///< Number of lines in the query grid, e.g., 406
     int        ncols;     ///< Number of columns in the query grid, e.g., 270
     int        inrows;    ///< Number of input dataset lines, e.g., 361.
     int        incols;    ///< Number of lnput dataset cols, e.g., 540.
     int        nsets;     ///< Number of datasets to map
     double   **dsetin;    ///< List of pointers to input datasets, inrows x incols
     double   **dsetout;   ///< List of pointers to output datasets, nlines x ncols
 } GridMapData;

/**
 * @brief meshgrid generates a 2D matrix of Lat/Lon values from two vectors
 * containing lats and lons respectively.
 *
 * @param latgrid[out] 2D matrix to receive the grid values
 * @param longrid[out] 2D matrix to receive the grid values
 * @param lats[in]     1D vector of lat values
 * @param lons[in]     1D vector of lon values
 */
void meshgrid(Matrix *latgrid, Matrix *longrid, Vector *lats, Vector *lons);

/**
 * @brief meshgrat performs the same funtion as the MatLab function. Generates a
 * 2D matrix of Lat/Lon values that span minLat,minLon to maxLat,Maxlon with a
 * row resolution of latRes and a column resolution of lonRes.
 *
 * @param latgrid[out] 2D matrix to receive the grid values
 * @param longrid[out] 2D matrix to receive the grid values
 * @param minLat    minimum latitude in the range -90.0, 90.0
 * @param minLon    minimum longitude in the range -180.0, 180.0
 * @param maxLat    maximum latitude 
 * @param maxLon    maximum longitude
 * @param latRes    latitude resolution for grid in degrees
 * @param lonRes    longitude resolution in degrees
 */
void meshgrat(Matrix *latgrid, Matrix *longrid, double minLat, double minLon, 
        double maxLat, double maxLon, double latRes, double lonRes);

/**
 * @brief Maps worldwide data to a grid of lat/lon positions. Multiple 
 * datasets can be interpolated to the lat/lon grid. The list of datasets 
 * must have at least one entry.
 *
 * @author  Robert Freepartner, JPL
 * @param   gmd     GridMapData struct pointer.
 * @return          0 = success, -1 = error
 */
int mapworldtogrid(GridMapData *gmd);

/**
 * @brief Maps a subset of worldwide data to a grid of lat/lon positions. 
 * Multiple datasets can be interpolated to the lat/lon grid. The list of 
 * datasets must have at least one entry.
 *
 * @author  Robert Freepartner, JPL
 * @param   gmd      GridMapData struct pointer.
 * @param   startLat Starting latitude of the input dataset(s)
 * @param   startLon Starting longitude of the input dataset(s)
 * @param   latres   Latitude resolution (e.g., 0.5)
 * @param   lonres   Longitude resolution (e.g., 0.625)
 * @return           0 = success, -1 = error
 */
int maptogrid(GridMapData *gmd, 
    double startLat, double startLon, double latres, double lonres);

/**
 * @brief Takes a grid of lat/lon query positions and a grid of
 * lat/lon sample positions. Using a 'nearest' method, it finds
 * the smaple nearest each query. Multiple datasets that match 
 * the sample grid size use the nearest row/col to map each
 * input dataset to a corresponding, initially undefined, 
 * output dataset. The number of datasets (GridMapData.nsets)
 * must be >= 1.
 *
 * @author  Robert Freepartner, JPL
 *
 * @param   gmd     GridMapData struct pointer.
 * @param   latsamp Sample latitude points [inrows x incols]
 * @param   lonsamp Sample Longitude points [inrows x incols]
 * @return          0 = success, -1 = error
 */
int mapnearest(GridMapData *gmd, double *latsamp, double *lonsamp);

