#pragma once
/** 
 * @file ASTER GED Emissivity and NDVI data is loaded and stored for 5 
 * bands (10, 11, 12, 13, 14) and lat/lon by 1/100th degree.
 * @author Robert Freepartner, JPL, January 2016
 * @author Dr. Glynn Hulley, JPL
 * Algorithm design by Glynn Hulley, JPL
 *
 * Copyright (c) JPL, All rights reserved
 */
 
#include <stdlib.h>
#include "bool.h"
#include "matrix.h"
 
#define ASTER_NAN -32768
#define ASTER_NBANDS 5
 
#define ASTER_GED_EMIS_DEFAULT 0.99
#define ASTER_GED_EMIS_NOMINAL 0.95
#define ASTER_GED_EMIS_OCEAN 0.99
#define ASTER_GED_EMIS_MIN 0.5
#define ASTER_GED_EMIS_MAX 1.0
#define ASTER_GED_EMIS_SCALE (1.0/1000.0)
 
#define ASTER_GED_NDVI_DEFAULT 0.0
#define ASTER_GED_NDVI_NOMINAL 0.5
#define ASTER_GED_NDVI_OCEAN 0.0
#define ASTER_GED_NDVI_MIN 0.0
#define ASTER_GED_NDVI_SCALE (1.0/100.0)
 
 typedef enum
 {
     MODIS,
     VIIRS,
     ECOSTRESS,
     N_SENSORS
 } Sensor;
 
/**
 * @struct AsterGed 
 * @brief This struct contains the ASTER GED data for a geographic area.
 *
 * ASTER GED data has a granularity of 1/100th degree. The data will always
 * be loaded to the start of a degree. The dimensions of the data loaded
 * are available from startLat, startLon, endLat, endLon. The number of
 * values in a longitude row is nLons. The number of latitude rows is nLats.
 * 
 */
typedef struct
{
    Sensor  sensor;     ///< Sensor type indicates adjustment
    double  startLat;   ///< Southmost latitude of the loaded data
    double  startLon;   ///< Westmost longitude of the loaded data
    double  endLat;     ///< Northmost latitude limit
    double  endLon;     ///< Eastmost longitude limit
    size_t  nLats;      ///< Number of 1/100th deg latitude rows
    size_t  nLons;      ///< Number of 1/100th degree longitude cols
    Matrix  emis;       ///< 2D Matrix of Emissivity/Mean by lat, lon
    Matrix  ndvi;       ///< 2D Matrix of NDVI/Mean by lat, lon
    bool    read_ndvi;  ///< Indicates that the NDVI data should be 
                        ///  retrieved from ASTER GED
    int     nFiles;     ///< Number of ASTER_GED files loaded 
    int     nDefault;   ///< Number of missing files replaced with default
    int     verbose;    ///< verbose level. <= 0 = no log output
} AsterGed;

/**
 * @brief Initialization (constructor) of the AsterGrid struct.
 *
 * Call this before making any other calls to ag_ functions.
 * @param ag    Pointer to AsterGrid struct.
 */
void ag_init(AsterGed *ag);

/**
 * @brief Clear (destructor) releasing allocated memory.
 * @param ag    Pointer to AsterGrid struct.
 */
void ag_clear(AsterGed *ag);

/* @brief Loads data from files in the asterdir directory to cover
 * the given lat/lon range.
 *
 * @param ag        Pointer to AsterGrid struct.
 * @param asterdir  Directory path containing ASTER GED data files.
 * @param minLat    Southmost latitude to be included
 * @param minLon    Westmost longitude
 * @param maxLat    Northmost latitude + 1
 * @param maxLon    Eastmost longitude + 1
 * @return          0=successs, -1=error
 */
int ag_load(AsterGed *ag, const char *asterdir, double minLat, double minLon,
            double maxLat, double maxLon);
/**
 * @brief Get the emissivity value for lat/lon.
 *
 * @param ag            Pointer to AsterGrid struct.
 * @param lat           Latitude, -90.0 to 90.0
 * @param lon           Longitude, -180.0 to 180.0
 * @return              Emissivity value, ASTER_NAN if params not in loaded area.
 */
double ag_get_emis(AsterGed *ag, double lat, double lon);

/**
 * @brief Get the NDVI value for lat/lon.
 *
 * @param ag            Pointer to AsterGrid struct.
 * @param lat           Latitude, -90.0 to 90.0
 * @param lon           Longitude, -180.0 to 180.0
 * @return              NDVI value, ASTER_NAN if params not in loaded area.
 */
double ag_get_ndvi(AsterGed *ag, double lat, double lon);
/**
 * @brief Reads ASTER GED data and maps it to the swath grid from GEO inputs.
 *
 * @param[in]   ASTER_dir       Directory that contains ASTER data
 * @param[in]   gran_lat        Latitude for the granule
 * @param[in]   gran_lon        Longitude for the granule
 * @param[out]  emis_aster      Aster EMIS data interpolated to the swath
 * @param[in]   oceanpix        Ocean map: 1=ocean, 0=land
 * @param[in]   snow_obs        Snow map: 1=snow, 0=not snow
 * @param[in]   ndvi_obs        NDVI for swath (NULL = do not produce NDVI)
 * @param[in]   adjust_aster    When true, performs NDVI adjustment to emis_aster
 * @param[in]   sensor_type     Type of sensor to compute for
 */
int read_aster_ged(const char *ASTER_dir, Matrix *gran_lat, Matrix *gran_lon,
    Matrix *emis_aster, MatUint8 *oceanpix, MatUint8 *snow_obs, 
    Matrix *ndvi_obs, bool adjust_aster, Sensor sensor_type);
