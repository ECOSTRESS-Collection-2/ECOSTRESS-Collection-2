#pragma once
/// @file Utility functions for the input of MERRA2 and other datasets
/// used by the TES (Temperature Emissivity Separation) algorithm.
//
// References: 
/// @author MATLAB code is developed by Dr. Tanvir Islam, JPL, [TI]
/// @author Initial C program is written by Simon Latyshev, Raytheon, [SL]
/// @author Robert Freepartner, JaDa Systems/Raytheon/JPL
//
// Copyright (c) 2016 JPL, All rights reserved

#include "matrix.h"

// --------------------------------------------------------------
// Define statements
// --------------------------------------------------------------

#define TES_STR_SIZE 1024	// string buffer size
#define TES_INT_NAN -32768	// CHECK: UPDATE: Consistent with HDF and/or MATLAB NaN
#define TES_REAL_NAN (0.0/0.0)
#define MAX_CHANNELS 5

/**
 * @brief Reduces cloud based on counting neighboring cloud pixels.
 *
 * @param[in,out]   cloud       Cloud mask where 0=not cloud
 * @param           n_neighbors Cloud must have > n_neighbors to remain cloud
 */
void constrict_cloud(MatUint8 *cloud, int n_neighbors);

/**
 * @brief Creates an extended cloud mask.
 *
 * @param[in]   cloud_in    Input cloud mask where 0=not cloud
 * @param[out]  cextend     Matrix for extended cloud
 * @param       cloud_extend    Radius to extend in pixels.
 * @param       extended_value  Value to put into extended pixels.
 */
void extend_cloud(MatUint8 *cloud_in, MatUint8 *cextend, int cloud_extend, 
    uint8 extended_value);

// --------------------------------------------------------------
// TesDoyDate structure for date/time containing day-of-year
// --------------------------------------------------------------
/// @brief contains a date/time to the minute
typedef struct 
{
    int year;   ///< range [0, 3000]
    int doy;    ///< range [1, 366]
    int hour;   ///< range [0, 23]
    int minute; ///< range [0, 59]
} TesDoyDate;

/// @brief Copy date, date = date_src
void tes_copy_doy_date(TesDoyDate *date, const TesDoyDate *date_src);

/// @brief Set date, initializing values
/// @return -1=invalid date, 0=valid
int tes_set_doy_date(TesDoyDate *date, int year, int doy, int hour, int minute);

// --------------------------------------------------------------
// TesDate structure for date/time containing month and day-of-month
// --------------------------------------------------------------
/// @brief contains a date/time to the minute
typedef struct 
{
    int year;   ///< range [0, 3000]
    int month;  ///< range [1, 12]
    int day;    ///< range [1, 31]
    int hour;   ///< range [0, 23]
    int minute; ///< range [0, 59]
} TesDate;

/// @brief Copy date, date = date_src
void tes_copy_date(TesDate *date, const TesDate *date_src);

/// @brief Set date, initializing values
/// @return -1=invalid date, 0=valid
int tes_set_date(TesDate *date, int _year, int _month, int _day, 
        int _hour, int _minute);

/**
 * @brief add minutes to a timestamp
 */

void tes_add_minutes(TesDate *date, int mins);

/* TBD
void tes_convert_date_to_doy(TesDoyDate *doydate, TesDate *mddate); 

void tes_convert_doy_to_mmdd(TesDate *mddate, TesDoyDate *doydate); 
*/

/**
 * @brief contains data for radiance, geolocation,
 * etc. from the GEO and RAD input files. Note:
 * ECOSTRESS has 6 radiance bands.
 */
typedef struct
{
    // From L1_GEO
    Matrix Lat;         ///< Latitude
    Matrix Lon;         ///< Longitude
    Matrix El;          ///< Height
    Matrix Satzen;      ///< SatelliteZenithAngle
    Matrix Solzen;      ///< Solar Zenith
    MatUint8 Water;     ///< land_fraction < 50% (1=water)
    char OrbitCorrectionPerformed[TES_STR_SIZE];
    // From L1_RAD
    Matrix SWIR;        	///< SWIR/swir_dn        1.6 μm SWIR Radiance
    Matrix Rad[MAX_CHANNELS];  ///< Radiance
    // Data Quality: 0=good 1=stripe 2=missing/bad
    MatInt32 DataQ[MAX_CHANNELS];    ///< Data_Quality_1 thru Data_Quality_5
    int nchannels; // 3 or 5
    char FieldOfViewObstruction[TES_STR_SIZE];
} RAD;

/**
 * @brief Initialize the RAD struct.
 *
 * @param rad   RAD struct pointer
 */
void init_rad(RAD *rad);

/**
 * @brief Finalize the ViirsRAD struct. Deallocates memory.
 *
 * @param rad   ViirsRAD struct pointer
 */
void clear_rad(RAD *rad);

/**
 * @brief Read and populate the RAD data from the L1B GEO file.
 *
 * Note: init_rad must have been called once for the RAD struct.
 *
 * @param rad           RAD struct pointer
 * @param GEO_filename  Pathname of the L1B GEO file.
 */
void tes_read_geo_data(RAD *rad, const char *GEO_filename);

/**
 * @brief Use BandSelection attribute to determine number of bands to use.
 *
 * @param RAD_filename  //< Pathname of the L1B RAD file.
 * @return number of bands for L2 processing
 */
int get_rad_nbands(const char *RAD_filename);

/**
 * @brief Read and populate the RAD data from the L1B RAD file.
 *
 * Note: init_rad must have been called once for the RAD struct.
 *
 * @param rad           //< RAD struct pointer
 * @param RAD_filename  //< Pathname of the L1B RAD file.
 * @param n_bands       //< Number of radiance bands to read
 * @param band_numbers  //< List of band indexes (0, 1, .. etc.)
 */
void tes_read_rad_data(RAD *rad, const char *RAD_filename,
        int n_bands, int *band_numbers);

// --------------------------------------------------------------
// TesATM: ATM atmospheric data structure 
// --------------------------------------------------------------

// MATLAB: ATM1 defined in read_atmos_merra

typedef struct 
{
    Vector lat;     // lat
    Vector lon;     // lon
    Vector lev;     // lev
    Mat3d t;        // t
    Mat3d q;        // q
    Matrix sp;      // sp
    Matrix skt;     // skt
    Matrix t2;      // t2
    Matrix q2;      // q2
    Matrix tcw;     // tcw
} TesATM;

// Clear, deallocate
void tes_atm_clear(TesATM *atm);

// Check, compare, matrix sizes
// void tes_atm_check_sizes(const TesATM *atm);

// Init all arrays and prepare for receiving values.
void tes_atm_init(TesATM *atm);

/**
 * @brief CropSize contains min/max lat/lon for cropping.
 */
typedef struct 
{
    double minLat;
    double maxLat;
    double minLon;
    double maxLon;
} CropSize;

/**
 * @brief Crop ATM to CropSize limits.
 *
 * @param[out]  outAtm  Cropped result
 * @param[in]   inATM   Uncropped source
 * @param[in]   crop    Lat/lon crop limits
 */
void crop_atm(TesATM *inATM, TesATM *outATM, CropSize *crop);

/**
 * @brief Read Lat, Lon, Height, SatZen and Senszen from 
 * the GEO input file.
 *
 * @param   rad         RAD data structure
 * @param   GEO_file    Pathname to L1B_GEO file
 */
void tes_read_geo_data(RAD *rad, const char *GEO_file);

/**
 * @brief Read ATM data from GEOS data
 *
 * @param[out]  atm      ATM data structure
 * @param       tdate    Date/time of run
 * @param       geos_dir Top level directory for GEOS data
 * @param       crop     Lat/Lon limits for cropped ATM
 */
void tes_read_interp_atmos_geos(TesATM *atm, const TesDate *tdate, 
        const char *geos_dir, CropSize *crop);

/**
 * @brief Read ATM data from MERRA data
 *
 * @param   atm     ATM data structure
 * @param   tdate   Date/time of run
 * @param   merra_dir   Top level directory for MERRA data
 */
void tes_read_interp_atmos_merra(TesATM *atm, const TesDate *tdate, 
        const char *merra_dir);

/**
 * @brief Read ATM data from NCEP data in GRIB1 format
 *
 * @param[out]  atm     ATM data structure
 * @param[in]   tdate   Granule timestamp
 * @param[in]   nwp_dir Pathname to the yop level NCEP directory
 */
void tes_read_interp_atmos_ncep(TesATM *atm, 
        const TesDate *tdate, const char *ncep_dir);

/**
 * @brief Get name of MERRA2 file name
 * 
 * @return Name of MERRA2 file.
 */
char *get_nwp_name();
