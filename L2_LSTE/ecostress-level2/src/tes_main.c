// tes_main.c
// 
/// @file This is the main source code for Land Surface Temperature and 
/// Emissivity processing.
//  % TES-WVS Algoritm (NPP/ECOSTRESS)
//  % MATLAB Version For MODAPS/SIPS
/// @author Dr. Glynn Hulley, 2016, NASA/JPL
/// @author Dr. Tanvir Islam, 2016, NASA/JPL
/// @author Dr. Nabin Malakar, 2016, NASA/JPL
//
//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  % Temperature/Emissivity Separation (TES) Algorithm:
//  %  - NWP input atmospheric profiles
//  %  - Auxiliary data
//  %  - RTTOV code
//  %  - Water Vapor Scaling (WVS) using 'EM/EMC/WVD (Islam)' coefficients
//  %
//  % Generates Land Surface Temperature and Emissivity
//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
/// @author Simon Latyshev, Raytheon -- Original RTTOV sub functions
/// @author Robert Freepartner, JaDa Systems/Raytheon/JPL -- VIIRS and EcoStress
/// implementations
//
/// @copyright (c) 2021, JPL, All rights reserved
//
//  2016-08-04 R. Freepartner Modified to perform TG_WVS
//  2016-10-10 R. Freepartner Modified to support NCEP data.
//  2017-07-28 R. Freepartner Modified to support GEOS data.
//  2018-01-25 R. Freepartner Add ASTER-GED input. Process data_quality flags.
//  2018-02-28 R. Freepartner Remaining updates for final L2_PGE pre-launch updates.
//  2018-04-24 R. Freepartner Corrections to metadata
//  2018-06-26 R. Freepartner Corrections to dataset metadata values
//  2018-07-09 R. Freepartner Implement option 1 to deal with bands that have missing scans
//  2018-07-24 R. Freepartner Corrected LocalGranuleID for L2_CLOUD
//  2018-07-27 R. Freepartner Fix GEOS5 interpolation and mapnearest
//  2018-07-30 R. Freepartner Implement Cloud algorithm revisions
//  2018-08-24 R. Freepartner Corrections to planck function interpolations
//  2018-08-29 R. Freepartner Added OrbitCorrectionPerformed to product metadata
//  2018-09-07 R. Freepartner Reduced Tmax by 30% in Cloud processing
//  2018-09-10 R. Freepartner Revised error estimate computations & metadata fixes
//  2019-02-26 R. Freepartner Cover L1B data_quality in all bands. T>380=nominal QC.
//  2019-04-10 R. Freepartner Add option to run with 3 bands based on orbit number.
//  2019-05-03 R. Freepartner Add BandSpecification metadata
//  2019-06-10 R. Freepartner L2 metadata changes for DAAC
//  2020-12-10 R. Freepartner Enable HDF5 file compression
//  2021-04-15 R. Freepartner BUILD 7 -- Implemented Cloud_v2
//  2021-07-13 R. Freepartner BUILD 7 -- New Cloud LUT file

// Make sure this version is updated in osp/PgeRunParameters.xml for both
// L2_LSTE and L2_CLOUD.
#define PGE_VERSION "2.0.1"

#include "asterged.h"
#include "cloud.h"
#include "config.h"
#include "error.h"
#include "fileio.h"
#include "interps.h"
#include "lste_lib.h"
#include "maptogrid.h"
#include "metadata.h"
#include "morphology.h"
#include "rttov_util.h"
#include "smooth2d.h"
#include "tes_util.h"
#include "tg_wvs.h"

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>

#define COMMENT
//#define DEBUG
#ifdef DEBUG
extern int debug_mapnearest;
#endif
//#define DEBUG_TES
//#define DEBUG_RTTOV

//#define DEBUG_CLOUD

//#define DEBUG_LINE 2600
//#define DEBUG_PIXEL 3000

//#define TESTING
#ifdef TESTING
#include "describe.h"
#endif

#define ENABLE_TIMING
#ifdef ENABLE_TIMING
#include "timer.h"
#endif

//#define PRODUCE_GEOLOCATION

static const double real_nan = REAL_NAN;
static const double nan1 = -32768;

// Number of channels, or bands, to process
#define MAX_CHANNELS 5
// Index for each band in Y
#define BAND_1 0
#define BAND_2 1
#define BAND_3 2
#define BAND_4 3
#define BAND_5 4

// With v1.2.0, the number of channels used is variable.
int n_channels = 5;
int band[MAX_CHANNELS] = {BAND_1, BAND_2, BAND_3, BAND_4, BAND_5};
// this is the reverse access, to get the correct index for a band number.
int band_index[MAX_CHANNELS] = {BAND_1, BAND_2, BAND_3, BAND_4, BAND_5};

// Select clearest band index for LST.
// Clearest band index for ECOSTRESS is band 4
int lst_band_index = 3;

// ========================================================
// RTTOV data structures
// ========================================================
// ATM structure, arrays for read-write RTTOV binary profile [SL]
typedef struct
{
    // Array sizes
    int imax;               // size of 1st dimension of 3d arrays
    int kmax;               // size of 2nd dimension of 3d arrays
    int jmax;               // size of 3rd dimension of 3d arrays
    int nmax;               // number of rows for 2d matrix P[nmax][imax]

    // Array size parameters
    int nprof;              // nprof = dm1*dm2
    int nlev;               // REPLACE [SL] nvel as 2d array

    // Arrays 2d
    double *T2;             // t2
    double *Q2;             // q2
    double *P;              // Pressure

    double *PSurf_sp;       // REPLACE [SL] sp
    double *HSurf_el;       // REPLACE [SL] el
    double *TSurf_skt;      // REPLACE [SL] skt

    double *ST;             // REPLACE [SL] ST, lsm
    double *Lat_atmos;      // Lat_atmos, lat
    double *Lon_atmos;      // Lon_atmos, lon
    double *SatZen;         // SatZen, Senszen
    double *Bemis;

    // Grid maps need original Lat/Lon sample
    double *Lat_samp;
    double *Lon_samp;

    // Arrays 3d
    double **T;        // t
    double **Q;        // q

    // Arrays 1d representing 3d Q, T arrays
    double *Tt;                            // ADD [SL] transposed, converted T 3d
    double *Qt;                            // ADD [SL] transposed, converted T 3d

    // TCW (Total Column Water) is used to produce PWV
    double *TCW;
} ATM_rttov;

// RAD structure, 2d arrays copied from 
// 1. Lat/Lon main arrays [SL]
//    Note: using main Lat, Lon arrays, not scaled by Write Atmos function
// 2. ATM->SatZen == main Senszen array
// MATLAB: structure sRAD from runRTTOVunixbin function

typedef struct
{
    // Array sizes
    int kmax;                // dm1, size of 1st dimension of 2d arrays
    int jmax;                // dm2, size of 2nd dimension of 2d arrays
    int nmax;                // ADD [SL] nmax = nprof = dm1*dm2 = kmax*jmax

    // Arrays 2d, copied from ATM structure
    double *Lat;             // sRAD.Lat = sATM1.lat
    double *Lon;             // sRAD.Lon = sATM1.lon
    double *SatZen;          // sRAD.SatZen = sATM1.SatZen

} RAD_rttov;

// RadOut, RTM 1d structure, 1d arrays with RTTOV data output [SL]
// MATLAB: variables of read_intep_rttov function
typedef struct
{
    int nprof;      // number of RTTOV profiles = nmax = kmax*jmax
    // 2D Arrays sized [n_channels][nprof]:
    double **TB_rttov;
    double **RadTot_rttov;
    double **RadUp_rttov;
    double **RadRefDn_rttov;
    double **Trans_rttov;
} RadOut_rttov;

// RTM structure, 2d arrays for read-interpret RTTOV results [SL]
// MATLAB: structure sRTM from read_intep_rttov function
typedef struct
{
    // Array sizes
    int kmax;                // dm1, size of 1st dimension of 2d arrays
    int jmax;                // dm2, size of 2nd dimension of 2d arrays
    int nmax;                // ADD [SL] nmax = nprof = dm1*dm2 = kmax*jmax

    // Arrays 2d, initialized with RadOut arrays
    //double **TB_rt;         // reshaped TB -- [BF] not in use
    //double **RadTot_rt;     // reshaped RadTot  -- [BF] not in use
    double **RadUp_rt;      // reshaped RadUp
    double **RadRefDn_rt;   // reshaped RadRefDn
    double **Trans_rt;      // reshaped Trans
    double *PWV;            // reshaped TCW

} RTM_rttov;

// ========================================================
// Function Prototypes
// ========================================================
void apply_tes_algorithm(
        // outputs
        Matrix *Ts, Mat3d *emisf, MatUint16 *QC,
        // inputs
        Mat3d *surfradi, Mat3d *skyr, MatUint8 *gp_water,
        Mat3d *t1r, Mat3d *Y);
void clear_ATM_rttov(ATM_rttov *ATM);
void init_ATM_rttov(ATM_rttov *ATM, int imax, int kmax, int jmax);
void set_rttov_atmos(
        ATM_rttov *ATM, ///< [out] result struct
        Mat3d *WV_mix,  ///< Q water vapor mix rate (PPMV)
        Mat3d *T_d,     ///< Temperature (K)
        double Pres[],  ///< Pressure at each level for T, Q (hPa)
        Matrix *Senszen,///< Sensor zenith (deg)
        Matrix *Psurf,  ///< Surface pressure (hPa)
        Matrix *Hsurf,  ///< Height (m)
        Matrix *Lat,    ///< Latitude (deg)
        Matrix *Lon,    ///< Longitude (deg)
        Matrix *tcw_in, ///< Total column water
        Matrix *skt,    ///< Skin temperature (K) -- may be empty
        Matrix *t2,     ///< 2 meter temperature (K) -- may be empty
        Matrix *q2);    ///< 2 meter water vapor (PPMV) -- may be empty
void clear_RAD_rttov(RAD_rttov *RAD);
void set_RAD_rttov(RAD_rttov *RAD, double *Lat, double *Lon, 
        double *SatZen, int kmax, int jmax);
void init_RadOut_rttov(RadOut_rttov *RadOut, int nprof);
void clear_RadOut_rttov(RadOut_rttov *RadOut);
void clear_RTM_rttov(RTM_rttov *RTM);
void init_RTM_rttov(RTM_rttov *RTM, int kmax, int jmax);
void write_rttov_atmos(const char *fpath_rttov_profile, 
        ATM_rttov *ATM, int wvs_case);
void write_rttov_profile(const char *fpath, ATM_rttov *ATM);
void run_rttov_script(const char *script, const char *exe_file, const char *coef_file);
void read_interp_rttov(const char* fpath_rttov_data, RTM_rttov *RTM, 
        RAD_rttov *RAD, ATM_rttov *ATM);
int NEM_planck(double, double [], double [], int, 
    double *, double *, double *, int *);
double getRuntimeParameter_f64(const char *name, double default_value);
int getRuntimeParameter_int(const char *name, int default_value);
char *getRuntimeParameter(const char *name);
// ========================================================
// Global data
// ========================================================
// Pathnames for this run
char log_filename[PATH_MAX] = "l2_pge.log";
char ASTER_dir[PATH_MAX] = "/project/test/ASTER_GED/AST_L3_1km";
char NWP_dir[PATH_MAX] = "";
char OSP_dir[PATH_MAX] = ".";
char RTTOV_dir[PATH_MAX] = ".";

// Input files
char GEO_filename[PATH_MAX] = "";
char RAD_filename[PATH_MAX] = "";
char input_pointer[2*PATH_MAX];

// Program output
// Filename format is: "ECOSTRESS_L2_TES_orbit_sceneTtime_build_version.h5"
char output_root[] = "ECOSTRESS_L2_LSTE_";
char output_filename[PATH_MAX] = 
    "ECOSTRESS_L2_LSTE_OOOOO_SSS_YYYYMMDDThhmmss_BBbb_VV.h5";
char cloud_root[] = "ECOSTRESS_L2_CLOUD_";
char cloud_filename[PATH_MAX] = 
    "ECOSTRESS_L2_CLOUD_OOOOO_SSS_YYYYMMDDThhmmss_BBbb_VV.h5";
char ecmwf_cutout_root[] = "L2_ECMWF_CUT";
char geos_cutout_root[] = "L2_GEOS_CUT";
char merra_cutout_root[] = "L2_MERRA_CUT";
char ncep_cutout_root[] = "L2_NCEP_CUT";
char cutout_filename[PATH_MAX] = 
    "L2_MERRA_CUT_OOOOO_SSS_YYYYMMDDThhmmss_BBbb_VV.h5";
char product_path[PATH_MAX];
char product_counter[CONFIG_VALUELEN];

// NWP keys indicate the type of NWP data being used.
// One of these must appear in the NWP_dir name.
char nwp_key_ecmwf[PATH_MAX] = "ECWMF";
char nwp_key_geos[PATH_MAX] = "GEOS";
char nwp_key_merra[PATH_MAX] = "MERRA";
char nwp_key_ncep[PATH_MAX] = "NCEP";

// Geometry Parameters
char orbit_number[CONFIG_VALUELEN] = "TBD";
char scene_id[CONFIG_VALUELEN];
char build_id[CONFIG_VALUELEN];

// Current line and pixel -- global for debug reference in subroutines
int iline, ipixel;

#ifdef TESTING
char debugfn[PATH_MAX];
#endif

// Operational Support Product (OSP) filenames
char L2_OSP_DIR[PATH_MAX] = ".";
char bt_lut_file[PATH_MAX] = "";
char cloud_btdiff_file[PATH_MAX] = "";
char cloud_lut_file[PATH_MAX] = "";
//char sky_tes_file[PATH_MAX] = "";
char rad_lut_file[PATH_MAX] = "";
char rttov_coef_filename[PATH_MAX] = "";
char rttov_radcon_filename[PATH_MAX] = "";
char wvs_coef_file[PATH_MAX] = "";

// RTTOV executable
char rttov_exe[PATH_MAX] = "";
char rttov_script[PATH_MAX] = "";

//% Specify physical coefficients
double emax_veg = 0.985;
double emax_bare = 0.97;
double co_veg[3] = {0.9950, 0.7264, 0.8002};
double co_bare[3] = {0.9950, 0.7264, 0.8002};
double c1 = 3.7418e-22; // Plank coeff 1
double c2 = 0.014388; // Plank coeff 2

// Other parameters
double elapsed_seconds = 300.0;
int it = 13;        // N iterations for Planck function
//double VI = 0.0005; // spectral variance threshold
//double NDV_thresh = 0.3;
//double NDVIvals[3] = {0.0, 0.4, 0.8};
//% Emissivities used for bare, transition, gray in QC error dataplanes for 
//  MODIS bands 29, 31, 32
//% einterp{1} = [0.8 0.85 0.97];
//% einterp{2} = [0.94 0.96 0.98];
//% einterp{3} = [0.96 0.97 0.98];
//
// TODO -- Get values for ECOSTRESS, or implement alternate method.
/*
double einterp[MAX_CHANNELS][3] = 
{
    {0.80, 0.85, 0.97},
    {0.94, 0.96, 0.98},
    {0.96, 0.97, 0.98},
    {0.96, 0.97, 0.98},
    {0.96, 0.97, 0.98}
};
*/

double PWVthresh1 = 2.5;
double PWVthresh2 = 2.0;
double PWVthresh = 2.0;
double TGthresh = 290.0;
double Terr_lims[3] = {1.0, 1.5, 2.5};
// % Use these for now, there are only 3, they are not band dependent
double emiserr_lims[3] = {0.013, 0.015, 0.017};  
//int smooth_radius = 10;
int smooth_scale1 = 750;
int smooth_scale2 = 150;
int smooth_scale = 150;
int cloud_extend = 5;
bool run_tgwvs = true;

//% Band model parameters
double g1 = 1.0;
double g2 = 0.7;
// NOTE: When the revised wavelength for band 4 is available,
//       these values will need to be revised. 
double bmp[5] = {1.2818, 1.5693, 1.6595, 1.8217, 1.8031};
int last_5_band_orbit = 3894;

// Values for bit mask operationns.
const uint16 bits[16] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80,
                   0x100, 0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x8000};

// Start date of processing
TesDate tes_date;

// Lookup Table (LUT)
int nlut_lines = 0;
#define LUT_COLS 6
double *lut[LUT_COLS];
const int COL_1 = 0;
const int COL_2 = 1;
const int COL_3 = 2;
const int COL_4 = 3;
const int COL_5 = 4;
const int COL_6 = 5;

// SKY TES parameters
// double c_sky[MAX_CHANNELS][3];

// PgeRunParameters.xml parsed inputs.
ConfigSet *run_params;

// Data structs
TesATM nwpATM;
RAD vRAD;

// Metadata
MetadataSet cloud_metadata;
MetadataSet lste_metadata;

// Constant coefficients

// Coefficients for emis_wb
double emis_wb_coeffs[] = {0.0715, 0.0657, 0.1970, 0.3384, 0.3703, -0.0443};

// Coefficients for emissivity error cals.
double xe[MAX_CHANNELS][3] = 
{
    {0.0153,    0.0155,   -0.0018},
    {0.0121,    0.0055,    0.0004},
    {0.0128,    0.0036,    0.0009},
    {0.0110,    0.0017,    0.0005},
    {0.0114,   -0.0038,    0.0025}
};
//double emiserr_coef0[MAX__CHANNELS] = {0.0185, 0.0113, 0.0107, 0.0098, 0.0041};
//double emiserr_coef1[MAX__CHANNELS] = {0.0049, 0.0045, 0.0042, 0.0023, 0.0058};

// LST coefficients
//double xtt[6] = {0.3590, 0.6488, -0.0012, 0.0006, -0.0619, 0.0001};
double xt[3] = {0.3842, 0.5307, 0.0055};
//double lsterr_coef0 = 0.5226;
//double lsterr_coef1 = 0.2630;
//double lsterr_coef2 = 0.0014;

// Compression parameters
int chunk_factor = 10;
int compression_deflate = 9;

// ========================================================
// tes_main processing
// ========================================================
int main (int argc, char *argv[])
{
    int stat;
    int i, j, k;
    int b; // for band index

#ifdef ENABLE_TIMING
    time_t startmain = getusecs();
#endif

    // ========================================================
    // Initialization
    // ========================================================

    constructMetadataSet(&cloud_metadata, 0);
    constructMetadataSet(&lste_metadata, 0);

    // Get run parameters and filenames

    // Check for valid number of command line parameters.
    if (argc < 2)
    {
        fprintf(stderr, "usage: %s <pcf_filename>\n", argv[0]);
        exit(1);
    }
    
    // Setup logging
    set_log_filename(log_filename);
    echo_log_to_stdout();

    // Perform popen commands here. popen uses fork() which 
    // duplicates the program space in memory. If these calls 
    // are made after a large number of memory allocations, the
    // fork() is likely to fail.
    //
    // Get the system time stamp.
    char timestamp[256];
    FILE *pipe = popen("date -u +'%Y-%m-%dT%H:%M:%S.%NZ'", "r");
    if (pipe == NULL)
    {
        if (errno == ENOMEM)
            ERREXIT(112, "This host machine is too low on memory to support "
                    "popen() and fork(). Run on another host, or enable "
                    "overcommit memory on this system.", NONE);
        else
            ERREXIT(112, "A call to popen() failed with errno=%d", errno);
    }
    char *strres = fgets(timestamp, 255, pipe);
    assert(strres != NULL);
    pclose(pipe);
    // Chop off nanoseconds to milliseconds
    char *findz = strrchr(timestamp, 'Z');
    if (findz != NULL)
        strcpy(findz-6, "Z");

    // Set the processing environment to "uname -a"
    char uname[256];
    pipe = popen("uname -a", "r");
    if (pipe == NULL)
    {
        if (errno == ENOMEM)
            ERREXIT(112, "This host machine is too low on memory to support "
                    "popen() and fork(). Run on another host, or enable "
                    "overcommit memory on this system.", NONE);
        else
            ERREXIT(112, "A call to popen() failed with errno=%d", errno);
    }
    strres = fgets(uname, 255, pipe);
    pclose(pipe);
    // Trim off trailing spaces and unprintable characters.
    char *trim = uname + strlen(uname) - 1;
    while (trim >= uname && *trim <= ' ')
        *trim-- = '\0';

    // ========================================================
    // Read Configuration Parameters
    // ========================================================
    const char *run_config_filename = argv[1];

    // Set the parameter file name
    setConfigVerbose(1); // Error messages only
    ConfigSet *run_config = parseConfig(run_config_filename);
    assert(run_config != NULL);
    if (run_config->status == -1)
        ERREXIT(2, "I/O error on file %s", run_config_filename);
    if (run_config->status == -2)
        ERREXIT(3, "XML parsing error in file %s", run_config_filename);
    ConfigElement *elem;

    // Get NWP_DIR
    elem = findGroupElement(run_config, 
        "StaticAncillaryFileGroup", "NWP_DIR");
    if (elem ==  NULL)
    {
        ERREXIT(4, "Unable to find NWP_DIR in group "
              "StaticAncillaryFileGroup in %s", run_config_filename);
    }
    strcpy(NWP_dir, getValue(elem, 0));

    // Get L2_OSP_DIR
    elem = findGroupElement(run_config, 
        "StaticAncillaryFileGroup", "L2_OSP_DIR");
    if (elem ==  NULL)
    {
        ERREXIT(4, "Unable to get L2_OSP_DIR from group "
              "StaticAncillaryFileGroup in %s", run_config_filename);
    }
    strcpy(OSP_dir, getValue(elem, 0));
    strcat(OSP_dir, "/");

    // Get the output file path
    elem = findGroupElement(run_config, "ProductPathGroup", "ProductPath");
    if (elem == NULL)
        ERREXIT(4, "Could not find ProductPath "
        "in group ProductPathGroup in %s", run_config_filename);
    char *paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(5, "Error reading value of ProductPath from %s", run_config_filename);
    strcpy(product_path, paramstr);

    // Get product counter, i.e., version number
    elem = findGroupElement(run_config, "ProductPathGroup", "ProductCounter");
    if (elem == NULL)
        ERREXIT(4, "Could not find ProductCounter "
        "in group ProductPathGroup in %s", run_config_filename);
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(5, "Error reading value of ProductCounter from %s", run_config_filename);
    strcpy(product_counter, paramstr);
    int ivv = atoi(product_counter);

    // L1B_GEO
    elem = findGroupElement(run_config, "InputFileGroup", "L1B_GEO");
    if (elem == NULL)
        ERREXIT(4, "Could not find L1B_GEO "
        "in group InputFileGroup in %s", run_config_filename);
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(5, "Error reading value of L1B_GEO from %s", run_config_filename);
    strcpy(GEO_filename, paramstr);

    // L1B_RAD
    elem = findGroupElement(run_config, "InputFileGroup", "L1B_RAD");
    if (elem == NULL)
        ERREXIT(4, "Could not find L1B_RAD "
        "in group InputFileGroup in %s", run_config_filename);
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(5, "Error reading value of L1B_RAD from %s", run_config_filename);
    strcpy(RAD_filename, paramstr);

    // Geometry: OrbitNumber
    elem = findGroupElement(run_config, "Geometry", "OrbitNumber");
    if (elem == NULL)
        ERREXIT(4, "Could not find OrbitNumber "
        "in group Geometry in %s", run_config_filename);
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(5, "Error reading value of OrbitNumber from %s", run_config_filename);
    strcpy(orbit_number, paramstr);
    int iorbit = atoi(orbit_number);
    
    // Geometry: SceneId
    elem = findGroupElement(run_config, "Geometry", "SceneID");
    if (elem == NULL)
        elem = findGroupElement(run_config, "Geometry", "SceneId");
    if (elem == NULL)
        ERREXIT(4, "Could not find SceneID "
        "in group Geometry in %s", run_config_filename);
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(5, "Error reading value of SceneID from %s", run_config_filename);
    strcpy(scene_id, paramstr);
    int iscene = atoi(scene_id);

    elem = findGroupElement(run_config, "PrimaryExecutable", "BuildID");
    if (elem == NULL)
    {
        ERREXIT(5, "Error reading value of BuildID from %s", run_config_filename);
    }
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(5, "Error reading value of BuildID from %s", run_config_filename);
    strcpy(build_id, paramstr);
    int ibuild = atoi(build_id);
    
    // ========================================================
    // Read Runtime Parameters
    // ========================================================
    // Parse the PgeRunParameters.xml file.
    char run_params_filename[PATH_MAX];
    strcpy(run_params_filename, OSP_dir);
    strcat(run_params_filename, "PgeRunParameters.xml");
    run_params = parseConfig(run_params_filename);
    assert(run_params != NULL);
    if (run_params->status == -1)
        ERREXIT(6, "I/O error on file %s", run_params_filename);
    if (run_config->status == -2)
        ERREXIT(7, "XML parsing error in file %s", run_params_filename);
    
    // Make sure the run parameters version matches the executable.
    elem = findGroupElement(run_params, "L2_LSTE_Metadata", "PGEVersion");
    if (elem == NULL)
        ERREXIT(8, "Could not find L2_LSTE_Metadata/PGEVersion in %s", 
            run_params_filename);
    if (strcmp(getValue(elem, 0), PGE_VERSION) != 0)
        ERREXIT(8, "L2_PGE is version %s, but %s is version %s",
            PGE_VERSION, run_params_filename, getValue(elem, 0));

    // Get RTTOV directory. If not found, leave set to "."
    elem = findGroupElement(run_params, "RTTOV", "RttovDir");
    if (elem != NULL)
    {
        strcpy(RTTOV_dir, getValue(elem, 0));
    }

    // Get RTTOV Coefficients file.
    elem = findGroupElement(run_params, "RTTOV", "RttovCoef");
    if (elem == NULL)
        ERREXIT(8, "Could not find RttovCoef "
        "in group RTTOV in %s", run_params_filename);
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(9, "Error reading value of RttovCoef from "
        "group RTTOV in %s", run_params_filename);
    strcpy(rttov_coef_filename, paramstr);

    // Get RTTOV executable.
    elem = findGroupElement(run_params, "RTTOV", "RttovExe");
    if (elem == NULL)
        ERREXIT(8, "Could not find RttovExe "
        "in group RTTOV in %s", run_params_filename);
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(9, "Error reading value of RttovExe from "
        "group RTTOV in %s", run_params_filename);
    strcpy(rttov_exe, paramstr);

    // Get RTTOV Lookup Table.
    elem = findGroupElement(run_params, "RTTOV", "RttovRadconLUT");
    if (elem == NULL)
        ERREXIT(8, "Could not find RttovRadconLUT "
        "in group RTTOV in %s", run_params_filename);
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(9, "Error reading value of RttovRadconLUT from "
        "group RTTOV %s", run_params_filename);
    strcpy(rttov_radcon_filename, paramstr);

    // LookupFiles: BT_LUT
    elem = findGroupElement(run_params, "StaticLookupFiles", "BT_LUT");
    if (elem == NULL)
        ERREXIT(8, "Could not find BT_LUT "
        "in group StaticLookupFiles in %s", run_params_filename);
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(9, "Error reading value of BT_LUT "
        "from group StaticLookupFiles in %s", run_params_filename);
    strcpy(bt_lut_file, OSP_dir);
    strcat(bt_lut_file, paramstr);
    strcpy(rad_lut_file, OSP_dir); // Uses same file
    strcat(rad_lut_file, paramstr);

    // LookupFiles: cloud btdiff parameter file
    elem = findGroupElement(run_params, "StaticLookupFiles", "cloud_BTdiff");
    if (elem == NULL)
        ERREXIT(8, "Could not find cloud_BTdiff "
        "from group StaticLookupFiles in %s", run_params_filename);
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(9, "Error reading value of cloud_BTdiff "
        "from group StaticLookupFiles in %s", run_params_filename);
    strcpy(cloud_btdiff_file, OSP_dir);
    strcat(cloud_btdiff_file, paramstr);

    // Cloud_LUT
    elem = findGroupElement(run_params, "StaticLookupFiles", "cloud_LUT");
    if (elem == NULL)
        ERREXIT(8, "Could not find cloud_LUT "
        "from group StaticLookupFiles in %s", run_params_filename);
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(9, "Error reading value of cloud_LUT "
        "from group StaticLookupFiles in %s", run_params_filename);
    strcpy(cloud_lut_file, OSP_dir);
    strcat(cloud_lut_file, paramstr);

/*    
    // SKY TES parameter file
    elem = findGroupElement(run_params, "StaticLookupFiles", "SKY_tes");
    if (elem == NULL)
        ERREXIT(8, "Could not find SKY_tes "
        "from group StaticLookupFiles in %s", run_params_filename);
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(9, "Error reading value of SKY_tes "
        "from group StaticLookupFiles in %s", run_params_filename);
    strcpy(sky_tes_file, OSP_dir);
    strcat(sky_tes_file, paramstr);
*/
    // LookupFiles: WvsCoefFile
    elem = findGroupElement(run_params, "StaticLookupFiles", "WvsCoefFile");
    if (elem == NULL)
        ERREXIT(8, "Could not find WvsCoefFile "
        "in group StaticLookupFiles in %s", run_params_filename);
    paramstr = getValue(elem, 0);
    if (!paramstr)
        ERREXIT(9, "Error reading value of WvsCoefFile "
        "from group StaticLookupFiles in %s", run_params_filename);
    strcpy(wvs_coef_file, OSP_dir);
    strcat(wvs_coef_file, paramstr);

    // Get NWP key values
    elem = findGroupElement(run_params, "NWP", "ecmwf");
    if (elem !=  NULL) strcpy(nwp_key_ecmwf, getValue(elem, 0));
    elem = findGroupElement(run_params, "NWP", "geos");
    if (elem !=  NULL) strcpy(nwp_key_geos, getValue(elem, 0));
    elem = findGroupElement(run_params, "NWP", "merra");
    if (elem !=  NULL) strcpy(nwp_key_merra, getValue(elem, 0));
    elem = findGroupElement(run_params, "NWP", "ncep");
    if (elem !=  NULL) strcpy(nwp_key_ncep, getValue(elem, 0));

    // Get the RunParameters from the PgeRunParameters.xml data.
    /* cloud_extend not used by Cloud_v2
    cloud_extend = getRuntimeParameter_int("CloudExtend", cloud_extend);
    */
    elapsed_seconds = getRuntimeParameter_f64("ElapsedSeconds", elapsed_seconds);
    emax_veg = getRuntimeParameter_f64("EmaxVeg", emax_veg);
    emax_bare = getRuntimeParameter_f64("EmaxBare", emax_bare);
    g1 = getRuntimeParameter_f64("G1", g1);
    g2 = getRuntimeParameter_f64("G2", g2);
    last_5_band_orbit =  getRuntimeParameter_int("Last5BandOrbit", last_5_band_orbit);
    //lsterr_coef0 = getRuntimeParameter_f64("LstErrCoef0", lsterr_coef0);
    //lsterr_coef1 = getRuntimeParameter_f64("LstErrCoef1", lsterr_coef1);
    //lsterr_coef2 = getRuntimeParameter_f64("LstErrCoef2", lsterr_coef2);
    it = getRuntimeParameter_int("NPlanckIterations", it);
    //NDV_thresh = getRuntimeParameter_f64("NDV", NDV_thresh);
    PWVthresh1 = getRuntimeParameter_f64("PwvThresh1", PWVthresh1);
    PWVthresh2 = getRuntimeParameter_f64("PwvThresh2", PWVthresh2);
    TGthresh = getRuntimeParameter_f64("TgThresh", TGthresh);
    //smooth_radius = getRuntimeParameter_int("SmoothRadius", smooth_radius);
    smooth_scale1 = getRuntimeParameter_int("SmoothScale1", smooth_scale1);
    smooth_scale2 = getRuntimeParameter_int("SmoothScale2", smooth_scale2);
    //VI = getRuntimeParameter_f64("VI", VI);

    paramstr = getRuntimeParameter("AsterDir");
    if (paramstr != NULL)
    {
        strcpy(ASTER_dir, paramstr);
    }

/*
    paramstr = getRuntimeParameter("EmisErrCoef0");
    if (paramstr != NULL)
    {
        sscanf(paramstr, "%lf %lf %lf %lf %lf", &emiserr_coef0[0], 
                &emiserr_coef0[1], &emiserr_coef0[2], &emiserr_coef0[3], 
                &emiserr_coef0[4]);
    }

    paramstr = getRuntimeParameter("EmisErrCoef1");
    if (paramstr != NULL)
    {
        sscanf(paramstr, "%lf %lf %lf %lf %lf", &emiserr_coef1[0], 
                &emiserr_coef1[1], &emiserr_coef1[2], &emiserr_coef1[3], 
                &emiserr_coef1[4]);
    }
*/

    paramstr = getRuntimeParameter("XE1");
    if (paramstr != NULL)
    {
        sscanf(paramstr, "%lf %lf %lf", &xe[0][0], &xe[0][1], &xe[0][2]);
    }

    paramstr = getRuntimeParameter("XE2");
    if (paramstr != NULL)
    {
        sscanf(paramstr, "%lf %lf %lf", &xe[1][0], &xe[1][1], &xe[1][2]);
    }

    paramstr = getRuntimeParameter("XE3");
    if (paramstr != NULL)
    {
        sscanf(paramstr, "%lf %lf %lf", &xe[2][0], &xe[2][1], &xe[2][2]);
    }

    paramstr = getRuntimeParameter("XE4");
    if (paramstr != NULL)
    {
        sscanf(paramstr, "%lf %lf %lf", &xe[3][0], &xe[3][1], &xe[3][2]);
    }

    paramstr = getRuntimeParameter("XE5");
    if (paramstr != NULL)
    {
        sscanf(paramstr, "%lf %lf %lf", &xe[4][0], &xe[4][1], &xe[4][2]);
    }
/*
    paramstr = getRuntimeParameter("XTT");
    if (paramstr != NULL)
    {
        sscanf(paramstr, "%lf %lf %lf %lf %lf %lf", 
                &xtt[0], &xtt[1], &xtt[2], &xtt[3], &xtt[4], &xtt[5]);
    }
*/
    paramstr = getRuntimeParameter("XTT");
    if (paramstr != NULL)
    {
        sscanf(paramstr, "%lf %lf %lf", 
                &xt[0], &xt[1], &xt[2]);
    }

    paramstr = getRuntimeParameter("EmisErrLims");
    if (paramstr != NULL)
    {
        sscanf(paramstr, "%lf %lf %lf", 
            &emiserr_lims[0], &emiserr_lims[1], &emiserr_lims[2]);
    }

/*
    paramstr = getRuntimeParameter("NDVIVALS");
    if (paramstr != NULL)
    {
        sscanf(paramstr, "%lf %lf %lf", 
            &NDVIvals[0], &NDVIvals[1], &NDVIvals[2]);
    }
*/
    paramstr = getRuntimeParameter("RunTgWvs");
    if (paramstr != NULL)
    {
        if (*paramstr == 'Y' || *paramstr == 'y')
        {
            run_tgwvs = true;
        }
        else if (*paramstr == 'N' || *paramstr == 'n')
        {
            run_tgwvs = false;
        }
        else
            WARNING("RunTimeParameters value for RunTgWvs should be Y or N, not %s",
                paramstr);
    }

    paramstr = getRuntimeParameter("TerrLims");
    if (paramstr != NULL)
    {
        sscanf(paramstr, "%lf %lf %lf", 
            &Terr_lims[0], &Terr_lims[1], &Terr_lims[2]);
    }

    // Get Date/Time from the GEO filename.
    // Sample filename: ECOSTRESS_L1B_GEO_80001_001_20120811T203500_0001_01.h5

    int yy, mo, dd, hh, mm;
    char *p = strstr(GEO_filename, "_20"); // Find the date/time
    if (p == NULL)
        ERREXIT(10, "Cannot find date/time in GEO filename %s", GEO_filename);
    sscanf(p, "_%4d%2d%2dT%2d%2d", &yy, &mo, &dd, &hh, &mm);
    printf("ECOSTRESS Processing Year=%04d, Month=%02d, Day=%02d, Time=%02d:%02d\n",
           yy, mo, dd, hh, mm);
    fflush(stdout);
    if (tes_set_date(&tes_date, yy, mo, dd, hh, mm) < 0)
        ERREXIT(10, "Invalid date/time from GEO filename at %s", p);
    
    // Create the output names from the input name.
    char *prad = strstr(RAD_filename, "RAD_");
    if (prad == NULL)
        ERREXIT(11, "L1B_RAD name %s does not contain 'RAD_' in %s",
                RAD_filename, run_config_filename);
    // Want the date "T" time string.
    char *T = strchr(prad, 'T');
    if (T == NULL)
        ERREXIT(11, "L1B_RAD name %s does not contain 'T' in the timestamp.",
                RAD_filename);
    char *dateTtime = strdup(T-8);
    dateTtime[15] = '\0';
    sprintf(output_filename, "%s/%s%05d_%03d_%s_%04d_%02d.h5",
        product_path, output_root, iorbit, iscene, dateTtime, ibuild, ivv);
    sprintf(cloud_filename, "%s/%s%05d_%03d_%s_%04d_%02d.h5",
        product_path, cloud_root, iorbit, iscene, dateTtime, ibuild, ivv);

    // Output all parameters
    INFO("LSTE Parameters:", NONE);
    INFO("ASTER_DIR    = %s", ASTER_dir);
    INFO("NWP_DIR      = %s", NWP_dir);
    INFO("ORBIT_NUMBER = %s", orbit_number);
    INFO("LAST_5_BAND_ORBIT = %d", last_5_band_orbit);
    INFO("SCENE_ID = %s", scene_id);
    INFO("BUILD_ID = %s", build_id);
    INFO("PRODUCT_COUNTER = %s", product_counter);
    INFO("PGE_OUTPUT_FILE = %s", output_filename);
    INFO("L2_OSP_DIR = %s", OSP_dir);
    INFO("RTTOV_DIR = %s", RTTOV_dir);
    INFO("RTTOV_EXE = %s", rttov_exe);
    INFO("RTTOV_COEF_ECOSTRESS = %s", rttov_coef_filename);
    INFO("RTTOV_RADCON_LUT_ECOSTRESS = %s", rttov_radcon_filename);
    INFO("BT_LUT_FILE =%s", bt_lut_file);
    INFO("CLOUD_BTDIFF_FILE =%s", cloud_btdiff_file);
    INFO("CLOUD_LUT_FILE =%s", cloud_lut_file);
    //INFO("SKY_TES_FILE =%s", sky_tes_file);
    INFO("WVS_COEFF_FILE = %s", wvs_coef_file);
    //INFO("EMISERR_COEF0 %.4f,%.4f,%.4f,%.4f,%.4f",
    //    emiserr_coef0[0], emiserr_coef0[1], emiserr_coef0[2], 
    //    emiserr_coef0[3], emiserr_coef0[4]);
    //INFO("EMISERR_COEF1 %.4f,%.4f,%.4f,%.4f,%.4f",
    //    emiserr_coef1[0], emiserr_coef1[1], emiserr_coef1[2], 
    //    emiserr_coef1[3], emiserr_coef1[4]);
    INFO("EMISERR_LIMS  %.4f,%.4f,%.4f", 
        emiserr_lims[0], emiserr_lims[1], emiserr_lims[2]);
    //INFO("LSTERR_COEF0  %.4f", lsterr_coef0);
    //INFO("LSTERR_COEF1  %.4f", lsterr_coef1);
    //INFO("LSTERR_COEF2  %.4f", lsterr_coef2);
    //INFO("NDVIVALS      %.4f,%.4f,%.4f", NDVIvals[0], NDVIvals[1], NDVIvals[2]);
    //INFO("NDV_THRESH    %.4f", NDV_thresh);
    INFO("PLANCK_ITERS  %d", it);
    INFO("PWV_THRESH1   %.4f", PWVthresh1);
    INFO("PWV_THRESH2   %.4f", PWVthresh2);
    INFO("SMOOTH_SCALE1 %d", smooth_scale1);
    INFO("SMOOTH_SCALE2 %d", smooth_scale2);
    INFO("TG_THRESH     %.4f", TGthresh);
    INFO("TERR_LIMS     %.4f,%.4f,%.4f", Terr_lims[0], Terr_lims[1], Terr_lims[2]);
    INFO("XE1           %.4f,%.4f,%.4f", xe[0][0], xe[0][1], xe[0][2]);
    INFO("XE2           %.4f,%.4f,%.4f", xe[1][0], xe[1][1], xe[1][2]);
    INFO("XE3           %.4f,%.4f,%.4f", xe[2][0], xe[2][1], xe[2][2]);
    INFO("XE4           %.4f,%.4f,%.4f", xe[3][0], xe[3][1], xe[3][2]);
    INFO("XE5           %.4f,%.4f,%.4f", xe[4][0], xe[4][1], xe[4][2]);
    INFO("XT            %.4f,%.4f,%.4f", xt[0], xt[1], xt[2]);

    // ========================================================
    // Determine channel set to use for this run.
    // ========================================================
    n_channels = get_rad_nbands(RAD_filename);
    
    if (n_channels > 3 && iorbit > last_5_band_orbit)
    {
        WARNING("L1B_RAD indicates more than 3 bands for orbit > %d.",
            last_5_band_orbit);
    }

    if (n_channels != 3 && n_channels != 5)
    {
        WARNING("L1B_RAD indicates %d bands. L2_PGE requires 3 or 5 "
                "bands (excluding SWIR band).", n_channels);
        n_channels = 3;
    }

    if (n_channels == 3)
    {
        band[0] = BAND_2;
        band[1] = BAND_4;
        band[2] = BAND_5;
        // Set up the reverse index. -1 = band not loaded.
        band_index[BAND_1] = -1;
        band_index[BAND_2] = 0;
        band_index[BAND_3] = -1;
        band_index[BAND_4] = 1;
        band_index[BAND_5] = 2;
        lst_band_index = band_index[BAND_4];
        paramstr = getRuntimeParameter("EmisWb3BandCoeffs");
        if (paramstr != NULL)
        {
            sscanf(paramstr, "%lf %lf %lf %lf", &emis_wb_coeffs[0], 
                    &emis_wb_coeffs[1], &emis_wb_coeffs[2], &emis_wb_coeffs[3]);
            INFO("EMIS_WB_COEF  %.4f,%.4f,%.4f,%.4f",
                emis_wb_coeffs[0], emis_wb_coeffs[1], emis_wb_coeffs[2], 
                emis_wb_coeffs[3]);
        }
        else
        {
            ERREXIT(9, "Error reading value of EmisWb3BandCoeffs from "
                    "group RuntimeParameters in %s", run_params_filename);
        }
        paramstr = getRuntimeParameter("CoVeg3Band");
        if (paramstr != NULL)
        {
            sscanf(paramstr, "%lf %lf %lf", &co_veg[0], &co_veg[1], &co_veg[2]);
        }
        else
        {
            ERREXIT(9, "Error reading value of CoVeg3Band from "
                    "group RuntimeParameters in %s", run_params_filename);
        }

        paramstr = getRuntimeParameter("CoBare3Band");
        if (paramstr != NULL)
        {
            sscanf(paramstr, "%lf %lf %lf", &co_bare[0], &co_bare[1], &co_bare[2]);
        }
        else
        {
            ERREXIT(9, "Error reading value of CoBare3Band from "
                    "group RuntimeParameters in %s", run_params_filename);
        }
        // Get RTTOV script.
        elem = findGroupElement(run_params, "RTTOV", "Rttov3BandScript");
        if (elem == NULL)
            ERREXIT(8, "Could not find Rttov3BandScript "
            "in group RTTOV in %s", run_params_filename);
        paramstr = getValue(elem, 0);
        if (!paramstr)
            ERREXIT(9, "Error reading value of Rttov3BandScript from "
            "group RTTOV in %s", run_params_filename);
        strcpy(rttov_script, paramstr);
    }
    else
    {
        paramstr = getRuntimeParameter("EmisWb5BandCoeffs");
        if (paramstr != NULL)
        {
            sscanf(paramstr, "%lf %lf %lf %lf %lf %lf", &emis_wb_coeffs[0], 
                    &emis_wb_coeffs[1], &emis_wb_coeffs[2], &emis_wb_coeffs[3], 
                    &emis_wb_coeffs[4], &emis_wb_coeffs[5]);
            INFO("EMIS_WB_COEF  %.4f,%.4f,%.4f,%.4f,%.4f,%.4f",
                emis_wb_coeffs[0], emis_wb_coeffs[1], emis_wb_coeffs[2], 
                emis_wb_coeffs[3], emis_wb_coeffs[4], emis_wb_coeffs[5]);
        }
        else
        {
            ERREXIT(9, "Error reading value of EmisWb5BandCoeffs from "
                    "group RuntimeParameters in %s", run_params_filename);
        }
        paramstr = getRuntimeParameter("CoVeg5Band");
        if (paramstr != NULL)
        {
            sscanf(paramstr, "%lf %lf %lf", &co_veg[0], &co_veg[1], &co_veg[2]);
        }
        else
        {
            ERREXIT(9, "Error reading value of CoVeg5Band from "
                    "group RuntimeParameters in %s", run_params_filename);
        }

        paramstr = getRuntimeParameter("CoBare5Band");
        if (paramstr != NULL)
        {
            sscanf(paramstr, "%lf %lf %lf", &co_bare[0], &co_bare[1], &co_bare[2]);
        }
        else
        {
            ERREXIT(9, "Error reading value of CoBare5Band from "
                    "group RuntimeParameters in %s", run_params_filename);
        }
        // Get RTTOV script.
        elem = findGroupElement(run_params, "RTTOV", "Rttov5BandScript");
        if (elem == NULL)
            ERREXIT(8, "Could not find Rttov5BandScript "
            "in group RTTOV in %s", run_params_filename);
        paramstr = getValue(elem, 0);
        if (!paramstr)
            ERREXIT(9, "Error reading value of Rttov5BandScript from "
            "group RTTOV in %s", run_params_filename);
        strcpy(rttov_script, paramstr);
    }
    INFO("CO_VEG       = %.4f,%.4f,%.4f", co_veg[0], co_veg[1], co_veg[2]);
    INFO("CO_BARE      = %.4f,%.4f,%.4f", co_bare[0], co_bare[1], co_bare[2]);
    INFO("EMAX_VEG     = %.4f", emax_veg);
    INFO("EMAX_BARE    = %.4f", emax_bare);
    INFO("RTTOV_SCRIPT = %s", rttov_script);
    INFO("Processing %d bands:", n_channels);
    for (b = 0; b < n_channels; b++)
        INFO("    Band %d", band[b] + 1);

    // ========================================================
    // Read ECOSTRESS GEO Data
    // ========================================================
#ifdef COMMENT
    printf("Reading GEO data.\n");
    fflush(stdout);
#endif
    init_rad(&vRAD);
    tes_read_geo_data(&vRAD, GEO_filename);
    int n_lines = vRAD.Lat.size1;
    int n_pixels = vRAD.Lat.size2;
    int channel_pix = n_lines * n_pixels; // Number of pixels per channel
    size_t channel_size = channel_pix * sizeof(double); // size for memcpy

    // ========================================================
    // Read ECOSTRESS Radiance Data
    // ========================================================
#ifdef COMMENT
    printf("Reading RAD data.\n");
    fflush(stdout);
#endif
    tes_read_rad_data(&vRAD, RAD_filename, n_channels, band);

    Mat3d Y;
    mat3d_init(&Y);
    mat3d_alloc(&Y, n_channels, n_lines, n_pixels);
    // Populate Y with Rad1 ... Radn
    for (b = 0; b < n_channels; b++)
        memcpy(&Y.vals[b * channel_pix], vRAD.Rad[b].vals, channel_size);

#ifdef TESTING
    // When testing, output all data into a debug file and
    // describe the data sets
    describe_matrix("Latitude", &vRAD.Lat);
    describe_matrix("Longitude", &vRAD.Lon);
    describe_matrix("Height", &vRAD.El);
    describe_matrix("Satzen", &vRAD.Satzen);
    describe_matrix("Solzen", &vRAD.Solzen);
    describe_mat3d("Y", &Y);

    strcpy(debugfn, product_path);
    strcat(debugfn, "/ECOSTRESSdebug1.h5");
    hid_t vdbg_id = create_hdf5(debugfn);
    if (vdbg_id == FAIL)
        ERROR("Could not create %s", debugfn);
    stat = 0;
    stat += hdf5write_mat(vdbg_id, "Latitude", &vRAD.Lat);
    stat += hdf5write_mat(vdbg_id, "Longitude", &vRAD.Lon);
    stat += hdf5write_mat(vdbg_id, "Height", &vRAD.El);
    stat += hdf5write_mat(vdbg_id, "Satzen", &vRAD.Satzen);
    stat += hdf5write_mat(vdbg_id, "Solzen", &vRAD.Solzen);
    if (stat < 0)
        WARNING("%d error(s) occured while writing %s", 
                -stat, debugfn);
    close_hdf5(vdbg_id);
#endif

#ifdef DEBUG_PIXEL
    for (b=0; b<n_channels; b++)
    {
        printf("Y(%d,%d,%d)=%.4f\n", b, DEBUG_LINE, DEBUG_PIXEL, 
               mat3d_get(&Y, b, DEBUG_LINE, DEBUG_PIXEL));
    }
#endif

    // Get crop size for subsequent processing.
    CropSize crop;
    double cext = 0.0; // Crop extend size in degrees
    mat_set_fillvalues_to_nans(&vRAD.Lat, -9999.0);
    mat_set_fillvalues_to_nans(&vRAD.Lon, -9999.0);
    crop.minLat = mat_min(&vRAD.Lat) - cext; if (crop.minLat < -90.0) crop.minLat = -90.0;
    crop.maxLat = mat_max(&vRAD.Lat) + cext; if (crop.maxLat > 90.0) crop.maxLat = 90.0;
    crop.minLon = mat_min(&vRAD.Lon) - cext; if (crop.minLon < -180.0) crop.minLon = -180.0;
    crop.maxLon = mat_max(&vRAD.Lon) + cext; if (crop.maxLon > 180.0) crop.maxLon = 180.0;

#ifdef TESTING
    printf("min(Lat)=%.2f max(Lat)=%.2f\n",  mat_min(&vRAD.Lat),  mat_max(&vRAD.Lat));
    printf("min(Lon)=%.2f max(Lon)=%.2f\n",  mat_min(&vRAD.Lon),  mat_max(&vRAD.Lon));
    printf("crop: minLat=%.2f maxLat=%.2f minLon=%.2f maxLon=%.2f\n", crop.minLat, crop.maxLat, crop.minLon, crop.maxLon);
#endif

    // ========================================================
    // ASTER GED Processing
    // ========================================================
#ifdef COMMENT
    if (run_tgwvs)
    {
        printf("Loading ASTER GED data.\n");
    }
#endif
#ifdef ENABLE_TIMING
    time_t startaster = getusecs();
#endif

    // Init emis_aster to get the results of ASTER GED processing.
    Matrix emis_aster;
    mat_init(&emis_aster);
    // Skip call to read ASTER if not running TG_WVS
    if (run_tgwvs)
    {
        bool adjust_aster = false;
        // ECOSTRESS doesn't currently have snow/water index. 
        // Create zero-filled matrix.
        MatUint8 swi;
        mat_uint8_init(&swi);
        mat_uint8_alloc(&swi, vRAD.Rad[0].size1, vRAD.Rad[0].size2);
        // Call ASTER GED loader to get data for the current granule.
        int aster_ged_result = 
            read_aster_ged(ASTER_dir, &vRAD.Lat, &vRAD.Lon, &emis_aster,
            &vRAD.Water, &swi, NULL, adjust_aster, ECOSTRESS);
        if (aster_ged_result < 0)
            ERREXIT(33, "Failed to input ASTER GED data", NONE);
        mat_uint8_clear(&swi);
    }

#ifdef ENABLE_TIMING
    if (run_tgwvs)
    {
        printf("Aster GED time: %.3f seconds\n",
            toseconds(elapsed(startaster)));
    }
#endif

#ifdef DEBUG_PIXEL
    if (run_tgwvs)
        printf("emis_aster(%d,%d) = %.4f\n", 
            DEBUG_LINE, DEBUG_PIXEL, mat_get(&emis_aster, DEBUG_LINE, DEBUG_PIXEL));
#endif

#ifdef TESTING
    describe_matrix("emis_aster", &emis_aster);
    strcpy(debugfn, product_path);
    strcat(debugfn, "/ASTER_debug.h5");
    vdbg_id = create_hdf5(debugfn);
    if (vdbg_id == FAIL)
        ERROR("Could not create %s", debugfn);
    stat  = hdf5write_mat(vdbg_id, "/EMIS_ASTER", &emis_aster);
    if (stat < 0)
        WARNING("unable to write to %s", debugfn);
    close_hdf5(vdbg_id);
#endif

    // ========================================================
    // Read NWP Data
    // ========================================================
    Vector Pres;
    vec_init(&Pres);
    tes_atm_init(&nwpATM);
    char *found_ecmwf = strstr(NWP_dir, nwp_key_ecmwf);
    char *found_geos = strstr(NWP_dir, nwp_key_geos);
    char *found_merra = strstr(NWP_dir, nwp_key_merra);
    char *found_ncep = strstr(NWP_dir, nwp_key_ncep);
    if (found_merra)
    {
        // ========================================================
        // Read MERRA Data
        // ========================================================
#ifdef COMMENT
        printf("Loading MERRA2 data.\n");
        fflush(stdout);
#endif
#ifdef ENABLE_TIMING
        time_t startmerra = getusecs();
#endif

        tes_read_interp_atmos_merra(&nwpATM, &tes_date, NWP_dir);

        //---------------------------------------------
        // Clamp Merra data
        //---------------------------------------------
        clamp3d(&nwpATM.t, hard_t_min, hard_t_max);
        clamp3d(&nwpATM.q, hard_v_min, hard_v_max);
        clamp2d(&nwpATM.sp, hard_p_min, hard_p_max);

        // Invert the Pressure (Height) matrix for use by RTTOV.
        vec_alloc(&Pres, nwpATM.lev.size);
        for(i=0; i<nwpATM.lev.size; i++) 
            Pres.vals[i] = nwpATM.lev.vals[nwpATM.lev.size-i-1];

        // Create the cutout data filename for MERRA.
        sprintf(cutout_filename, "%s/%s_%05d_%03d_%s_%04d_%02d.h5",
            product_path, merra_cutout_root, iorbit, iscene, dateTtime, ibuild, ivv);

        // Designate MERRA2 NWP data 
        updateValue(&lste_metadata, "NWPSource",
                DFNT_CHAR8, strlen("MERRA2"), strdup("MERRA2"));

#ifdef ENABLE_TIMING
        printf("Merra time: %.3f seconds\n",
                toseconds(elapsed(startmerra)));
        fflush(stdout);
#endif
    }
    else if (found_geos)
    {
        // ========================================================
        // Read GEOS Data
        // ========================================================
#ifdef COMMENT
        printf("Loading GEOS data.\n");
        fflush(stdout);
#endif
#ifdef ENABLE_TIMING
        time_t startgeos = getusecs();
#endif

        tes_read_interp_atmos_geos(&nwpATM, &tes_date, NWP_dir, &crop);

        // Invert the Pressure (Height) matrix for use by RTTOV.
        vec_alloc(&Pres, nwpATM.lev.size);
        for(i=0; i<nwpATM.lev.size; i++) 
            Pres.vals[i] = nwpATM.lev.vals[nwpATM.lev.size-i-1];

        // Create the cutout data filename for GEOS.
        sprintf(cutout_filename, "%s/%s_%05d_%03d_%s_%04d_%02d.h5",
            product_path, geos_cutout_root, iorbit, iscene, dateTtime, ibuild, ivv);

        // Designate GEOS5 NWP data 
        updateValue(&lste_metadata, "NWPSource",
                DFNT_CHAR8, strlen("GEOS5"), strdup("GEOS5"));

#ifdef ENABLE_TIMING
        printf("GEOS time: %.3f seconds\n",
                toseconds(elapsed(startgeos)));
        fflush(stdout);
#endif
    }
    else if (found_ncep)
    {
        // ========================================================
        // Read NCEP Data
        // ========================================================
#ifdef COMMENT
        printf("Loading NCEP data.\n");
        fflush(stdout);
#endif
#ifdef ENABLE_TIMING
        time_t startncep = getusecs();
#endif

        TesATM ncepATM;
        tes_atm_init(&ncepATM);
        tes_read_interp_atmos_ncep(&ncepATM, &tes_date, NWP_dir);

// Expand data from NCEP to double size
#ifdef COMMENT
        printf("Expanding NCEP data.\n");
#endif
        // Alloc the output matrix sizes to double the NCEP results.
        vec_alloc(&nwpATM.lat, ncepATM.lat.size * 2);
        vec_alloc(&nwpATM.lon, ncepATM.lon.size * 2);
        vec_alloc(&nwpATM.lev, ncepATM.lev.size);
        mat3d_alloc(&nwpATM.t, nwpATM.lev.size, nwpATM.lat.size, nwpATM.lon.size);
        mat3d_alloc(&nwpATM.q, nwpATM.lev.size, nwpATM.lat.size, nwpATM.lon.size);
        mat_alloc(&nwpATM.sp, nwpATM.lat.size, nwpATM.lon.size);
        mat_alloc(&nwpATM.skt, nwpATM.lat.size, nwpATM.lon.size);
        mat_alloc(&nwpATM.t2, nwpATM.lat.size, nwpATM.lon.size);
        mat_alloc(&nwpATM.q2, nwpATM.lat.size, nwpATM.lon.size);
        mat_alloc(&nwpATM.tcw, nwpATM.lat.size, nwpATM.lon.size);

        // Create new Lat/Lon for doubled arrays.
        int inwp, jnwp;
        int inwp1, jnwp1;
        int incep, jncep, k;
        double nv0;
        // Fill in the latitude values.
        for (incep = 0; incep < ncepATM.lat.size; incep++)
        {
          inwp = incep * 2;
          inwp1 = inwp + 1; 
          if (inwp1 >= nwpATM.lat.size)
              inwp1 = inwp;
          nv0 = vec_get(&ncepATM.lat, incep);
          vec_set(&nwpATM.lat, inwp, nv0);
          vec_set(&nwpATM.lat, inwp1, nv0 + 0.5);
        }

        // Fill in the longitude values.
        for (jncep = 0; jncep < ncepATM.lon.size; jncep++)
        {
          jnwp = jncep * 2;
          jnwp1 = jnwp + 1; 
          if (jnwp1 >= nwpATM.lon.size)
              jnwp1 = jnwp;
          nv0 = vec_get(&ncepATM.lon, jncep);
          vec_set(&nwpATM.lon, jnwp, nv0);
          vec_set(&nwpATM.lon, jnwp+1, nv0 + 0.5);
        }

        // Copy the pressule level values.
        for (k = 0; k < ncepATM.lev.size; k++)
            nwpATM.lev.vals[k] = ncepATM.lev.vals[k];

        // Create lat/lon meshgrids for the NCEP inputs.
        Matrix ncepsamplemeshLat;
        mat_init(&ncepsamplemeshLat);
        mat_alloc(&ncepsamplemeshLat, ncepATM.lat.size, ncepATM.lon.size);
        Matrix ncepsamplemeshLon;
        mat_init(&ncepsamplemeshLon);
        mat_alloc(&ncepsamplemeshLon, ncepATM.lat.size, ncepATM.lon.size);
        int inceprow, incepcol;
        for (inceprow = 0; inceprow < ncepATM.lat.size; inceprow++)
        {
            for (incepcol = 0; incepcol < ncepATM.lon.size; incepcol++)
            {
                mat_set(&ncepsamplemeshLat, inceprow, incepcol, 
                        vec_get(&ncepATM.lat, inceprow));
                mat_set(&ncepsamplemeshLon, inceprow, incepcol, 
                        vec_get(&ncepATM.lon, incepcol));
            }
        }

        // Create lat/lon meshgrids for the NWP outputs.
        Matrix nwpquerymeshLat;
        mat_init(&nwpquerymeshLat);
        mat_alloc(&nwpquerymeshLat, nwpATM.lat.size, nwpATM.lon.size);
        Matrix nwpquerymeshLon;
        mat_init(&nwpquerymeshLon);
        mat_alloc(&nwpquerymeshLon, nwpATM.lat.size, nwpATM.lon.size);
        int inwprow, inwpcol;
        for (inwprow = 0; inwprow < nwpATM.lat.size; inwprow++)
        {
            for (inwpcol = 0; inwpcol < nwpATM.lon.size; inwpcol++)
            {
                mat_set(&nwpquerymeshLat, inwprow, inwpcol, vec_get(&nwpATM.lat, inwprow));
                mat_set(&nwpquerymeshLon, inwprow, inwpcol, vec_get(&nwpATM.lon, inwpcol));
            }
        }

        // Set array of data pointers for 2D data to be interpolated.
        // Treat each pressure level in T and QV as a separate 2D dataset.
        int nsets = ncepATM.lev.size * 2 + 5;
        double *sampsets[nsets];
        double *outsets[nsets];
        int iset = 0;
        for (k = 0; k < ncepATM.lev.size; k++)
        {
            sampsets[iset] = &ncepATM.t.vals[k * ncepATM.t.size2 * ncepATM.t.size3];
            outsets[iset] = &nwpATM.t.vals[k * nwpATM.t.size2 * nwpATM.t.size3];
            iset++;
            sampsets[iset] = &ncepATM.q.vals[k * ncepATM.q.size2 * ncepATM.q.size3];
            outsets[iset] = &nwpATM.q.vals[k * nwpATM.q.size2 * nwpATM.q.size3];
            iset++;
        }
        sampsets[iset] = ncepATM.sp.vals;
        outsets[iset] = nwpATM.sp.vals;
        iset++;
        sampsets[iset] = ncepATM.skt.vals;
        outsets[iset] = nwpATM.skt.vals;
        iset++;
        sampsets[iset] = ncepATM.t2.vals;
        outsets[iset] = nwpATM.t2.vals;
        iset++;
        sampsets[iset] = ncepATM.q2.vals;
        outsets[iset] = nwpATM.q2.vals;
        iset++;
        sampsets[iset] = ncepATM.tcw.vals;
        outsets[iset] = nwpATM.tcw.vals;

        // Interpolate the NCEP data onto the expanded NWP grid.
        multi_interp2(
                ncepsamplemeshLat.vals, ncepsamplemeshLon.vals, 
                ncepATM.lat.size, ncepATM.lon.size, 
                nwpquerymeshLat.vals, nwpquerymeshLon.vals,
                nwpATM.lat.size * nwpATM.lon.size,
                sampsets, outsets, nsets);

        mat_clear(&ncepsamplemeshLat);
        mat_clear(&ncepsamplemeshLon);
        mat_clear(&nwpquerymeshLat);
        mat_clear(&nwpquerymeshLon);

#ifdef TESTING
        describe_mat3d("ATM.t", &nwpATM.t);
        describe_mat3d("ATM.q", &nwpATM.q);
        describe_matrix("ATM.sp", &nwpATM.sp);
#endif

        //---------------------------------------------
        // Clamp NCEP data
        //---------------------------------------------
        clamp3d(&nwpATM.t, hard_t_min, hard_t_max);
        clamp3d(&nwpATM.q, hard_v_min, hard_v_max);
        clamp2d(&nwpATM.sp, hard_p_min, hard_p_max);
        clamp2d(&nwpATM.skt, hard_t_min, hard_t_max);
        clamp2d(&nwpATM.t2, hard_t_min, hard_t_max);
        clamp2d(&nwpATM.q2, hard_v_min, hard_v_max);

        // Copy the Pressure (Height) matrix for use by RTTOV.
        vec_alloc(&Pres, nwpATM.lev.size);
        for(i=0; i<nwpATM.lev.size; i++) 
            Pres.vals[i] = nwpATM.lev.vals[i];

        // Create the cutout data filename for NCEP.
        sprintf(cutout_filename, "%s/%s_%05d_%03d_%s_%04d_%02d.h5",
            product_path, ncep_cutout_root, iorbit, iscene, dateTtime, ibuild, ivv);

        // Designate NCEP NWP data 
        updateValue(&lste_metadata, "NWPSource",
                DFNT_CHAR8, strlen("NCEP"), strdup("NCEP"));

#ifdef ENABLE_TIMING
        printf("NCEP time: %.3f seconds\n",
                toseconds(elapsed(startncep)));
        fflush(stdout);
#endif
#ifdef TESTING
        describe_mat3d("ATM.t", &nwpATM.t);
        describe_mat3d("ATM.q", &nwpATM.q);
        describe_matrix("ATM.sp", &nwpATM.sp);
        describe_matrix("ATM.tcw", &nwpATM.tcw);
#endif
    }
    else if (found_ecmwf)
    {
        // ========================================================
        // Read ECMWF Data
        // ========================================================
        ERROR("ECMWF processing is not implemented yet.", NONE);
#ifdef COMMENT
        printf("Loading ECMWF data.\n");
        fflush(stdout);
#endif
#ifdef ENABLE_TIMING
        time_t startecmwf = getusecs();
#endif

        // TODO: tes_read_interp_atmos_ecmwf(&nwpATM, &tes_date, NWP_dir);

        // Create the cutout data filename for NCEP.
        sprintf(cutout_filename, "%s/%s_%05d_%03d_%s_%04d_%02d.h5",
            product_path, ecmwf_cutout_root, iorbit, iscene, dateTtime, ibuild, ivv);

        // Copy the Pressure (Height) matrix for use by RTTOV.
        vec_alloc(&Pres, nwpATM.lev.size);
        for(i=0; i<nwpATM.lev.size; i++) 
            Pres.vals[i] = nwpATM.lev.vals[i];

        // Designate ECMWF NWP data 
        updateValue(&lste_metadata, "NWPSource",
                DFNT_CHAR8, strlen("ECMWF"), strdup("ECMWF"));

#ifdef ENABLE_TIMING
        printf("ECMWF time: %.3f seconds\n",
                toseconds(elapsed(startecmwf)));
        fflush(stdout);
#endif
    }
    else
    {
        ERREXIT(36, "%s does not designate an NWP input type. "
            "Expected %s, %s, %s, or %s in the name.",
            NWP_dir, nwp_key_merra, nwp_key_ncep, nwp_key_ecmwf, nwp_key_geos);
    }
 
    // ========================================================
    // RTTOV Input Profile Data Preparation
    // ========================================================
#ifdef COMMENT
    printf("Mapping data for RTTOV.\n");
    fflush(stdout);
#endif
#ifdef ENABLE_TIMING
    time_t startrttov = getusecs();
#endif

    // If the NWP input did not provide TCW, create the TCW matrix 
    // from the q data.
    if (nwpATM.tcw.size1 == 0)
    {
        int px, py, pz;
        double cumsum, dp, dq, q0, q1, sp_val;
        mat_alloc(&nwpATM.tcw, nwpATM.q.size2, nwpATM.q.size3);

        const double PPMV_to_g_per_kg = 1.0/(1000.0*(28.966/18.015));

        for (px = 0; px < nwpATM.q.size2; px++)
        {
            for (py = 0; py < nwpATM.q.size3; py++)
            {
                sp_val = mat_get(&nwpATM.sp, px, py);
                
                cumsum = 0;
                for (pz = 0; pz < nwpATM.q.size1 - 1; pz++)
                {
                    // Filter out levels where Pressure >= surface pressure.
                    if (Pres.vals[pz] >= sp_val)
                        continue;
                    if (Pres.vals[pz+1] >= sp_val)
                        continue;
                    q0 = mat3d_get(&nwpATM.q, pz, px, py);
                    q1 = mat3d_get(&nwpATM.q, pz+1, px, py);
                    if (is_nan(q0) || is_nan(q1))
                    {
                        // If NaN, do not add anything into cumsum
                        continue;
                    }
                    dq = (q0 + q1) / 2.0;
                    dp = Pres.vals[pz+1] - Pres.vals[pz]; // diff between height pressures
                    cumsum += dq*PPMV_to_g_per_kg * dp;
                }
                mat_set(&nwpATM.tcw, px, py, cumsum / 100.0 / 9.8);
            }
        }
    }

#ifdef DEBUG_PIXEL
    printf("DEBUG -- NWP Data (%d,%d): sp=%.4f, tcw=%.4f\n",
           DEBUG_LINE, DEBUG_PIXEL, 
           mat_get(&nwpATM.sp, DEBUG_LINE, DEBUG_PIXEL),
           mat_get(&nwpATM.tcw, DEBUG_LINE, DEBUG_PIXEL));
#endif

    // Determine the Lat/Lon range for the granule.
    double minLat = crop.minLat;
    double maxLat = crop.maxLat;
    double minLon = crop.minLon;
    double maxLon = crop.maxLon;

    // Create lat/lon meshgrids for Merra data
    int llrow, llcol, lllev, llix;
    int startLat;
    int startLon;
    Matrix meshgridLat;
    mat_init(&meshgridLat);
    mat_alloc(&meshgridLat, nwpATM.lat.size, nwpATM.lon.size);
    Matrix meshgridLon;
    mat_init(&meshgridLon);
    mat_alloc(&meshgridLon, nwpATM.lat.size, nwpATM.lon.size);
    for (llrow = 0; llrow < nwpATM.lat.size; llrow++)
    {
        for (llcol = 0; llcol < nwpATM.lon.size; llcol++)
        {
            mat_set(&meshgridLat, llrow, llcol, vec_get(&nwpATM.lat, llrow));
            mat_set(&meshgridLon, llrow, llcol, vec_get(&nwpATM.lon, llcol));
        }
    }

    // Crop the NWP data to a subgrid that overlays the granule.
    GridMapData griddata;
    double *insets[2];  // Input datasets to map
    double *outsets[2]; // Output of mapping
    int nLats, nLons;
    Matrix cropLat;
    mat_init(&cropLat);
    Matrix cropLon;
    mat_init(&cropLon);

    // GEOS data is already cropped.
    if (found_geos)
    {
        startLat = 0;
        startLon = 0;
        nLats = nwpATM.lat.size;
        nLons = nwpATM.lon.size;
        mat_alloc(&cropLat, nLats, nLons);
        mat_alloc(&cropLon, nLats, nLons);
        griddata.lat = cropLat.vals;     // Query grid latitudes
        griddata.lon = cropLon.vals;     // Query grid longitudes
        griddata.nlines = nwpATM.lat.size;  // n lines in query data
        griddata.ncols = nwpATM.lon.size;   // n pixels per line in query data
        griddata.inrows = vRAD.Lat.size1;   // NWP data dimension: nlines
        griddata.incols = vRAD.Lat.size2;   // NWP data dimension: ncols
    }
    else // not GEOS
    {
        const double latextend = 2.0;
        const double lonextend = 2.0;
        startLat = nwpATM.lat.size;
        startLon = nwpATM.lon.size;
        int endLat = 0;
        int endLon = 0;
        double llval;
        for (llix = 0; llix < nwpATM.lat.size; llix++)
        {
            llval = vec_get(&nwpATM.lat, llix);
            if (llval >= (minLat - latextend) && llval <= (maxLat + latextend))
            {
                if (llix < startLat)
                {
                    startLat = llix;
                }
                if (llix > endLat)
                {
                    endLat = llix;
                }
            }
        }

        for (llix = 0; llix < nwpATM.lon.size; llix++)
        {
            llval = vec_get(&nwpATM.lon, llix);
            if (llval >= (minLon - lonextend) && llval <= (maxLon + lonextend))
            {
                if (llix < startLon)
                {
                    startLon = llix;
                }
                if (llix > endLon)
                {
                    endLon = llix;
                }
            }
        }

        endLat++; // Point end values 1 beyond last usable entry.
        endLon++;
        nLats = endLat - startLat;
        nLons = endLon - startLon;


        // Map granule data onto Merra cropped grid.
        //   vRAD.Satzen
        //   Hsurf
        mat_alloc(&cropLat, nLats, nLons);
        mat_alloc(&cropLon, nLats, nLons);
        griddata.lat = cropLat.vals;     // Query grid latitudes
        griddata.lon = cropLon.vals;     // Query grid longitudes
        griddata.nlines = nLats;         // n lines in query data
        griddata.ncols = nLons;          // n pixels per line in query data
        griddata.inrows = vRAD.Lat.size1;  // Merra data dimension: nlines
        griddata.incols = vRAD.Lat.size2;  // Merra data dimension: ncols
    }

    // Transfer NWP data to cropped arrays
    Mat3d Pmesh; // Pmesh is for the cutout file.
    mat3d_init(&Pmesh);
    mat3d_alloc(&Pmesh, nwpATM.t.size1, nLats, nLons);
    Mat3d cropT;
    mat3d_init(&cropT);
    mat3d_alloc(&cropT, nwpATM.t.size1, nLats, nLons);
    Mat3d cropQ;
    mat3d_init(&cropQ);
    mat3d_alloc(&cropQ, nwpATM.q.size1, nLats, nLons);
    Matrix cropSP;
    mat_init(&cropSP);
    mat_alloc(&cropSP, nLats, nLons);
    Matrix cropTCW;
    mat_init(&cropTCW);
    mat_alloc(&cropTCW, nLats, nLons);
    Matrix cropskt;
    mat_init(&cropskt);
    if (nwpATM.skt.size1 != 0)
        mat_alloc(&cropskt, nLats, nLons);
    Matrix cropt2;
    mat_init(&cropt2);
    if (nwpATM.t2.size1 != 0)
        mat_alloc(&cropt2, nLats, nLons);
    Matrix cropq2;
    mat_init(&cropq2);
    if (nwpATM.q2.size1 != 0)
        mat_alloc(&cropq2, nLats, nLons);

    for (llrow = 0; llrow < nLats; llrow++)
    {
        for (llcol = 0; llcol < nLons; llcol++)
        {
            // Crop 2D matrices
            mat_set(&cropLat, llrow, llcol,
                mat_get(&meshgridLat, startLat + llrow, startLon + llcol));
            mat_set(&cropLon, llrow, llcol,
                mat_get(&meshgridLon, startLat + llrow, startLon + llcol));
            mat_set(&cropSP, llrow, llcol, 
                mat_get(&nwpATM.sp, startLat + llrow, startLon + llcol));
            mat_set(&cropTCW, llrow, llcol, 
                mat_get(&nwpATM.tcw, startLat + llrow, startLon + llcol));
            if (cropskt.size1 != 0)
                mat_set(&cropskt, llrow, llcol, 
                    mat_get(&nwpATM.skt, startLat + llrow, startLon + llcol));
            if (cropt2.size1 != 0)
                mat_set(&cropt2, llrow, llcol, 
                    mat_get(&nwpATM.t2, startLat + llrow, startLon + llcol));
            if (cropq2.size1 != 0)
                mat_set(&cropq2, llrow, llcol, 
                    mat_get(&nwpATM.q2, startLat + llrow, startLon + llcol));
            // Crop 3D matrices
            for (lllev = 0; lllev < nwpATM.t.size1; lllev++)
            {
                mat3d_set(&cropT, lllev, llrow, llcol, 
                    mat3d_get(&nwpATM.t, lllev, startLat + llrow, startLon + llcol));
                mat3d_set(&cropQ, lllev, llrow, llcol, 
                    mat3d_get(&nwpATM.q, lllev, startLat + llrow, startLon + llcol));
                mat3d_set(&Pmesh, lllev, llrow, llcol, Pres.vals[lllev]);
            }
        }
    }
    Matrix cropHsurf;
    mat_init(&cropHsurf);
    mat_alloc(&cropHsurf, nLats, nLons);
    Matrix cropSatZen;
    mat_init(&cropSatZen);
    mat_alloc(&cropSatZen, nLats, nLons);
    insets[0] = vRAD.Satzen.vals;
    insets[1] = vRAD.El.vals;
    outsets[0] = cropSatZen.vals;
    outsets[1] = cropHsurf.vals;
    griddata.nsets = 2;                 // Each level acts like a separate set
    griddata.dsetin = insets;
    griddata.dsetout = outsets;

#ifdef TESTING
    time_t startmap = getusecs();
#endif

    mapnearest(&griddata, vRAD.Lat.vals, vRAD.Lon.vals);

#ifdef TESTING
    printf("mapnearest time: %.3f seconds\n",
            toseconds(elapsed(startmap)));
    fflush(stdout);
#endif

    // Get lowest level of T for Tsurf
    Matrix cropTsurf;
    mat_init(&cropTsurf);
    mat_alloc(&cropTsurf, nLats, nLons);
    Matrix cropQsurf;
    mat_init(&cropQsurf);
    mat_alloc(&cropQsurf, nLats, nLons);
    Matrix cropST;
    mat_init(&cropST);
    mat_alloc(&cropST, nLats, nLons);
    int tx, ty, tz;
    tz = cropT.size1 - 1;
    for (tx = 0; tx < nLats; tx++)
    {
        for (ty = 0; ty < nLons; ty++)
        {
            mat_set(&cropTsurf, tx, ty, 
                mat3d_get(&cropT, tz, tx, ty));
            mat_set(&cropQsurf, tx, ty, 
                mat3d_get(&cropQ, tz, tx, ty));
        }
    }

    // Select Tsurf matrix if skt is not available.
    Matrix *Tskin = &cropskt;
    if (cropskt.size1 == 0)
        Tskin = &cropTsurf;

    // Create Bemis for output, filled with 1e-6.
    Mat3d Bemis;
    mat3d_init(&Bemis);
    mat3d_alloc(&Bemis, n_channels, nLats, nLons);
    for (b = 0; b < n_channels * nLats * nLons; b++)
        Bemis.vals[b] = 1e-6;

#ifdef TESTING
    describe_matrix("cropSatZen", &cropSatZen);
    describe_matrix("cropHsurf", &cropHsurf);
    describe_matrix("cropQsurf", &cropQsurf);
    describe_matrix("cropskt", &cropskt);
    describe_matrix("cropt2", &cropt2);
    describe_matrix("cropq2", &cropq2);
    describe_matrix("cropTCW", &cropTCW);
#endif
    // ========================================================
    // Create cutout file to contain RTTOV inputs
    // ========================================================
#ifdef CREATE_CUTOUT
#ifdef COMMENT
    printf("Creating RTTOV cutout file %s\n", cutout_filename);
    fflush(stdout);
#endif
#ifdef ENABLE_TIMING
    time_t startcutout = getusecs();
#endif
    // Write the cutout file
    hid_t cutout = create_hdf5(cutout_filename);
    if (cutout == -1)
        ERREXIT(12, "Unable to create file %s", cutout_filename);
    hid_t cutout_geo_group = H5Gcreate2(cutout, "/GEO", 
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (cutout_geo_group == -1)
        ERREXIT(13, "Could not create group /GEO in %s", cutout_filename);

    hid_t cutout_data_group = H5Gcreate2(cutout, "/DATA", 
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (cutout_data_group == -1)
        ERREXIT(13, "Could not create group /DATA in %s", cutout_filename);

    stat = hdf5write_mat3d(cutout, "/GEO/P", &Pmesh);
    if (stat == -1)
        ERREXIT(13, "Unable to write dataset P to %s", cutout_filename);
    double minP = 0.0;
    double maxP = 1100.0;
    stat = writeDatasetMetadataHdf5(cutout, "/GEO/P",
            "Pressure Meshgrid", "hPa", 
            DFNT_FLOAT64, &minP, &maxP, 
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(13, "Failed to write attributes to dataset /GEO/P in %s",
            output_filename);

    stat = hdf5write_mat3d(cutout, "/DATA/T", &cropT);
    if (stat == -1)
        ERREXIT(13, "Unable to write dataset T to %s", cutout_filename);
    double minT = 150.0;
    double maxT = 350.0;
    stat = writeDatasetMetadataHdf5(cutout, "/DATA/T",
            "Temperature", "K", 
            DFNT_FLOAT64, &minT, &maxT, 
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        WARNING("Failed to write attributes to dataset /DATA/T in %s",
            output_filename);

    double minQV = 0.0;
    double maxQV = 1e9;
    stat = hdf5write_mat3d(cutout, "/DATA/QV", &cropQ);
    if (stat == -1)
        ERREXIT(13, "Unable to write dataset QV to %s", cutout_filename);
    stat = writeDatasetMetadataHdf5(cutout, "/DATA/QV",
            "Water Vapor", "ppmv", 
            DFNT_FLOAT64, &minQV, &maxQV,
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        WARNING("Failed to write attributes to dataset /DATA/QV in %s",
            output_filename);

    stat = hdf5write_mat(cutout, "/DATA/Psurf", &cropSP);
    if (stat == -1)
        ERREXIT(13, "Unable to write dataset Psurf to %s", cutout_filename);
    stat = writeDatasetMetadataHdf5(cutout, "/DATA/Psurf",
            "Surface Pressure", "hPa", 
            DFNT_FLOAT64, &minP, &maxP,
            DFNT_NONE, "n/a", 1.0, 0.0);
    if (stat < 0)
        WARNING("Failed to write attributes to dataset /DATA/Psurf in %s",
            output_filename);

    stat = hdf5write_mat(cutout, "/DATA/Qsurf", &cropQsurf);
    if (stat == -1)
        ERREXIT(13, "Unable to write dataset Qsurf to %s", cutout_filename);
    stat = writeDatasetMetadataHdf5(cutout, "/DATA/Qsurf",
            "Surface Water Vapor", "ppmv", 
            DFNT_FLOAT64, &minQV, &maxQV,
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        WARNING("Failed to write attributes to dataset /DATA/Qsurf in %s",
            output_filename);

    stat = hdf5write_mat(cutout, "/DATA/Tskin", Tskin);
    if (stat == -1)
        ERREXIT(13, "Unable to write dataset Tskin to %s", cutout_filename);
    stat = writeDatasetMetadataHdf5(cutout, "/DATA/Tskin",
            "Skin Temperature", "K", 
            DFNT_FLOAT64, &minT, &maxT,
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        WARNING("Failed to write attributes to dataset /DATA/Tskin in %s",
            output_filename);

    if (cropt2.size1 != 0)
    {
        stat = hdf5write_mat(cutout, "/DATA/T2", &cropt2);
        if (stat == -1)
            ERREXIT(13, "Unable to write dataset T2 to %s", cutout_filename);
        stat = writeDatasetMetadataHdf5(cutout, "/DATA/T2",
                "2 Meter Temperature", "K", 
                DFNT_FLOAT64, &minT, &maxT,
                DFNT_NONE, NULL, 1.0, 0.0);
        if (stat < 0)
            WARNING("Failed to write attributes to dataset /DATA/T2 in %s",
                output_filename);
    }

    if (cropq2.size1 != 0)
    {
        stat = hdf5write_mat(cutout, "/DATA/Q2", &cropq2);
        if (stat == -1)
            ERREXIT(13, "Unable to write dataset Q2 to %s", cutout_filename);
        stat = writeDatasetMetadataHdf5(cutout, "/DATA/Q2",
                "2 Meter Water Vapor", "ppmv", 
                DFNT_FLOAT64, &minQV, &maxQV,
                DFNT_NONE, NULL, 1.0, 0.0);
        if (stat < 0)
            WARNING("Failed to write attributes to dataset /DATA/Q2 in %s",
                output_filename);
    }

    stat = hdf5write_mat(cutout, "/DATA/ST", &cropST);
    if (stat == -1)
        ERREXIT(13, "Unable to write dataset ST to %s", cutout_filename);
    double minST = 0;
    double maxST = 1;
    stat = writeDatasetMetadataHdf5(cutout, "/DATA/ST",
            "Surface Type", "0=land 1=water",
            DFNT_FLOAT64, &minST, &maxST,
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        WARNING("Failed to write attributes to dataset /DATA/ST in %s",
            output_filename);

    stat = hdf5write_mat(cutout, "/DATA/Hsurf", &cropHsurf);
    if (stat == -1)
        ERREXIT(13, "Unable to write dataset Hsurf to %s", cutout_filename);
    double minElev = -500.0;
    double maxElev = 9000.0;
    stat = writeDatasetMetadataHdf5(cutout, "/DATA/Hsurf",
            "Surface Elevation", "km", 
            DFNT_FLOAT64, &minElev, &maxElev,
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        WARNING("Failed to write attributes to dataset /DATA/Hsurf in %s",
            output_filename);

    stat = hdf5write_mat(cutout, "/DATA/SatZen", &cropSatZen);
    if (stat == -1)
        ERREXIT(13, "Unable to write dataset SatZen to %s", cutout_filename);
    double minZen = 0.0;
    double maxZen = 90.0;
    stat = writeDatasetMetadataHdf5(cutout, "/DATA/SatZen",
            "Satellite Zenith", "degrees", 
            DFNT_FLOAT64, &minZen, &maxZen,
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        WARNING("Failed to write attributes to dataset /DATA/SatZen in %s",
            output_filename);

    stat = hdf5write_mat(cutout, "/GEO/Lat", &cropLat);
    if (stat == -1)
        ERREXIT(13, "Unable to write dataset Lat to %s", cutout_filename);
    double mincutlat = -90.0;
    double maxcutlat = 90.0;
    stat = writeDatasetMetadataHdf5(cutout, "/GEO/Lat",
            "Latitude", "degrees", 
            DFNT_FLOAT64, &mincutlat, &maxcutlat,
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        WARNING("Failed to write attributes to dataset /DATA/Lat in %s",
            output_filename);

    stat = hdf5write_mat(cutout, "/GEO/Lon", &cropLon);
    if (stat == -1)
        ERREXIT(13, "Unable to write dataset Lon to %s", cutout_filename);
    double mincutlon = -180.0;
    double maxcutlon = 180.0;
    stat = writeDatasetMetadataHdf5(cutout, "/GEO/Lon",
            "Longitude", "degrees",
            DFNT_FLOAT64, &mincutlon, &maxcutlon,
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        WARNING("Failed to write attributes to dataset /GEO/Lon in %s",
            output_filename);

    stat = hdf5write_mat3d(cutout, "/DATA/Bemis", &Bemis);
    if (stat == -1)
        ERREXIT(13, "Unable to write dataset Bemis to %s", cutout_filename);
    double minB = 0.0;
    double maxB = 1.0;
    stat = writeDatasetMetadataHdf5(cutout, "/DATA/Bemis",
            "Emissivity", "ratio",
            DFNT_FLOAT64, &minB, &maxB,
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        WARNING("Failed to write attributes to dataset /DATA/Bemis in %s",
            output_filename);
    close_hdf5(cutout);
#ifdef ENABLE_TIMING
    printf("Cutout file creation time: %.3f seconds\n",
            toseconds(elapsed(startcutout)));
#endif

#endif

    // ========================================================
    // Run RTTOV Simulation
    // ========================================================
#ifdef COMMENT
    printf("Running RTTOV.\n");
    fflush(stdout);
#endif

    // Initialize all the structs for RTTOV processing
    ATM_rttov rATM;
    init_ATM_rttov(&rATM, nwpATM.t.size1, nLats, nLons);
    set_rttov_atmos(&rATM, &cropQ, &cropT, Pres.vals, &cropSatZen, 
            &cropSP, &cropHsurf, &cropLat, &cropLon, &cropTCW,
            Tskin, &cropt2, &cropq2);

    RAD_rttov rRAD;
    set_RAD_rttov(&rRAD, vRAD.Lat.vals, vRAD.Lon.vals, vRAD.Satzen.vals, 
        vRAD.Lat.size1, vRAD.Lat.size2);

    RTM_rttov RTM1, RTM2;
    init_RTM_rttov(&RTM1, vRAD.Lat.size1, vRAD.Lat.size2); // RTM1 for the 1st RTTOV run

    init_RTM_rttov(&RTM2, vRAD.Lat.size1, vRAD.Lat.size2); // RTM2 for the 2nd RTTOV run

    //---------------------------------------------
    // 1st RTTOV run [SL]
    //---------------------------------------------

    int wvs_case = 0;    // reset WVS case for the 1st RTTOV run [SL]

    char fpath_rttov_profile[FILENAME_MAX];
    strcpy(fpath_rttov_profile, RTTOV_dir);
    strcat(fpath_rttov_profile, "/prof_in.bin");
    write_rttov_atmos(fpath_rttov_profile, &rATM, wvs_case);

    run_rttov_script(rttov_script, rttov_exe, rttov_coef_filename);

    char fpath_radout_dat[FILENAME_MAX];
    strcpy(fpath_radout_dat, RTTOV_dir);
    strcat(fpath_radout_dat, "/rad_out.dat");
    read_interp_rttov(fpath_radout_dat, &RTM1, &rRAD, &rATM);

    // Convert PWV to a Matrix.
    Matrix pwv;
    mat_init(&pwv);
    mat_alloc(&pwv, RTM1.kmax, RTM1.jmax);
    size_t n_pwv_bytes = RTM1.kmax * RTM1.jmax * sizeof(double);
    memcpy(pwv.vals, RTM1.PWV, n_pwv_bytes);

#ifdef DEBUG_PIXEL
    printf("PWV(%d,%d) = %.4f\n", 
        DEBUG_LINE, DEBUG_PIXEL, mat_get(&pwv, DEBUG_LINE, DEBUG_PIXEL));
#endif

#ifdef TESTING
    strcpy(debugfn, product_path);
    strcat(debugfn, "/ECOSTRESSdebug2.h5");
    vdbg_id = create_hdf5(debugfn);
    if (vdbg_id == FAIL)
        ERROR("Could not create %s", debugfn);

    describe_mat3d("NWP T", &nwpATM.t);
    describe_mat3d("NWP QV", &nwpATM.q);
    describe_matrix("NWP SP", &nwpATM.sp);
    describe_vector("NWP Pressure", &Pres);
    describe_matrix("NWP PWV", &nwpATM.tcw);
    stat  = hdf5write_mat3d(vdbg_id, "NWP_T", &nwpATM.t);
    stat += hdf5write_mat3d(vdbg_id, "NWP_QV", &nwpATM.q);
    stat += hdf5write_mat(vdbg_id, "NWP_SP", &nwpATM.sp);
    stat += hdf5write_vec(vdbg_id, "NWP_Pressure", &Pres);
    stat += hdf5write_mat(vdbg_id, "NWP_PWV", &nwpATM.tcw);
    //stat += hdfwrite_mat(vdbg_id, "cropSatZen", &cropSatZen);
    if (stat < 0)
        WARNING("unable to write to %s", debugfn);
    close_hdf5(vdbg_id);
#endif

    //---------------------------------------------    
    // 2nd RTTOV run [SL]
    //---------------------------------------------
    wvs_case = 1;    // reset WVS case for the 2nd RTTOV run [SL]

    write_rttov_atmos(fpath_rttov_profile, &rATM, wvs_case);

    run_rttov_script(rttov_script, rttov_exe, rttov_coef_filename);

    read_interp_rttov(fpath_radout_dat, &RTM2, &rRAD, &rATM);

    tes_atm_clear(&nwpATM);

    // ========================================================
    //  Process RTTOV Radiance Data
    // ========================================================
#ifdef COMMENT
    printf("Processing RTTOV results.\n");
    fflush(stdout);
#endif
    
    Mat3d t1r, t2r, pathr, skyr;
    mat3d_init(&t1r);
    mat3d_alloc(&t1r, n_channels, RTM1.kmax, RTM1.jmax);
    mat3d_init(&t2r);
    mat3d_alloc(&t2r, n_channels, RTM1.kmax, RTM1.jmax);
    mat3d_init(&pathr);
    mat3d_alloc(&pathr, n_channels, RTM1.kmax, RTM1.jmax);
    mat3d_init(&skyr);
    mat3d_alloc(&skyr, n_channels, RTM1.kmax, RTM1.jmax);
    double t1_val, t2_val, path_val;

    int valpix = 0;

    for (k = 0; k < RTM1.kmax; k++)
    {
        for (j = 0; j < RTM1.jmax; j++)
        {
            int kj = k * RTM1.jmax + j;
            for (b = 0; b < n_channels; b++)
            {
                if (RTM1.Trans_rt[b][kj] > 0.0) valpix++;
                mat3d_set(&t1r, b, k, j, RTM1.Trans_rt[b][kj]);
                mat3d_set(&t2r, b, k, j, RTM2.Trans_rt[b][kj]);
                mat3d_set(&pathr, b, k, j, RTM1.RadUp_rt[b][kj]);
                mat3d_set(&skyr, b, k, j, RTM1.RadRefDn_rt[b][kj]);
#ifdef DEBUG_RTTOV
                if (RTM1.RadUp_rt[b][kj] > 10.0)
                {
                    printf("pathr(%d,%d,%d)=RTM1.RadUp_rt[%d][%d] = %.4f\n",
                           b, k, j, b, kj, RTM1.RadUp_rt[b][kj]);
                }
                if (RTM1.RadRefDn_rt[b][kj] > 10.0)
                {
                    printf("skyr(%d,%d,%d)=RTM1.RadRefDn_rt[%d][%d] = %.4f\n",
                           b, k, j, b, kj, RTM1.RadRefDn_rt[b][kj]);
                }
#endif
            }
        }
    }

#ifdef ENABLE_TIMING
    printf("RTTOV time: %.3f seconds\n",
            toseconds(elapsed(startrttov)));
#endif

#ifdef TESTING
    describe_mat3d("t1r", &t1r);
    describe_mat3d("t2r", &t2r);
    describe_mat3d("pathr", &pathr);
    describe_mat3d("skyr", &skyr);

    strcpy(debugfn, product_path);
    strcat(debugfn, "/ECOSTRESS_rttov_debug.h5");
    vdbg_id = create_hdf5(debugfn);
    if (vdbg_id == FAIL)
        ERROR("Could not create %s", debugfn);

    stat = hdf5write_mat(vdbg_id, "pwv", &pwv);
    if (stat < 0)
        WARNING("Couldn't write pwv to %s", debugfn);

    stat = hdf5write_mat3d(vdbg_id, "t1r", &t1r);
    if (stat < 0)
        WARNING("Couldn't write t1r to %s", debugfn);

    stat = hdf5write_mat3d(vdbg_id, "t2r", &t2r);
    if (stat < 0)
        WARNING("Couldn't write t2r to %s", debugfn);

    stat = hdf5write_mat3d(vdbg_id, "pathr", &pathr);
    if (stat < 0)
        WARNING("Couldn't write pathr to %s", debugfn);

    stat = hdf5write_mat3d(vdbg_id, "skyr", &skyr);
    close_hdf5(vdbg_id);
    if (stat < 0)
        WARNING("Couldn't write skyr to %s", debugfn);
#endif

#ifdef COMMENT
    printf("Valid pixel count = %d\n", valpix);
    fflush(stdout);
#endif

    if (valpix < 1)
    {
        INFO("GRANULE CONTAINS NO VALID PIXELS", NONE);
        INFO("No valid RTTOV data was generated for this granule. It should be skipped.", NONE);
        exit(__LINE__);
    }

    // Release memory 
    mat_clear(&meshgridLat);
    mat_clear(&meshgridLon);
    mat_clear(&cropSP);
    mat_clear(&cropLat);
    mat_clear(&cropLon);
    mat_clear(&cropTCW);
    mat_clear(&cropskt);
    mat_clear(&cropt2);
    mat_clear(&cropq2);
    mat_clear(&cropHsurf);
    mat_clear(&cropSatZen);
    mat_clear(&cropQsurf);
    mat_clear(&cropTsurf);
    mat_clear(&cropST);
    mat3d_clear(&Pmesh);
    mat3d_clear(&cropT);
    mat3d_clear(&cropQ);
    mat3d_clear(&Bemis);

    // ========================================================
    // ECOSTRESS Cloud Data
    // ========================================================

#ifdef COMMENT
    printf("Processing L2_CLOUD\n");
    fflush(stdout);
#endif
#ifdef ENABLE_TIMING
    time_t startcloud = getusecs();
#endif
    // Load the LookUp Table (LUT) for Radiance to brightness conversion.

    // Open the lookup table file.
    FILE *flut = fopen(rad_lut_file, "r");
    if (!flut)
    {
        perror(rad_lut_file);
        ERREXIT(14, "Unable to open lookup table from %s", rad_lut_file);
    }

    // Count the number of lines in the LUT file.
    nlut_lines = 0;
    int ch;
    while (EOF != (ch=getc(flut)))
    {
        if (ch == '\n') nlut_lines++;
    }
    rewind(flut);

    // Interpolation requires minimum two lines of values.
    if (nlut_lines < 2)
        ERREXIT(15, "File %s is invalid. Must have a minimum of two lines of values.",
            rad_lut_file);

    // Create space for the LUT data.
    for (i = 0; i < LUT_COLS; i++)
    {
        lut[i] = (double*)malloc(nlut_lines * sizeof(double));
    }

    // Read the LUT data. The interpolation will generally go faster
    // if the file is read into memory in descending order. Note:
    // a binary search was implemented in the interp1d_npts function
    // so it no longer matters wheter the LUT lines are stored in
    // ascending or descending order.
    int fin;
    int numin;
    int scanok;
    for (fin = nlut_lines-1; fin >= 0; fin--)
    {
        for (numin = 0; numin < LUT_COLS; numin++)
        {
            scanok = fscanf(flut, "%lf", &lut[numin][fin]);
            if (!scanok)
                ERREXIT(15, "Error reading file %s at line %d, column %d",
                    rad_lut_file, fin, numin);
        }
    }
    fclose(flut);

    // Specify compression for the output product
    set_hdf5_2d_compression_on(
        n_lines / chunk_factor, n_pixels / chunk_factor, compression_deflate);

    // Process cloud
    process_cloud(&vRAD, cloud_lut_file, cloud_btdiff_file,
                  lut, nlut_lines, cloud_filename);

    // Read Cloud_test_final to get cloud mask.
    MatUint8 cloud;
    mat_uint8_init(&cloud);
    hid_t cloudid = open_hdf5(cloud_filename, H5F_ACC_RDONLY);
    if (cloudid < 0)
        ERREXIT(17, "Error opening %s for read.", cloud_filename);
    stat = hdf5read_mat_uint8(cloudid, "/SDS/Cloud_test_final", &cloud);
    if (stat == -1)
        ERREXIT(17, "Error reading /SDS/Cloud_test_final from %s", cloud_filename);
    close_hdf5(cloudid);

    // Compute statistical attribute values
    double slapse = 5.0;          // Surface temperature lapse rate (~5K/1000m)
    double BT_land_day = 280.0;
    double BT_land_night = 270.0;
    double BT_ocean = 270.0;      // Based on MODIS thresholds 
    slapse = getRuntimeParameter_f64("slapse", slapse);
    BT_land_day = getRuntimeParameter_f64("BT_land_day", BT_land_day);
    BT_land_night = getRuntimeParameter_f64("BT_land_night", BT_land_night);

    // Night or Day?
    double BT_land = BT_land_day;
    char *DayNightFlag = "Day"; // for file-level metadata
    double mean_solzen = 
        mean(vRAD.Solzen.vals, vRAD.Solzen.size1 * vRAD.Solzen.size2);
    if (mean_solzen > 85.0)
    {
        BT_land = BT_land_night;
        DayNightFlag = "Night";
    }

    Vector btvec;
    vec_init(&btvec);
    vec_alloc(&btvec, cloud.size1 * cloud.size2);
    size_t btcount = 0;
    double btmin = 9999.0; // initial values for finding minimum
    double btmax = 0.0;   // and maximumm values.
    double Hgt_km, Tba, BT11;
    for (i = 0; i < cloud.size1 * cloud.size2; i++)
    {
        if (cloud.vals[i] == 1)
        {
            Hgt_km = vRAD.El.vals[i];
            Tba = Hgt_km * slapse;
            if (vRAD.Water.vals[i] != 0)
            {
                BT11 = BT_ocean - Tba;
            }
            else
            {
                BT11 = BT_land - Tba;
            }
            btvec.vals[btcount] = BT11;
            btcount++;
            if (BT11 > btmax) btmax = BT11;
            if (BT11 < btmin) btmin = BT11;
        }
    }
    double btmean = mean(btvec.vals, btcount);
    double btstddev = stddev(btvec.vals, btcount);
    int32 cloud_percent = (int32)(0.5 + 100.0 *
            (double)btcount / 
            (double)(vRAD.Rad[0].size1 * vRAD.Rad[0].size2));
    updateValue(&cloud_metadata, "QAPercentCloudCover",
        DFNT_INT32, 1, dupint32(cloud_percent));
    updateValue(&cloud_metadata, "CloudMeanTemperature",
        DFNT_FLOAT64, 1, dupfloat64(btmean));
    updateValue(&cloud_metadata, "CloudMaxTemperature",
        DFNT_FLOAT64, 1, dupfloat64(btmax));
    updateValue(&cloud_metadata, "CloudMinTemperature",
        DFNT_FLOAT64, 1, dupfloat64(btmin));
    updateValue(&cloud_metadata, "CloudSDevTemperature",
        DFNT_FLOAT64, 1, dupfloat64(btstddev));
    // Add parameter from L1B_GEO product metadata
    updateValue(&cloud_metadata, "OrbitCorrectionPerformed",
        DFNT_CHAR8, strlen(vRAD.OrbitCorrectionPerformed), 
        strdup(vRAD.OrbitCorrectionPerformed));

    updateValue(&lste_metadata, "QAPercentCloudCover",
        DFNT_INT32, 1, dupint32(cloud_percent));
    updateValue(&lste_metadata, "CloudMeanTemperature",
        DFNT_FLOAT64, 1, dupfloat64(btmean));
    updateValue(&lste_metadata, "CloudMaxTemperature",
        DFNT_FLOAT64, 1, dupfloat64(btmax));
    updateValue(&lste_metadata, "CloudMinTemperature",
        DFNT_FLOAT64, 1, dupfloat64(btmin));
    updateValue(&lste_metadata, "CloudSDevTemperature",
        DFNT_FLOAT64, 1, dupfloat64(btstddev));
    // Add parameter from L1B_GEO product metadata to LSTE also.
    updateValue(&lste_metadata, "OrbitCorrectionPerformed",
        DFNT_CHAR8, strlen(vRAD.OrbitCorrectionPerformed), 
        strdup(vRAD.OrbitCorrectionPerformed));

    vec_clear(&btvec);

#ifdef ENABLE_TIMING
    printf("L2_CLOUD time: %.3f seconds\n",
            toseconds(elapsed(startcloud)));
    fflush(stdout);
#endif

    // ========================================================
    //  Run TES via STD or WVS
    // ========================================================
    Mat3d g, gi;
    mat3d_init(&g);
    mat3d_alloc(&g, n_channels, n_lines, n_pixels);
    mat3d_init(&gi);
    mat3d_alloc(&gi, n_channels, n_lines, n_pixels);
    // For the initial pass, fill all gi = 1.0.
    mat3d_fill(&gi, 1.0);
    // Initialize the output matrices.
    Matrix Ts;
    mat_init(&Ts);
    Mat3d emisf;
    mat3d_init(&emisf);
    MatUint16 QC;
    mat_uint16_init(&QC);

    // For first pass, all gp_water == 1
    MatUint8 gp_water;
    mat_uint8_init(&gp_water);
    mat_uint8_alloc(&gp_water, n_lines, n_pixels);
    mat_uint8_fill(&gp_water, 1);

    Mat3d Tg;
    mat3d_init(&Tg);

    //Mat3d skyi;
    //mat3d_init(&skyi);
    //mat3d_alloc(&skyi, n_channels, n_lines, n_pixels);
    //double skyi_val;

    Mat3d surfradi;
    mat3d_init(&surfradi);
    mat3d_alloc(&surfradi, n_channels, n_lines, n_pixels);
    double surfradi_val;
/*
    // Read the SKY TES parameter file
    FILE *sky_tes = fopen(sky_tes_file, "r");
    if (!sky_tes)
    {
        perror(sky_tes_file);
        ERREXIT(14, "Unable to open parameter file %s", sky_tes_file);
    }

    // Read the SKY data. 
    int isky, jsky;
    for (isky = 0; isky < n_channels; isky++)
    {
        for (jsky = 0; jsky < 3; jsky++)
        {
            scanok = fscanf(sky_tes, "%lf", &c_sky[isky][jsky]);
            if (!scanok)
                ERREXIT(15, "Error reading file %s at line %d, column %d",
                    sky_tes_file, isky+1, jsky+1);
        }
    }
    fclose(sky_tes);

#ifdef TESTING
    printf("%s data:\n", sky_tes_file);
    for (isky = 0; isky < n_channels; isky++)
    {
        printf("%d: ", isky+1);
        for (jsky = 0; jsky < 3; jsky++)
        {
            printf("%.4f ", c_sky[isky][jsky]);
        }
        printf("\n");
    }
    printf("\n");
#endif
*/
    // ========================================================
    //  TG WVS
    // ========================================================
    if (run_tgwvs)
    {
#ifdef COMMENT
        printf("Calculate Tg and WVS parameters\n");
        fflush(stdout);
#endif

#ifdef ENABLE_TIMING
        time_t start_tg_wvs = getusecs();
#endif

        // Get Tg surface brightness temps
        stat = tg_wvs(&Tg, &pwv, &vRAD.Satzen, &emis_aster, &Y, 
            &vRAD.Solzen, wvs_coef_file, lut, nlut_lines);
        if (stat < 0)
            ERREXIT(34, "tg_wvs failed to complete successfully.", NONE);

        // Interpolate brightness from Tg.
        Mat3d B;
        mat3d_init(&B);
        mat3d_alloc(&B, n_channels, n_lines, n_pixels);
        double Tg_val;
        for(b = 0; b < n_channels; b++) 
        {
            for(iline = 0; iline < n_lines; iline++) 
            {
                for(ipixel = 0; ipixel < n_pixels; ipixel++) 
                {
                    // to get nan in high sensor zenith angle pixels
                    mat3d_set(&g, b, iline, ipixel, real_nan); 

                    Tg_val = mat3d_get(&Tg, b, iline, ipixel);
                    if(is_nan(Tg_val))
                    {
                        mat3d_set(&B, b, iline, ipixel, real_nan);
                    }
                    else
                    {
                        // Interpolate the B value from the lut data.
                        mat3d_set(&B, b, iline, ipixel, 
                            interp1d_npts(lut[0], lut[band[b]+1], Tg_val, nlut_lines));
                    }
#ifdef DEBUG_PIXEL
                    if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
                    {
                        printf("(%d,%d) b=%d Tg=%.1f B=%.4f using LUT band %d\n",
                                iline, ipixel, b, Tg_val, 
                                mat3d_get(&B, b, iline, ipixel), band[b]+1);
                    }
#endif
                }
            }
        }

#ifdef TESTING
        describe_mat3d("Tg", &Tg);
        describe_mat3d("B", &B);
#endif

        // Compute gamma terms
        double gf; // exponent per band
        double numerator, denominator;
        Mat3d term_1, term_2t, term_3;
        mat3d_init(&term_1);
        mat3d_alloc(&term_1,  n_channels, n_lines, n_pixels);
        double term_1_val;
        mat3d_init(&term_2t);
        mat3d_alloc(&term_2t, n_channels, n_lines, n_pixels);
        double term_2t_val;
        mat3d_init(&term_3);
        mat3d_alloc(&term_3,  n_channels, n_lines, n_pixels);
        double term_3_val;

        for(b=0; b<n_channels; b++) 
        {
            // compute gf exponent for band b
            gf = pow(g2,bmp[band[b]]);
            for(iline=0; iline < n_lines; iline++) 
            {
                for(ipixel = 0; ipixel < n_pixels; ipixel++) 
                {
                    if (is_nan(mat3d_get(&Y, b, iline, ipixel)) ||
                        is_nan(mat3d_get(&B, b, iline, ipixel)))
                    {
                        mat3d_set(&term_1, b, iline, ipixel, real_nan);
                        mat3d_set(&term_2t, b, iline, ipixel, real_nan);
                        mat3d_set(&term_3, b, iline, ipixel, real_nan);
                    }
                    else
                    {
                        t1_val = mat3d_get(&t1r, b, iline, ipixel);
                        t2_val = mat3d_get(&t2r, b, iline, ipixel);
                        if (is_nan(t1_val))
                        {
#ifdef TESTING
                            WARNING("t1r[%d,%d,%d] == nan", b, iline, ipixel);
#endif
                            mat3d_set(&term_1, 0, iline, ipixel, real_nan);
                            continue;
                        }

                        if(t1_val < 0.0) 
                        {
#ifdef TESTING
                            WARNING(
                                "t1r[%d,%d,%d] < 0, which will make term_1 undefined\n",
                                b,iline,ipixel);
#endif
                            mat3d_set(&term_1, 0, iline, ipixel, real_nan);
                            continue;
                        }

                        // term_1 = t2r / t1r^gf
                        mat3d_set(&term_1, b, iline, ipixel,
                            t2_val / pow(t1_val, gf));
                        // term_2t = (B - pathr/(1 - t1r)) / (Y - pathr/(1 - t1r))
                        numerator = mat3d_get(&B, b, iline, ipixel) -
                            mat3d_get(&pathr, b, iline, ipixel) /
                            (1.0 - t1_val);
                        denominator = mat3d_get(&Y, b, iline, ipixel) -
                            mat3d_get(&pathr, b, iline, ipixel) /
                            (1.0 - t1_val);
                        mat3d_set(&term_2t, b, iline, ipixel, numerator/denominator);
                        // term_3 = t2r / t1r
                        mat3d_set(&term_3, b, iline, ipixel, t2_val / t1_val);
                    }
                }
            } 
        }

        //% Find where term2 is negative in any band (computing g will lead to nan
        //% if negative values not found)
        bool nan_or_neg;
        for(iline = 0; iline < n_lines; iline++)
        {
            for(ipixel = 0; ipixel < n_pixels; ipixel++)
            {
                nan_or_neg = false;
                for (b = 0; b < n_channels; b++)
                {
                    term_2t_val = mat3d_get(&term_2t, b, iline, ipixel);
                    if (is_nan(term_2t_val))
                    {
                        nan_or_neg = true;
                        break;
                    }
                    if (term_2t_val < 0.0)
                    {
                        nan_or_neg = true;
                        break;
                    }
                }
                // If a nan or negative value was found in any band, 
                // set the terms to nan in all bands.
                if (nan_or_neg)
                {
                    for(b=0; b<n_channels; b++) 
                    {
                        mat3d_set(&term_1, b, iline, ipixel, real_nan);
                        mat3d_set(&term_2t, b, iline, ipixel, real_nan);
                        mat3d_set(&term_3, b, iline, ipixel, real_nan);
                    }
                }
            }
        }

        double g_val;
        double term_2;
        // for each of n_channels bands
        for(b=0; b<n_channels; b++) 
        {
            gf = pow(g2,bmp[band[b]]);
            // for each pixel
            for(iline=0; iline<n_lines; iline++) 
            {
                for(ipixel = 0; ipixel < n_pixels; ipixel++) 
                {
                    term_1_val = mat3d_get(&term_1, b, iline, ipixel);
                    if(is_nan(term_1_val)) continue;
                    term_2t_val = mat3d_get(&term_2t, b, iline, ipixel);
                    if(is_nan(term_2t_val)) continue;
                    term_3_val = mat3d_get(&term_3, b, iline, ipixel);
                    if(is_nan(term_3_val)) continue;
                    term_2 = pow(term_2t_val,g1-gf);
                    //if((term_1_val * term_2) < 0.0)
                    //    continue;
                    g_val = 
                        log(term_1_val * term_2) / log(term_3_val);
                    // if either log() is complex, log returns nan.
                    // Replace nan gamma values with 1.0.
                    if (is_nan(g_val)) g_val = 1.0;
                    mat3d_set(&g, b, iline, ipixel, g_val);
#ifdef DEBUG_PIXEL
                    if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
                    {
                        printf("{%d}(%d, %d) term_1=%.4f term_2t=%.4f term_2=%.4f "
                                "term_3=%.4f B=%.4f g=%.4f\n",
                                b, iline, ipixel, 
                                mat3d_get(&term_1, b, iline, ipixel),
                                mat3d_get(&term_2t, b, iline, ipixel),
                                term_2,
                                mat3d_get(&term_3, b, iline, ipixel),
                                mat3d_get(&B, b, iline, ipixel),
                                g_val);
                    }
#endif
                }
            }
        }

        mat3d_clear(&term_1);
        mat3d_clear(&term_2t);
        mat3d_clear(&term_3);

        // Final steps to gamma (g) modification
        uint8 gpin;
        //int count_gpin = 0;
        //double PWV_total = 0.0;
        //double PWV_count = 0;
        double TG_total = 0.0;
        double TG_count = 0;
        bool is_cloud;
        double PWV_val;
        for(iline=0; iline < n_lines; iline++) 
        {
            for(ipixel = 0; ipixel < n_pixels; ipixel++) 
            {
                /* 2018 removed PWVmean
                // PWV Low Test
                PWV_val = mat_get(&pwv, iline, ipixel);
                if (!is_nan(PWV_val))
                {
                    PWV_total += PWV_val;
                    PWV_count++;
                }
                */

                // Get mean Tg over land.
                if (mat_uint8_get(&vRAD.Water, iline, ipixel) == 0)
                {
                    Tg_val = mat3d_get(&Tg, band_index[BAND_4], iline, ipixel);
                    if (!is_nan(Tg_val))
                    {
                        TG_total += Tg_val;
                        TG_count++;
                    }
                }

                //%% Get Graybody pixels (veg, ocean, ice, snow)
                gpin = mat_uint8_get(&gp_water, iline, ipixel);

                //% Over bare pixels (gpin == 1) use LST band g value (g{BAND_5}) 
                //  for all bands.
                if (gpin == 1)
                {
                    // ECOSTRESS use BAND_5 gpin pixels.
                    g_val = mat3d_get(&g, band_index[BAND_5], iline, ipixel);
                    for (b = 0; b < n_channels; b++)
                    {
                        //count_gpin++;
                        mat3d_set(&g, b, iline, ipixel, g_val);
                    }
                }

                // Prepare gi for smoothing.
                is_cloud = (mat_uint8_get(&cloud, iline, ipixel) & 3) == 3;
                for (b = 0; b < n_channels; b++)
                {
                    if (is_cloud)
                        g_val = real_nan;
                    else
                    {
                        g_val = mat3d_get(&g, b, iline, ipixel);
                        // Clamp to -2.0 .. 3.0
                        if (g_val > 3.0)
                            g_val = 3.0;
                        else if (g_val < -2.0)
                            g_val = -2.0;
                    }
                    mat3d_set(&gi, b, iline, ipixel, g_val);
                }
            }
        }

        // Threshold check
        //double PWVmean = PWV_total / (double) PWV_count;
        double TGmean = TG_total / (double) TG_count;

        if (TGmean < TGthresh)
        {
            smooth_scale = smooth_scale1;
            PWVthresh = PWVthresh1;
        }
        else
        {
            smooth_scale = smooth_scale2;
            PWVthresh = PWVthresh2;
        }
#ifdef COMMENT
        printf("TGmean=%.1f TGthresh=%.1f PWVthresh=%.1f smooth_scale=%d\n",
                TGmean, TGthresh, PWVthresh, smooth_scale);
#endif

        // Use PWV to set gi=1.0 when PWV is low.
        for (iline = 0; iline < n_lines; iline++)
        {
            for (ipixel = 0; ipixel < n_pixels; ipixel++)
            {
                PWV_val = mat_get(&pwv, iline, ipixel);
                if (PWV_val < PWVthresh)
                {
                    // pwvlow
                    for (b = 0; b < n_channels; b++)
                    {
                        mat3d_set(&gi, b, iline, ipixel, 1.0);
                    }
                }
            }
        }

#ifdef TESTING
        time_t startsmooth2d = getusecs();
#endif

        // Perform smoothing per level (channel) on gi.
        for (b = 0; b < n_channels; b++)
        {
            smooth2d(&gi.vals[b * n_pixels * n_lines], 
                n_lines, n_pixels, smooth_scale, smooth_scale);
        }

#ifdef TESTING
        printf("smooth2d time: %.3f seconds\n",
                toseconds(elapsed(startsmooth2d)));
#endif

        // After smoothing, set gi for cloud pixels to 1.0
        for (iline = 0; iline < n_lines; iline++)
        {
            for (ipixel = 0; ipixel < n_pixels; ipixel++)
            {
                is_cloud = (mat_uint8_get(&cloud, iline, ipixel) & 3) == 3;
                if (is_cloud)
                {
                    for (b = 0; b < n_channels; b++)
                    {
                        mat3d_set(&gi, b, iline, ipixel, 1.0);
                    }
                }
            }
        }

        // Apply WVS to atmospheric parameters
        //const double to_rads = M_PI / 180.0;
        //double path_z;
        double gi_val;
        double ti;
        double pathi;

        for (b = 0; b < n_channels; b++)
        {
            gf = pow(g2, bmp[band[b]]);
            for (iline = 0; iline < n_lines; iline++)
            {
                for (ipixel = 0; ipixel < n_pixels; ipixel++)
                {
                    // ti{b} = t1{b}.^((gi{b} - gf)./(1 - gf)).*t2{b}.^((1 - gi{b})./(1 - gf));
                    gi_val = mat3d_get(&gi, b, iline, ipixel);
                    t1_val = mat3d_get(&t1r, b, iline, ipixel);
                    t2_val = mat3d_get(&t2r, b, iline, ipixel);
                    ti = pow(t1_val, (gi_val - gf) / (1.0 - gf))
                           * pow(t2_val, ((1.0 - gi_val) / (1.0 - gf)));
                    // pathi{b} = path{b}.*((1 - ti{b})./(1 - t1{b}));
                    path_val = mat3d_get(&pathr, b, iline, ipixel);
                    pathi = path_val * (1.0 - ti) / (1.0 - t1_val);
/*
                    // path_z = pathi{b}.*((1 - ti{b}.^cos(SatZen.*pi/180))./(1 - ti{b}));
                    satzen_rad = mat_get(&vRAD.Satzen, iline, ipixel) * to_rads;
                    path_z = pathi * 
                        (1.0 - pow(ti, cos(satzen_rad)) / (1.0 - ti));
                    // skyi_c = sc{b}(1) + sc{b}(2).*path_z + sc{b}(3).*path_z.^2;
                    skyi_val = c_sky[b][0] +
                               c_sky[b][1] * path_z + 
                               c_sky[b][2] * path_z * path_z;
*/
                    // surfradi_c = (Rad{b} - pathi{b})./ti{b};
                    surfradi_val = 
                        (mat3d_get(&Y, b, iline, ipixel) - pathi) / ti;
                    if (surfradi_val < 0.0)
                    {
                        surfradi_val = real_nan;
                        //skyi_val = real_nan;
                    }
                    //mat3d_set(&skyi, b, iline, ipixel, skyi_val);
                    mat3d_set(&surfradi, b, iline, ipixel, surfradi_val);
#ifdef DEBUG_PIXEL
                    if (ipixel == DEBUG_PIXEL && iline == DEBUG_LINE)
                    {
                        printf("DEBUG {%d}(%d,%d) gi=%.4f t1r=%.4f t2r=%.4f skyr=%.4f "
                               "pathr=%.4f pathi=%.4f Rad=%.4f ti=%.4f surfradi=%.4f\n",
                            b, iline, ipixel, gi_val, t1_val, t2_val,
                            mat3d_get(&skyr, b, iline, ipixel), path_val, pathi, 
                            mat3d_get(&Y, b, iline, ipixel), ti, surfradi_val);
                    }
#endif
               }
            }
        }
#ifdef ENABLE_TIMING
        printf("TG_WVS time: %.3f seconds\n",
                toseconds(elapsed(start_tg_wvs)));
        fflush(stdout);
#endif
#ifdef TESTING
      char tg_wvs_debug_name[PATH_MAX];
      sprintf(tg_wvs_debug_name, "%s/tg_wvs_debug.h5", product_path);
      hid_t tgwvsdbg = create_hdf5(tg_wvs_debug_name);
      if (tgwvsdbg == -1)
          WARNING("Could not create output file %s", tg_wvs_debug_name);
      stat = hdf5write_mat3d(tgwvsdbg, "Tg", &Tg);
      stat += hdf5write_mat3d(tgwvsdbg, "g", &g);
      stat += hdf5write_mat3d(tgwvsdbg, "gi", &gi);
      close_hdf5(tgwvsdbg);
      if (stat < 0)
          WARNING("%d error(s) writing tg_wvs_debug.h5", -stat);
#endif

    }
    // ========================================================
    //  Run STD only -- no WVS
    // ========================================================
    else
    {
#ifdef COMMENT
        printf("Calculate Surface Radiances\n");
        fflush(stdout);
#endif

#ifdef ENABLE_TIMING
        time_t start_surfrad = getusecs();
#endif

        for (b = 0; b < n_channels; b++)
        {
            for (iline = 0; iline < n_lines; iline++)
            {
                for (ipixel = 0; ipixel < n_pixels; ipixel++)
                {
                    // surfradi{b} = (Rad{b} - path{b})./t1{b};
                    mat3d_set(&surfradi, b, iline, ipixel,
                        (mat3d_get(&Y, b, iline, ipixel) -
                         mat3d_get(&pathr, b, iline, ipixel))
                        / mat3d_get(&t1r, b, iline, ipixel));
#ifdef DEBUG_PIXEL
                    if (ipixel == DEBUG_PIXEL && iline == DEBUG_LINE)
                    {
                        printf("DEBUG {%d}(%d,%d) t1r=%.4f "
                               "pathr=%.4f Rad=%.4f surfradi=%.4f\n",
                            b, iline, ipixel, 
                            mat3d_get(&t1r, b, iline, ipixel), 
                            mat3d_get(&pathr, b, iline, ipixel), 
                            mat3d_get(&Y, b, iline, ipixel),
                            mat3d_get(&surfradi, b, iline, ipixel));
                    }
#endif
                }
            }
        }

#ifdef ENABLE_TIMING
        printf("Calc Surfrad time: %.3f seconds\n",
                toseconds(elapsed(start_surfrad)));
        fflush(stdout);
#endif
    }
    
    // ========================================================
    //  TES
    // ========================================================
#ifdef COMMENT
    printf("Run TES\n");
    fflush(stdout);
#endif

#ifdef ENABLE_TIMING
    time_t start_tes = getusecs();
#endif

    // Apply the TES algorithm, using the default gi=1.0.
    apply_tes_algorithm(&Ts, &emisf, &QC, &surfradi, &skyr, &gp_water, &t1r, &Y);

#ifdef ENABLE_TIMING
    printf("TES time: %.3f seconds\n",
            toseconds(elapsed(start_tes)));
    fflush(stdout);
#endif

#ifdef TESTING
    strcpy(debugfn, product_path);
    strcat(debugfn, "/ECOSTRESSdebug_tes.h5");
    vdbg_id = create_hdf5(debugfn);
    if (vdbg_id == FAIL)
        ERROR("Could not create %s", debugfn);
    describe_mat3d("surfradi", &surfradi);
    describe_mat3d("skyr", &skyr);
    describe_matrix("Ts", &Ts);
    describe_mat3d("emisf", &emisf);
    describe_mat_u16("QC", &QC);
    stat  = hdf5write_mat(vdbg_id, "Ts", &Ts);
    stat += hdf5write_mat3d(vdbg_id, "surfradi", &surfradi);
    stat += hdf5write_mat3d(vdbg_id, "skyr", &skyr);
    stat += hdf5write_mat3d(vdbg_id, "emisf", &emisf);
    stat += hdf5write_mat_uint16(vdbg_id, "QC", &QC);
    if (stat < 0)
        WARNING("%d errors writing debug data to %s.", -stat, debugfn);
    close_hdf5(vdbg_id);
#endif

    // ========================================================
    // Prepare the EmisWB data
    // ========================================================
#ifdef COMMENT
    printf("Prepare EmisWB\n");
    fflush(stdout);
#endif

#ifdef ENABLE_TIMING
    time_t start_emiswb = getusecs();
#endif

    Matrix emis_wb;
    mat_init(&emis_wb);
    mat_alloc(&emis_wb, n_lines, n_pixels);
    double emis_val;
    double sum_emis_wb;
#ifdef TESTING
    int nonzero[MAX_CHANNELS] = {0,0,0,0,0};
#endif
    for (iline = 0; iline < n_lines; iline++)
    {
        for(ipixel=0; ipixel < n_pixels; ipixel++)
        {
            // Initialize EmisWB sum for Wide Band Emissivity.
            sum_emis_wb = emis_wb_coeffs[n_channels];
            // Use emissivity from each channel.
            for(b=0; b<n_channels; b++)
            {
                emis_val = mat3d_get(&emisf, b, iline, ipixel);
#ifdef TESTING
                if (emis_val != 0.0) ++nonzero[b];
#endif
#ifdef DEBUG_PIXEL
                if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
                {
                    printf("DEBUG emis%d(%d, %d)=%.4f EmisWbCoef=%.4f\n", 
                            band[b] + 1, iline, ipixel, emis_val, emis_wb_coeffs[b]);
                }
#endif
                // Accumulate the sum of coefficient*emisf[b] for emis_wb.
                sum_emis_wb += emis_val * emis_wb_coeffs[b];
            }
            // Save the computed sum as the emis_wb value.
            mat_set(&emis_wb, iline, ipixel, sum_emis_wb);
#ifdef DEBUG_PIXEL
            if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
            {
                printf("DEBUG emis_wb(%d, %d)=%.4f\n", iline, ipixel, sum_emis_wb);
            }
#endif
        }
    }

#ifdef ENABLE_TIMING
    printf("EmisWB time: %.3f seconds\n",
            toseconds(elapsed(start_emiswb)));
    fflush(stdout);
#endif

    // ========================================================
    // Error Estimation Processing
    // ========================================================
#ifdef COMMENT
    printf("Performing error estimation.\n");
    fflush(stdout);
#endif
#ifdef ENABLE_TIMING
    time_t starterrest = getusecs();
#endif

    // Create de and dT for error estimates
    Mat3d de;
    mat3d_init(&de);
    mat3d_alloc(&de, n_channels, emisf.size2, emisf.size3);
    Matrix dT;
    mat_init(&dT);
    mat_alloc(&dT, emisf.size2, emisf.size3);
    
    // Variables for error estimate loop
    double TCW; // Total Column Water from pwv
    double SVA; // Sensor Viewing Angle from vRAD.Satzen
    double det;
    double dTt;
    double se; // sum squared emis error

    // Get PWVavg for use in error estimation
    //double PWVavg = mean(pwv.vals, pwv.size1 * pwv.size2);
#ifdef COMMENT
    //if (PWVavg <= 2.0)
    //    printf("Mean PWV of %.4f is <= 2.0 threshold.\n", PWVavg);
    //fflush(stdout);
#endif

    uint16 qc_val;
    for (iline = 0; iline < n_lines; iline++)
    {
        for (ipixel = 0; ipixel < n_pixels; ipixel++)
        {
            qc_val = mat_uint16_get(&QC, iline, ipixel);

            // If QC already set to "not processed", skip further settings.
            if ((qc_val & 0x03) == 0x03)
                continue;

            // JAN 2018 -- cloud is anything with bits 0 & 1 set.
            if ((mat_uint8_get(&cloud, iline, ipixel) & 3) == 3)
            {
                // Input cloudy pixel
                // 2018-02-28 Bits 5&4 not set by L2 for Cloud/Ocean flag
                //qc_val |= (bits[4]+bits[5]);
                qc_val &= 0xFFFC;
                qc_val |= bits[1];
                mat_uint16_set(&QC, iline, ipixel, qc_val);
                //continue; // ECOSTRESS: produce error estimates over cloudy pixels.
            }

/* JAN 2018 -- not using extended cloud
            else if (mat_uint8_get(&cextend, iline, ipixel) == 4)
            {
                // Extended cloudy pixel
                //qc_val |= bits[5];
                qc_val &= 0xFFFC;
                qc_val |= bits[0];
            }
*/
            TCW = mat_get(&pwv, iline, ipixel);
            if (is_nan(TCW))
            {
                qc_val |= (bits[12]+bits[13]+bits[14]+bits[15]);
                mat_uint16_set(&QC, iline, ipixel, qc_val);
                continue;
            }

            SVA = mat_get(&vRAD.Satzen, iline, ipixel);

            // Estimate Emis error for each band.
            se = 0.0;
            for (b = 0; b < n_channels; b++)
            {
                //det = emiserr_coef0[b] + emiserr_coef1[b] * TCW;
                det = xe[band[b]][0] + xe[band[b]][1] * TCW + xe[band[b]][2] * TCW * TCW;
                mat3d_set(&de, b, iline, ipixel, det);
                se += det*det;
            }

            // Estimate LST error.
/*
            if (PWVavg <= 2.0)
            {
                dTt = lsterr_coef0 
                    + lsterr_coef1 * TCW 
                    + lsterr_coef2 * SVA;
            }
            else
            {
*/
/*
                dTt = xtt[0]+xtt[1]*TCW+xtt[2]*SVA+
                    xtt[3]*TCW*SVA+xtt[4]*TCW*TCW+
                    xtt[5]*SVA*SVA;
*/
                dTt = xt[0] + xt[1] * TCW + xt[2] * SVA;
/*
            }
*/
            mat_set(&dT, iline, ipixel, dTt);

            // Set quality flag bits to reflect emissivity error level.
            // Using RMSE of all emissivity bands.
            det = sqrt(se / n_channels);
            if(det >= emiserr_lims[2]) {    // leave QC bits 0
            } else if(det > emiserr_lims[1]) {
                qc_val |= bits[12];
            } else if(det > emiserr_lims[0]) {
                qc_val |= bits[13];
            } else 
                qc_val |= (bits[12]+bits[13]);

            // Set quality flag bits to reflect land surface temperature 
            // error level.
            if(dTt >= Terr_lims[2]) {    // leave QC bits 0 
            } else if(dTt > Terr_lims[1]) {
                qc_val |= bits[14];
            } else if(dTt > Terr_lims[0]) {
                qc_val |= bits[15];
            } else {
                qc_val |= (bits[14]+bits[15]);
            }
            mat_uint16_set(&QC, iline, ipixel, qc_val);
        }
    }

    // Missing Scan Processing
    //
    // Account for additional uncertainty due to missing scan lines that 
    // were filled with neural net.
    // Missing scan: Rad_quality in bands 1/5 = 1
    bool MissScanFilled;
    bool MissScanNotFilled;
    bool MissingBad;
    bool NotSeen;
    const double LST_ScanERR = 0.34;
    const double Emis_ScanErr[MAX_CHANNELS] =  
        {0.0099, 0.0053, 0.0050, 0.0047, 0.0023};

    int32 dq_val;
    for (iline = 0; iline < n_lines; iline++)
    {
        for (ipixel = 0; ipixel < n_pixels; ipixel++)
        {
            MissScanFilled = false;
            MissScanNotFilled = false;
            MissingBad = false;
            NotSeen = false;
            for (b = 0; b < n_channels; b++)
            {
                dq_val = mat_int32_get(&vRAD.DataQ[b], iline, ipixel);
                if (dq_val == 1) MissScanFilled = true;
                if (dq_val == 2) MissScanNotFilled = true;
                if (dq_val == 3) MissingBad = true;
                if (dq_val == 4) NotSeen = true;
            }
            if (MissScanNotFilled || MissingBad || NotSeen)
            {
                // Missing/Bad or NotSeen reflected as "not produced"
                mat_uint16_set(&QC, iline, ipixel, 
                    mat_uint16_get(&QC, iline, ipixel) | 
                            bits[3]|bits[2]|bits[1]|bits[0]);
            }
            else if (MissScanFilled)
            {
                // Missing Scan reflected in QC
                qc_val = mat_uint16_get(&QC, iline, ipixel);
                // Indicate missing scan in bits 3&2 = 01
                qc_val |= bits[2];
                // If QC 1&0 if cloud or not-produced, do not override. 
                // Else, set 1&0 = 01
                if ((qc_val & 3) == 0)
                    qc_val |= bits[0];
                mat_uint16_set(&QC, iline, ipixel, qc_val);
                // Increase error estimates for missing scan pixels.
                mat_set(&dT, iline, ipixel,
                    mat_get(&dT, iline, ipixel) + LST_ScanERR);
                for (b = 0; b < n_channels; b++)
                {
                    mat3d_set(&de, b, iline, ipixel,
                        mat3d_get(&de, b, iline, ipixel) + 
                        Emis_ScanErr[band[b]]);
                }
            }
        }
    }

    // Determine the "good" totals for metadata
    int count_good = 0;
    double tot_good_lst = 0.0;
    double tot_good_emis[n_channels];
    for (b = 0; b < n_channels; b++) tot_good_emis[b] = 0.0;
    for (iline = 0; iline < n_lines; iline++) 
    {
        for(ipixel=0; ipixel < n_pixels; ipixel++) 
        {
            // "good" is definied by bits 0,1 of QC being clear.
            qc_val = mat_uint16_get(&QC, iline, ipixel);
            if ((qc_val & 3) == 0)
            {
                // Mandatory QA is "Good"
                double Ts_val = mat_get(&Ts, iline, ipixel);
                if (is_nan(Ts_val))
                {
                    // QC should not be "Good" if Ts is NaN
                    qc_val |= 0x000F; // set "not produced, missing/bad L1B data"
                    mat_uint16_set(&QC, iline, ipixel, qc_val);
                }
                else if (Ts_val > 380.0)
                {
                    // This is too high a value to be "good"
                    // Set Data Quality = 11 (missing/bad L1B) and 
                    // Mandatory QA = 01 (nominal).
                    qc_val &= 0xFFF0;
                    qc_val |= 0x000D;
                    mat_uint16_set(&QC, iline, ipixel, qc_val);
                }
                else
                {
                    count_good++;
                    tot_good_lst += mat_get(&Ts, iline, ipixel);
                    for (b = 0; b < n_channels; b++)
                    {
                        tot_good_emis[b] += 
                            mat3d_get(&emisf, b, iline, ipixel);
                    }
                }
            }
            // Be sure that "not produced" is set to fill value
            else if ((qc_val & 3) == 3) // not produced
            {
                mat_set(&Ts, iline, ipixel, REAL_NAN);
                for (b = 0; b < n_channels; b++)
                {
                    mat3d_set(&emisf, b, iline, ipixel, REAL_NAN);
                }
            }
        }
    }

#ifdef ENABLE_TIMING
    printf("Error estimate time: %.3f seconds\n",
            toseconds(elapsed(starterrest)));
    fflush(stdout);
#endif

#ifdef TESTING
    describe_mat_u16("QC", &QC);
    describe_mat3d("de", &de);
    describe_matrix("dT", &dT);
    strcpy(debugfn, product_path);
    strcat(debugfn, "/ECOSTRESSdebug3.h5");
    vdbg_id = create_hdf5(debugfn);
    if (vdbg_id == FAIL)
        ERROR("Could not create %s", debugfn);
    stat = hdf5write_mat_uint16(vdbg_id, "QC", &QC);
    stat += hdf5write_mat3d(vdbg_id, "de", &de);
    stat += hdf5write_mat(vdbg_id, "dT", &dT);
    if (stat < 0)
        WARNING("%d errors writing error estimate results to %s", 
                -stat, debugfn);
    close_hdf5(vdbg_id);
#endif

    // ========================================================
    //  Create the Output Data Product
    // ========================================================
#ifdef COMMENT
    printf("Writing output datasets to %s.\n", output_filename);
    fflush(stdout);
#endif

#ifdef ENABLE_TIMING
    time_t startoutput = getusecs();
#endif

    // Output datasets
    
    // Land Surface Temperature
    MatUint16 LST;
    mat_uint16_init(&LST);
    mat_uint16_alloc(&LST, emisf.size2, emisf.size3);

    // Emissivity bands
    MatUint8 Emis[MAX_CHANNELS];
    // Error estimates for emisf bands and LST.
    MatUint16 Emis_err[MAX_CHANNELS];
    for (b = 0; b < n_channels; b++)
    {
        mat_uint8_init(&Emis[b]);
        mat_uint8_alloc(&Emis[b], emisf.size2, emisf.size3);
        mat_uint16_init(&Emis_err[b]);
        mat_uint16_alloc(&Emis_err[b], emisf.size2, emisf.size3);
    }
    MatUint8 EmisWB;
    mat_uint8_init(&EmisWB);
    mat_uint8_alloc(&EmisWB, emisf.size2, emisf.size3);

    MatUint8 LST_err;
    mat_uint8_init(&LST_err);
    mat_uint8_alloc(&LST_err, emisf.size2, emisf.size3);

    // Other scaled datasets
    MatUint16 PWV;
    mat_uint16_init(&PWV);
    mat_uint16_alloc(&PWV, pwv.size1, pwv.size2);
    //MatInt16 NDVI;
    //mat_int16_init(&NDVI);
    //mat_int16_alloc(&NDVI, ndvi.size1, ndvi.size2);
/* 2018-02-28 Removed oceanpix and View_angle which can be obtained from L1B_GEO
    MatUint8 View_angle;
    mat_uint8_init(&View_angle);
    mat_uint8_alloc(&View_angle, vRAD.Satzen.size1, vRAD.Satzen.size2);
*/
    // Double values need to be converted by scale and offset into int/uint outputs.
    const double Scale_LST = 0.02;
    const double Scale_EMIS = 0.002;
    const double Offset_EMIS = 0.49;
    const double Scale_LST_err = 0.04;
    const double Scale_EMIS_err = 1.0e-4;
    const double Scale_PWV = 0.001;
    //const double Scale_NDVI = 1.0e-4;
    //const double Scale_VA = 0.5;
    const float32 Fill_LAT_LON = -9999.0;
    const uint16 Fill_Value_16 = 0;
    const uint8 Fill_Value_8 = 0;
    double d;

    // Convert internal f64 to external scaled ints.
    for (iline = 0; iline < n_lines; iline++)
    {
        for (ipixel = 0; ipixel < n_pixels; ipixel++)
        {
            d = mat_get(&Ts, iline, ipixel);
            if (is_nan(d) || d < 0.0) 
                mat_uint16_set(&LST, iline, ipixel, Fill_Value_16);
            else
                mat_uint16_set(&LST, iline, ipixel, (uint16)(d / Scale_LST + 0.5));
            
            d = mat_get(&dT, iline, ipixel);
            if (is_nan(d)) 
                mat_uint8_set(&LST_err, iline, ipixel, Fill_Value_8);
            else
                mat_uint8_set(&LST_err, iline, ipixel, (uint8)(d / Scale_LST_err + 0.5));

            d = mat_get(&pwv, iline, ipixel);
            if (is_nan(d)) 
                mat_uint16_set(&PWV, iline, ipixel, Fill_Value_16);
            else
                mat_uint16_set(&PWV, iline, ipixel, (uint16)(d / Scale_PWV + 0.5));

/* 2018-02-28 Removed oceanpix and View_angle which can be obtained from L1B_GEO
            d = mat_get(&vRAD.Satzen, iline, ipixel);
            if (is_nan(d)) 
                mat_uint8_set(&View_angle, iline, ipixel, Fill_Value_8);
            else
                mat_uint8_set(&View_angle, iline, ipixel, 
                    (uint8)(d / Scale_VA + 0.5));
*/

            d = mat_get(&vRAD.Lat, iline, ipixel);
            if (is_nan(d) || d < -90.0 || d > 90.0)
                mat_set(&vRAD.Lat, iline, ipixel, Fill_LAT_LON);

            d = mat_get(&vRAD.Lon, iline, ipixel);
            if (is_nan(d) || d < -180.0 || d > 180.0)
                mat_set(&vRAD.Lon, iline, ipixel, Fill_LAT_LON);

            for (b = 0; b < n_channels; b++)
            {
                d = mat3d_get(&emisf, b, iline, ipixel) - Offset_EMIS;
                if (is_nan(d) || d < 0.0) 
                    mat_uint8_set(&Emis[b], iline, ipixel, Fill_Value_8);
                else
                    mat_uint8_set(&Emis[b], iline, ipixel, 
                        (uint8)(d / Scale_EMIS + 0.5));
#ifdef DEBUG_PIXEL
                if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
                {
                    printf("DEBUG emisf(%d, %d, %d)=%.4f Emis=%d result=%.4f\n",
                            b, iline, ipixel, d + Offset_EMIS, 
                            mat_uint8_get(&Emis[b], iline, ipixel),
                            (double)mat_uint8_get(&Emis[b], iline, ipixel) * 0.002 + 0.49);
                }
#endif
                d = mat3d_get(&de, b, iline, ipixel);
                if (is_nan(d)) 
                    mat_uint16_set(&Emis_err[b], iline, ipixel, Fill_Value_16);
                else
                    mat_uint16_set(&Emis_err[b], iline, ipixel, 
                        (uint16)(d / Scale_EMIS_err + 0.5));
            }
            
            d = mat_get(&emis_wb, iline, ipixel) - Offset_EMIS;
            if (isnan(d) || d < 0.0)
                mat_uint8_set(&EmisWB, iline, ipixel, Fill_Value_8);
            else
                mat_uint8_set(&EmisWB, iline, ipixel,
                    (uint8)(d / Scale_EMIS + 0.5));
#ifdef DEBUG_PIXEL
            if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
            {
                printf("DEBUG emis_wb(%d, %d)=%.4f\n",
                        iline, ipixel, d + Offset_EMIS);
            }
#endif
        }
    }

    hid_t h5out = create_hdf5(output_filename);
    if (h5out == -1)
        ERREXIT(18, "Could not create output file %s", output_filename);
    //hid_t top_group = H5Gcreate2(h5out, "/LSTE", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t top_group = H5Gcreate2(h5out, "/SDS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (top_group == -1)
        ERREXIT(19, "Could not create group /SDS in %s", output_filename);
    H5Gclose(top_group);
/*
#ifdef PRODUCE_GEOLOCATION
    hid_t group1 = H5Gcreate2(h5out, "/LSTE/GeolocationFields", 
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group1 == -1)
        ERREXIT(19, "Could not create group /LSTE/GeolocationFields in %s", output_filename);
    H5Gclose(group1);
#endif
    hid_t group2 = H5Gcreate2(h5out, "/LSTE/DataFields", 
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group2 == -1)
        ERREXIT(19, "Could not create group /LSTE/DataFields in %s", output_filename);
    H5Gclose(group2);
*/
    stat = 0;
#ifdef PRODUCE_GEOLOCATION
    stat += hdf5write_mat_asfloat32(h5out, 
        "/LSTE/GeolocationFields/Latitude", &vRAD.Lat);
    stat += hdf5write_mat_asfloat32(h5out, 
        "/LSTE/GeolocationFields/Longitude", &vRAD.Lon);
#endif
    stat += hdf5write_mat_uint16(h5out, "/SDS/QC", &QC);
    stat += hdf5write_mat_uint16(h5out, "/SDS/LST", &LST);
    stat += hdf5write_mat_uint8(h5out, "/SDS/LST_err", &LST_err);
    char emis_ds[32];
    for (b = 0; b < n_channels; b++)
    {
        sprintf(emis_ds, "/SDS/Emis%d", band[b]+1);
        stat += hdf5write_mat_uint8(h5out, emis_ds, &Emis[b]);
        strcat(emis_ds, "_err");
        stat += hdf5write_mat_uint16(h5out, emis_ds, &Emis_err[b]);
    }
    if (n_channels == 3)
    {
        // Add dummy Emis1, Emis1_err, Emis3, Emis_err.
        MatUint8 dummy8;
        mat_uint8_init(&dummy8);
        mat_uint8_alloc(&dummy8, emisf.size2, emisf.size3);
        stat += hdf5write_mat_uint8(h5out, "/SDS/Emis1", &dummy8);
        stat += hdf5write_mat_uint8(h5out, "/SDS/Emis3", &dummy8);
        mat_uint8_clear(&dummy8);
        MatUint16 dummy16;
        mat_uint16_init(&dummy16);
        mat_uint16_alloc(&dummy16, emisf.size2, emisf.size3);
        stat += hdf5write_mat_uint16(h5out, "/SDS/Emis1_err", &dummy16);
        stat += hdf5write_mat_uint16(h5out, "/SDS/Emis3_err", &dummy16);
        mat_uint16_clear(&dummy16);
    }
    stat += hdf5write_mat_uint8(h5out, "/SDS/EmisWB", &EmisWB);
    stat += hdf5write_mat_uint16(h5out, "/SDS/PWV", &PWV);
    //stat += hdf5write_mat_int16(h5out, "/SDS/NDVI", &NDVI);
/* 2018-02-28 Removed oceanpix and View_angle which can be obtained from L1B_GEO
    stat += hdf5write_mat_uint8(h5out, 
        "/SDS/oceanpix", &vRAD.Water);
    stat += hdf5write_mat_uint8(h5out, 
        "/SDS/View_angle", &View_angle);
*/
    // Output cloud_mask.
    stat += hdf5write_mat_uint8(h5out, "/SDS/cloud_mask", &cloud);

    // Output water_mask (formerly oceanix).
    stat += hdf5write_mat_uint8(h5out, "/SDS/water_mask", &vRAD.Water);

    if (stat < 0)
    {
        ERREXIT(19, "%d error(s) occured while writing 20 datasets to %s",
            -stat, output_filename);
    }

#ifdef PRODUCE_GEOLOCATION
    
    // Add GeolocationFields to Cloud output
    cloudid = open_hdf5(cloud_filename, H5F_ACC_RDWR);
    stat = hdf5write_mat_asfloat32(cloudid, 
        "/SDS/GeolocationFields/Latitude", &vRAD.Lat);
    stat += hdf5write_mat_asfloat32(cloudid, 
        "/SDS/GeolocationFields/Longitude", &vRAD.Lon);
    if (stat < 0)
    {
        ERREXIT(17, "%d error(s) occured while writing Geolocation data to %s",
            -stat, cloud_filename);
    }
    close_hdf5(cloudid);
#endif

    // ========================================================
    //  Write output file metadata
    // ========================================================
#ifdef COMMENT
    printf("Writing metadata.\n");
    fflush(stdout);
#endif

    char *newval;
    int newlen;
    MetadataSet lste_set;
    constructMetadataSet(&lste_set, 0);
    MetadataSet cloud_set;
    constructMetadataSet(&cloud_set, 0);
    int e;

    // Copy the metadata from the GEO file
    hid_t meta_src_id = open_hdf5(GEO_filename, H5F_ACC_RDONLY);
    if (meta_src_id == FAIL)
        ERREXIT(25, "Unable to open file %s", GEO_filename);
    stat = loadMetadataGroup(&lste_set, meta_src_id, "StandardMetadata");
    if (stat < 0)
        ERREXIT(26, "Failed to read StandardMetadata from %s", GEO_filename);
    stat = loadMetadataGroup(&cloud_set, meta_src_id, "StandardMetadata");
    if (stat < 0)
        ERREXIT(26, "Failed to read StandardMetadata from %s", GEO_filename);
    close_hdf5(meta_src_id);

    // Add or replace elements with items from the Metadata group.
    char mnamebuf[CONFIG_VALUELEN];
    ConfigGroup *metadataGroup = getGroup(run_params, "Metadata");
    if (metadataGroup == NULL)
        ERREXIT(8, "Missing <group> Metadata in %s", run_params_filename);
    ConfigElement *metadataElem = metadataGroup->head;
    while (metadataElem != NULL)
    {
        for (e = 0; e < metadataElem->count; e++)
        {
            // If multiple values, make unique names for each.
            if (e == 0)
                strcpy(mnamebuf, metadataElem->name);
            else
                sprintf(mnamebuf, "%s%d", metadataElem->name, e);
            newlen = strlen(metadataElem->value[e]);
            newval = strdup(metadataElem->value[e]);
            // Add or update value
            updateValue(&lste_set, mnamebuf, DFNT_CHAR8, newlen, newval);
            newval = strdup(metadataElem->value[e]);
            updateValue(&cloud_set, mnamebuf, DFNT_CHAR8, newlen, newval);
        }
        metadataElem = metadataElem->next;
    }
    
    // Find the L2_CLOUD_Metadata group
    metadataGroup = getGroup(run_params, "L2_CLOUD_Metadata");
    if (metadataGroup == NULL)
        ERREXIT(8, "Missing <group> L2_CLOUD_Metadata in %s", run_params_filename);
    metadataElem = metadataGroup->head;
    while (metadataElem != NULL)
    {
        for (e = 0; e < metadataElem->count; e++)
        {
            // If multiple values, make unique names for each.
            if (e == 0)
                strcpy(mnamebuf, metadataElem->name);
            else
                sprintf(mnamebuf, "%s%d", metadataElem->name, e);
            newlen = strlen(metadataElem->value[e]);
            newval = strdup(metadataElem->value[e]);
            // Add or update value
            updateValue(&cloud_set, mnamebuf, DFNT_CHAR8, newlen, newval);
        }
        metadataElem = metadataElem->next;
    }
    
    // Find the L2_LSTE_Metadata group
    metadataGroup = getGroup(run_params, "L2_LSTE_Metadata");
    if (metadataGroup == NULL)
        ERREXIT(8, "Missing <group> L2_LSTE_Metadata in %s", run_params_filename);
    metadataElem = metadataGroup->head;
    while (metadataElem != NULL)
    {
        for (e = 0; e < metadataElem->count; e++)
        {
            // If multiple values, make unique names for each.
            if (e == 0)
                strcpy(mnamebuf, metadataElem->name);
            else
                sprintf(mnamebuf, "%s%d", metadataElem->name, e);
            newlen = strlen(metadataElem->value[e]);
            newval = strdup(metadataElem->value[e]);
            // Add or update value
            updateValue(&lste_set, mnamebuf, DFNT_CHAR8, newlen, newval);
        }
        metadataElem = metadataElem->next;
    }

    // ProductionDateTime from piped output of date command.
    newval = strdup(timestamp);
    newlen = strlen(newval);
    updateValue(&lste_set, "ProductionDateTime",
            DFNT_CHAR8, newlen, newval);
    newval = strdup(timestamp);
    updateValue(&cloud_set, "ProductionDateTime",
            DFNT_CHAR8, newlen, newval);

    // ProcessingEnvironment from piped output of uname -a
    newval = strdup(uname);
    newlen = strlen(newval);
    updateValue(&lste_set, "ProcessingEnvironment",
        DFNT_CHAR8, newlen, newval);
    newval = strdup(uname);
    updateValue(&cloud_set, "ProcessingEnvironment",
        DFNT_CHAR8, newlen, newval);

    // Add bounding coordinates
    updateValue(&lste_set, "NorthBoundingCoordinate",
            DFNT_FLOAT64, 1, dupfloat64(maxLat));
    updateValue(&lste_set, "SouthBoundingCoordinate",
            DFNT_FLOAT64, 1, dupfloat64(minLat));
    updateValue(&lste_set, "EastBoundingCoordinate",
            DFNT_FLOAT64, 1, dupfloat64(maxLon));
    updateValue(&lste_set, "WestBoundingCoordinate",
            DFNT_FLOAT64, 1, dupfloat64(minLon));
    updateValue(&cloud_set, "NorthBoundingCoordinate",
            DFNT_FLOAT64, 1, dupfloat64(maxLat));
    updateValue(&cloud_set, "SouthBoundingCoordinate",
            DFNT_FLOAT64, 1, dupfloat64(minLat));
    updateValue(&cloud_set, "EastBoundingCoordinate",
            DFNT_FLOAT64, 1, dupfloat64(maxLon));
    updateValue(&cloud_set, "WestBoundingCoordinate",
            DFNT_FLOAT64, 1, dupfloat64(minLon));

    // Get BuildID.
    elem = findGroupElement(run_config, "PrimaryExecutable", "BuildID");
    if (elem == NULL)
    {
	    paramstr = "Unspecified";
    }
    else
    {
        paramstr = getValue(elem, 0);
        if (!paramstr)
            ERREXIT(5, "Error reading value of BuildID from %s", run_config_filename);
    }
    newval = strdup(paramstr);
    newlen = strlen(newval);
    updateValue(&lste_set, "BuildID", DFNT_CHAR8, newlen, newval);
    newval = strdup(paramstr);
    updateValue(&cloud_set, "BuildID", DFNT_CHAR8, newlen, newval);

    // SceneID
    newval = strdup(scene_id);
    newlen = strlen(newval);
    updateValue(&lste_set, "SceneID", DFNT_CHAR8, newlen, newval);
    newval = strdup(scene_id);
    updateValue(&cloud_set, "SceneID", DFNT_CHAR8, newlen, newval);

    // orbit_number
    newval = strdup(orbit_number);
    newlen = strlen(newval);
    updateValue(&lste_set, "StartOrbitNumber", DFNT_CHAR8, newlen, newval);
    newval = strdup(orbit_number);
    updateValue(&cloud_set, "StartOrbitNumber", DFNT_CHAR8, newlen, newval);
    // TODO Determine if StopOrbitNumber is incremented or not.
    newval = strdup(orbit_number);
    updateValue(&lste_set, "StopOrbitNumber", DFNT_CHAR8, newlen, newval);
    newval = strdup(orbit_number);
    updateValue(&cloud_set, "StopOrbitNumber", DFNT_CHAR8, newlen, newval);

    // Determine if the input metadata has RangeBeginning time. If it does, 
    // use that and add elapsed seconds for ending time. If not, then use 
    // hour and minute from tes_date.
    TesDate beg_date;
    tes_copy_date(&beg_date, &tes_date);
    int beg_secs = 0;
    int beg_usec = 0;
    MetadataAttribute *rbtime = findAttribute(&lste_set, "RangeBeginningTime");
    if (rbtime)
    {
        sscanf(rbtime->value, "%02d:%02d:%02d.%d",
                &beg_date.hour, &beg_date.minute, &beg_secs, &beg_usec);
    }
    TesDate end_date;
    tes_copy_date(&end_date, &beg_date);
    int add_mins = (int)(elapsed_seconds / 60);
    int add_secs = (int)fmod(elapsed_seconds, 60.0);
    int end_secs = beg_secs + add_secs;
    if (end_secs > 60)
    {
        end_secs %= 60;
        add_mins++;
    }
    tes_add_minutes(&end_date, add_mins);

    // RangeBeginningData -- RangeEndingDate
    char dtbuf[256];
    sprintf(dtbuf, "%4d-%02d-%02d", 
        beg_date.year, beg_date.month, beg_date.day);
    newval = strdup(dtbuf);
    newlen = strlen(newval);
    updateValue(&lste_set, "RangeBeginningDate", DFNT_CHAR8, newlen, newval);
    newval = strdup(dtbuf);
    updateValue(&cloud_set, "RangeBeginningDate", DFNT_CHAR8, newlen, newval);
    sprintf(dtbuf, "%4d-%02d-%02d", 
        end_date.year, end_date.month, end_date.day);
    newval = strdup(dtbuf);
    newlen = strlen(newval);
    updateValue(&lste_set, "RangeEndingDate", DFNT_CHAR8, newlen, newval);
    newval = strdup(dtbuf);
    updateValue(&cloud_set, "RangeEndingDate", DFNT_CHAR8, newlen, newval);

    // RangeEndingDate -- RangeEndingTime
    sprintf(dtbuf, "%02d:%02d:%02d.%06d", beg_date.hour, beg_date.minute,
        beg_secs, beg_usec);
    newval = strdup(dtbuf);
    newlen = strlen(newval);
    updateValue(&lste_set, "RangeBeginningTime", DFNT_CHAR8, newlen, newval);
    newval = strdup(dtbuf);
    updateValue(&cloud_set, "RangeBeginningTime", DFNT_CHAR8, newlen, newval);
    sprintf(dtbuf, "%02d:%02d:%02d.%06d", end_date.hour, end_date.minute,
        end_secs, beg_usec);
    newval = strdup(dtbuf);
    newlen = strlen(newval);
    updateValue(&lste_set, "RangeEndingTime", DFNT_CHAR8, newlen, newval);
    newval = strdup(dtbuf);
    updateValue(&cloud_set, "RangeEndingTime", DFNT_CHAR8, newlen, newval);

    // If the Metadata doesn't include RegionID, make it up.
    MetadataAttribute *regionid = findAttribute(&lste_set, "RegionID");
    if (regionid == NULL)
    {
        newval = strdup("Unspecified");
        newlen = strlen(newval);
        updateValue(&lste_set, "RegionID", DFNT_CHAR8, newlen, newval);
        newval = strdup("Unspecified");
        updateValue(&cloud_set, "RegionID", DFNT_CHAR8, newlen, newval);
    }

    // Set other attributes from the PCF file.
    // Set the Granule ID using the output file name.
    // Begin looking from last slash if name contains slash(es).
    char *froot = strrchr(output_filename, '/');
    if (froot != NULL)
        froot++;
    else
        froot = output_filename;
    /* Change granule ID to full filename
    // Attempt to find the part that comes right after LSTE_.
    char *granid = strstr(froot, "LSTE_");
    if (granid != NULL)
        granid += 5;
    else
        granid = froot;
	granid = strdup(granid);; //Copy output_filename before modifying it.	
    // Strip off the filetype, e.g., ".h5".
    char *ftypstr = strrchr(granid, '.');
    if (ftypstr != NULL)
        *ftypstr = '\0';
    newval = strdup(granid);
*/
    newval = strdup(froot);
    newlen = strlen(newval);
    updateValue(&lste_set, "LocalGranuleID", DFNT_CHAR8, newlen, newval);
    froot = strrchr(cloud_filename, '/');
    if (froot != NULL)
        froot++;
    else
        froot = cloud_filename;
    newval = strdup(froot);
    newlen = strlen(newval);
    updateValue(&cloud_set, "LocalGranuleID", DFNT_CHAR8, newlen, newval);

    // DayNightFlag
    newval = strdup(DayNightFlag);
    newlen = strlen(newval);
    updateValue(&lste_set, "DayNightFlag", DFNT_CHAR8, newlen, newval);
    newval = strdup(DayNightFlag);
    updateValue(&cloud_set, "DayNightFlag", DFNT_CHAR8, newlen, newval);

    // Set the actual line/pixel size
    int32 image_lines = n_lines;
    int32 image_pixels = n_pixels;
    updateValue(&lste_set, "ImageLines",
        DFNT_INT32, 1, dupint32(image_lines));
    updateValue(&cloud_set, "ImageLines",
        DFNT_INT32, 1, dupint32(image_lines));
    updateValue(&lste_set, "ImagePixels",
        DFNT_INT32, 1, dupint32(image_pixels));
    updateValue(&cloud_set, "ImagePixels",
        DFNT_INT32, 1, dupint32(image_pixels));

    // Determine the line and pixel spacing
    double lat1 = mat_get(&vRAD.Lat, 0, 0);
    double lon1 = mat_get(&vRAD.Lon, 0, 0);
    double lat2 = mat_get(&vRAD.Lat, 0, n_pixels - 1);
    double lon2 = mat_get(&vRAD.Lon, 0, n_pixels - 1);
    double pixel_spacing = distance(lat1, lon1, lat2, lon2) / n_pixels;
    lat2 = mat_get(&vRAD.Lat, n_lines - 1, 0);
    lon2 = mat_get(&vRAD.Lon, n_lines - 1, 0);
    double line_spacing = distance(lat1, lon1, lat2, lon2) / n_lines;
    updateValue(&lste_set, "ImageLineSpacing", DFNT_FLOAT64, 1, 
            dupfloat64(line_spacing));
    updateValue(&cloud_set, "ImageLineSpacing", DFNT_FLOAT64, 1, 
            dupfloat64(line_spacing));
    updateValue(&lste_set, "ImagePixelSpacing", DFNT_FLOAT64, 1, 
            dupfloat64(pixel_spacing));
    updateValue(&cloud_set, "ImagePixelSpacing", DFNT_FLOAT64, 1, 
            dupfloat64(pixel_spacing));

    // Create the input pointer from GEO and RAD filenames.
    strcpy(input_pointer, GEO_filename);
    strcat(input_pointer, ", ");
    strcat(input_pointer, RAD_filename);
    newlen = strlen(input_pointer);
    newval = strdup(input_pointer);
    updateValue(&lste_set, "InputPointer",  DFNT_CHAR8, newlen, newval);
    newval = strdup(input_pointer);
    updateValue(&cloud_set, "InputPointer",  DFNT_CHAR8, newlen, newval);

    // SIS refers to the PSD that describes the L2 data products.
    const char* SISName = "Level 2 Product Specification Document (JPL D-94635)";
    const char* SISVersion = "Version 1";
    newlen = strlen(SISName);
    newval = strdup(SISName);
    updateValue(&lste_set, "SISName",  DFNT_CHAR8, newlen, newval);
    newval = strdup(SISName);
    updateValue(&cloud_set, "SISName",  DFNT_CHAR8, newlen, newval);
    newlen = strlen(SISVersion);
    newval = strdup(SISVersion);
    updateValue(&lste_set, "SISVersion",  DFNT_CHAR8, newlen, newval);
    newval = strdup(SISVersion);
    updateValue(&cloud_set, "SISVersion",  DFNT_CHAR8, newlen, newval);

    // Remove some metadata elements that are not listed in the PSD.
    MetadataAttribute* delId;
    delId = findAttribute(&lste_set, "ProcessingQAAttribute");
    if (delId != NULL) deleteAttribute(delId);
    delId = findAttribute(&cloud_set, "ProcessingQAAttribute");
    if (delId != NULL) deleteAttribute(delId);
    delId = findAttribute(&lste_set, "ProcessingQADescription");
    if (delId != NULL) deleteAttribute(delId);
    delId = findAttribute(&cloud_set, "ProcessingQADescription");
    if (delId != NULL) deleteAttribute(delId);

    // Some test data had BuildId -- delete that
    delId = findAttribute(&lste_set, "BuildId");
    if (delId != NULL) deleteAttribute(delId);
    delId = findAttribute(&cloud_set, "BuildId");
    if (delId != NULL) deleteAttribute(delId);
    updateValue(&lste_set, "FieldOfViewObstruction",
        DFNT_CHAR8, strlen(vRAD.FieldOfViewObstruction), 
        strdup(vRAD.FieldOfViewObstruction));
    updateValue(&cloud_set, "FieldOfViewObstruction",
        DFNT_CHAR8, strlen(vRAD.FieldOfViewObstruction), 
        strdup(vRAD.FieldOfViewObstruction));

    // Add LSTE only metadata
    updateValue(&lste_metadata, "AncillaryNWP",
            DFNT_CHAR8, strlen(get_nwp_name()), get_nwp_name());
    double good_frac = (int32)(0.5 + 100.0 *
            (double)count_good / (double)(n_lines * n_pixels));
    updateValue(&lste_metadata, "QAPercentGoodQuality",
            DFNT_INT32, 1, dupint32(good_frac));
    double avg = tot_good_lst / count_good;
    updateValue(&lste_metadata, "LSTGoodAvg",
            DFNT_FLOAT64, 1, dupfloat64(avg));
    char metadata_name[32];
    for (b = 0; b < n_channels; b++)
    {
        sprintf(metadata_name, "Emis%dGoodAvg", band[b] + 1);
        avg = tot_good_emis[b] / count_good;
        updateValue(&lste_metadata, metadata_name,
            DFNT_FLOAT64, 1, dupfloat64(avg));
    }

    // If only 3 channels, add zero good averages for other channels.
    if (n_channels == 3)
    {
        avg = 0.0;
        updateValue(&lste_metadata, "Emis1GoodAvg",
            DFNT_FLOAT64, 1, dupfloat64(avg));
        updateValue(&lste_metadata, "Emis3GoodAvg",
            DFNT_FLOAT64, 1, dupfloat64(avg));
    }

    // Add 3-band additions to metadata.
    updateValue(&lste_metadata, "NumberOfBands",
            DFNT_UINT8, 1, dupuint8((uint8)n_channels));
    float32 wavelength[6] = {0.0, 8.28, 8.78, 9.07, 10.52, 12.0};
    if (n_channels == 3)
    {
        wavelength[1] = 0.0;
        wavelength[3] = 0.0;
    }
    char *bsa = (char*)malloc(6*sizeof(float32));
    memcpy(bsa, wavelength, 6*sizeof(float32));
    updateValue(&lste_metadata, "BandSpecification",
            DFNT_FLOAT32, 6, bsa);

    // Create Metadata groups in L2_LSTE
    hid_t stdmet_group = H5Gcreate2(h5out, "/StandardMetadata", 
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (stdmet_group < 0)
        ERREXIT(19, "Failed to create group StandardMetadata in L2_LSTE", NONE);
    H5Gclose(stdmet_group);

    hid_t prodmet_group = H5Gcreate2(h5out, "/L2 LSTE Metadata", 
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (prodmet_group < 0)
        ERREXIT(19, "Failed to create group /L2 LSTE Metadata in L2_LSTE", NONE);
    H5Gclose(prodmet_group);

    // Create Metadata groups in L2_CLOUD
    hid_t cloudout = open_hdf5(cloud_filename, H5F_ACC_RDWR);
    stdmet_group = H5Gcreate2(cloudout, "/StandardMetadata", 
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (stdmet_group < 0)
        ERREXIT(17, "Failed to create group StandardMetadata in L2_CLOUD", NONE);
    H5Gclose(stdmet_group);

    prodmet_group = H5Gcreate2(cloudout, "/L2 CLOUD Metadata", 
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (prodmet_group < 0)
        ERREXIT(19, "Failed to create group /L2 CLOUD Metadata in L2_CLOUD", NONE);
    H5Gclose(prodmet_group);

    // Output the metadata to the output file StandardMetadata group.
    stat = writeMetadataGroup(&lste_set, h5out, "StandardMetadata");
    if (stat < 0)
        ERREXIT(19, "Failed to write StandardMetadata into L2_LSTE", NONE);
    stat = writeMetadataGroup(&cloud_set, cloudout, "StandardMetadata");
    if (stat < 0)
        ERREXIT(17, "Failed to write StandardMetadata into L2_CLOUD", NONE);

    // Write product specific metadata
    stat = writeMetadataGroup(&lste_metadata, h5out, "L2 LSTE Metadata");
    if (stat < 0)
        ERREXIT(19, "Failed to write product specific metadata into L2_LSTE", NONE);
    stat = writeMetadataGroup(&cloud_metadata, cloudout, "L2 CLOUD Metadata");
    if (stat < 0)
        ERREXIT(17, "Failed to write product specific metadata into L2_CLOUD", NONE);

    // Write metadata attributes to each dataset.

    // LST
    //int8 filli8 = 0;
    uint8 fillu8 = 0;
    uint16 fillu16 = 0;

    uint16 minLST = 7500;
    uint maxLST = 65535;
    stat = writeDatasetMetadataHdf5(h5out, "/SDS/LST",
            "Land Surface Temperature", "K", 
            DFNT_UINT16, &minLST, &maxLST,
            DFNT_UINT16, &fillu16, 0.02, 0.0);
    if (stat < 0)
        ERREXIT(19, "Failed to write attributes to dataset LST in %s",
            output_filename);

    // LST_err
    uint8 minLST_err = 1;
    uint8 maxLST_err = 255;
    stat = writeDatasetMetadataHdf5(h5out, "/SDS/LST_err",
            "Land Surface Temperature error", "K",
            DFNT_UINT8, &minLST_err, &maxLST_err,
            DFNT_UINT8, &fillu8, 0.04, 0.0);
    if (stat < 0)
        ERREXIT(19, "Failed to write attributes to dataset LST_err in %s",
            output_filename);

    // QC
    uint16 u16_0 = 0;
    uint16 u16_65535 = 65535;
    stat = writeDatasetMetadataHdf5(h5out, "/SDS/QC",
            "Quality Control for LST and emissivity", "n/a",
            DFNT_UINT16, &u16_0, &u16_65535,
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(19, "Failed to write attributes to dataset QC in %s",
            output_filename);

    // Emis*
    uint8 minEmis = 1;
    uint8 maxEmis = 255;
    char emis_ln[32];
    for (b = 0; b < n_channels; b++)
    {
        sprintf(emis_ds, "/SDS/Emis%d", band[b] + 1);
        sprintf(emis_ln, "Band %d Emissivity", band[b] + 1);
        stat = writeDatasetMetadataHdf5(h5out, emis_ds,
            emis_ln, "n/a", 
            DFNT_UINT8, &minEmis, &maxEmis,
            DFNT_UINT8, &fillu8, 0.002, 0.49);
        if (stat < 0)
            ERREXIT(19, "Failed to write attributes to dataset %s in %s",
                emis_ds, output_filename);
    }
    if (n_channels == 3)
    {
        // Add metadata to the dummy datasets
        stat = writeDatasetMetadataHdf5(h5out, "/SDS/Emis1",
            "Band 1 Emissivity", "n/a", 
            DFNT_UINT8, &minEmis, &maxEmis,
            DFNT_UINT8, &fillu8, 0.002, 0.49);
        if (stat < 0)
            ERREXIT(19, "Failed to write attributes to dataset /SDS/Emis1 in %s",
                emis_ds, output_filename);
        stat = writeDatasetMetadataHdf5(h5out, "/SDS/Emis3",
            "Band 3 Emissivity", "n/a", 
            DFNT_UINT8, &minEmis, &maxEmis,
            DFNT_UINT8, &fillu8, 0.002, 0.49);
        if (stat < 0)
            ERREXIT(19, "Failed to write attributes to dataset /SDS/Emis3 in %s",
                emis_ds, output_filename);
    }

    // EmisWB
    stat = writeDatasetMetadataHdf5(h5out, "/SDS/EmisWB",
            "Wide Band Emissivity", "n/a",
            DFNT_UINT8, &minEmis, &maxEmis,
            DFNT_UINT8, &fillu8, 0.002, 0.49);
    if (stat < 0)
        ERREXIT(19, "Failed to write attributes to dataset EmisWB in %s",
            output_filename);

    // Emis*_err
    uint16 minEmis_err = 1;
    uint16 maxEmis_err = 65535;
    for (b = 0; b < n_channels; b++)
    {
        sprintf(emis_ds, "/SDS/Emis%d_err", band[b] + 1);
        sprintf(emis_ln, "Band %d Emissivity error", band[b] + 1);
        stat = writeDatasetMetadataHdf5(h5out, emis_ds,
            emis_ln, "n/a",
            DFNT_UINT16, &minEmis_err, &maxEmis_err,
            DFNT_UINT16, &fillu16, 1.0e-4, 0.0);
        if (stat < 0)
            ERREXIT(19, "Failed to write attributes to dataset %s in %s",
                emis_ds, output_filename);
    }
    if (n_channels == 3)
    {
        // Add metadata to dummy _err datasets.
        stat = writeDatasetMetadataHdf5(h5out, "/SDS/Emis1_err",
            "Band 1 Emissivity error", "n/a",
            DFNT_UINT16, &minEmis_err, &maxEmis_err,
            DFNT_UINT16, &fillu16, 1.0e-4, 0.0);
        if (stat < 0)
            ERREXIT(19, "Failed to write attributes to dataset /SDS/Emis1_err in %s",
                emis_ds, output_filename);
        stat = writeDatasetMetadataHdf5(h5out, "/SDS/Emis3_err",
            "Band 3 Emissivity error", "n/a",
            DFNT_UINT16, &minEmis_err, &maxEmis_err,
            DFNT_UINT16, &fillu16, 1.0e-4, 0.0);
        if (stat < 0)
            ERREXIT(19, "Failed to write attributes to dataset /SDS/Emis3_err in %s",
                emis_ds, output_filename);
    }

    // PWV
    uint16 u16_1 = 1;
    stat = writeDatasetMetadataHdf5(h5out, "/SDS/PWV",
            "Precipitable Water Vapor", "cm",
            DFNT_UINT16, &u16_1, &u16_65535,
            DFNT_UINT16, &u16_0, 0.001, 0.0);
    if (stat < 0)
        ERREXIT(19, "Failed to write attributes to dataset PWV in %s",
            output_filename);

    // Add attributes to the output dataset.
    uint8 minCM = 0;
    uint8 maxCM = 1;
    uint8 fillCM = 255;
    stat = writeDatasetMetadataHdf5(h5out, "/SDS/cloud_mask",
            "Final cloud mask",
            "1 = cloudy",
            DFNT_UINT8, &minCM, &maxCM,
            DFNT_UINT8, &fillCM, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(19, "Failed to write attributes to dataset "
            "/SDS/cloud_mask in %s", output_filename);

    // water_mask
    uint8 pix0 = 0;
    uint8 pix1 = 1;
    stat = writeDatasetMetadataHdf5(h5out, "/SDS/water_mask",
            "water mask", "0=land 1=water",
            DFNT_UINT8, &pix0, &pix1,
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(19, "Failed to write attributes to dataset water_mask in %s",
            output_filename);

/* 2018-02-28 Removed oceanpix and View_angle which can be obtained from L1B_GEO
    // oceanpix
    uint8 pix0 = 0;
    uint8 pix1 = 1;
    stat = writeDatasetMetadataHdf5(h5out, "/SDS/oceanpix",
            "ocean pixels", "0=land 1=ocean",
            DFNT_UINT8, &pix0, &pix1,
            DFNT_NONE, NULL, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(19, "Failed to write attributes to dataset oeanpix in %s",
            output_filename);

    // View_angle
    uint8 minva = 0;
    uint8 maxva = 180;
    uint8 fillva = 255;
    stat = writeDatasetMetadataHdf5(h5out, "/SDS/View_angle",
            "Sensor Zenith", "degrees",
            DFNT_UINT8, &minva, &maxva,
            DFNT_UINT8, &fillva, 0.5, 0.0);
    if (stat < 0)
        ERREXIT(19, "Failed to write attributes to dataset View_angle in %s",
            output_filename);
*/

#ifdef PRODUCE_GEOLOCATION
    // Latitude
    float minlat_f32 = -90.0;
    float maxlat_f32 =  90.0;
    float fillf32 = -9999.0;
    stat = writeDatasetMetadataHdf5(h5out, "/LSTE/GeolocationFields/Latitude",
            "Latitude of each pixel", "degrees",
            DFNT_FLOAT32, &minlat_f32, &maxlat_f32,
            DFNT_FLOAT32, &fillf32, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(19, "Failed to write attributes to dataset Latitude in %s",
            output_filename);
    stat = writeDatasetMetadataHdf5(cloudout, "/CLOUD/GeolocationFields/Latitude",
            "Latitude of each pixel", "degrees",
            DFNT_FLOAT32, &minlat_f32, &maxlat_f32,
            DFNT_FLOAT32, &fillf32, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(17, "Failed to write attributes to dataset Latitude in %s",
            cloud_filename);

    // Longitude
    float minlon_f32 = -180.0;
    float maxlon_f32 =  180.0;
    stat = writeDatasetMetadataHdf5(h5out, "/LSTE/GeolocationFields/Longitude",
            "Longitude of each pixel", "degrees",
            DFNT_FLOAT32, &minlon_f32, &maxlon_f32,
            DFNT_FLOAT32, &fillf32, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(19, "Failed to write attributes to dataset Longitude in %s",
            output_filename);
    stat = writeDatasetMetadataHdf5(cloudout, "/CLOUD/GeolocationFields/Longitude",
            "Longitude of each pixel", "degrees",
            DFNT_FLOAT32, &minlon_f32, &maxlon_f32,
            DFNT_FLOAT32, &fillf32, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(17, "Failed to write attributes to dataset Longitude in %s",
            cloud_filename);
#endif

    // ========================================================
    //  Cleanup and Exit
    // ========================================================

    close_hdf5(cloudout);
    close_hdf5(h5out);

#ifdef ENABLE_TIMING
    printf("Output time: %.3f seconds\n",
            toseconds(elapsed(startoutput)));
    fflush(stdout);
#endif

    clear_rad(&vRAD);

#ifdef ENABLE_TIMING
    double runsecs = toseconds(elapsed(startmain));
    int runmins = runsecs / 60;
    runsecs -= runmins * 60;
    printf("L2_PGE runtime: %d minutes %.1f seconds\n",
            runmins, runsecs);
    fflush(stdout);
#endif

    INFO("ECOSTRESS L2 completed successfully.", NONE);
    return 0;
}

//--------------------------------------------------------
//  ATM RTTOV structure [SL]
//--------------------------------------------------------

// Initialize ATM structure, allocating arrays [SL]
void init_ATM_rttov(ATM_rttov *ATM, int imax, int kmax, int jmax)
{
    assert(ATM != NULL);

    // Set array sizes
    ATM->imax = imax;
    ATM->kmax = kmax;
    ATM->jmax = jmax;
    ATM->nmax = kmax * jmax;

    // Set size parameters to default values
    ATM->nprof = 0;
    ATM->nlev = 0;

    // Allocate pressure level vectors.
    ATM->T = (double **)calloc(ATM->imax, sizeof(double*));
    ATM->Q = (double **)calloc(ATM->imax, sizeof(double*));

    // Allocate 2d arrays, setting values to zeros

    ATM->T2 = (double *)calloc(ATM->kmax * ATM->jmax, sizeof(double));	// t2 = zeros(dm1,dm2)
    ATM->Q2 = (double *)calloc(ATM->kmax * ATM->jmax, sizeof(double));	// q2 = zeros(dm1,dm2)
    ATM->P =  (double *)calloc(ATM->nmax * ATM->imax, sizeof(double));	// REPLACE [SL] 
    // FROM: p = zeros(dm1,dm2,dm3)
    //       P = reshape(sATM.p,dm1*dm2,[],1)
    // TO: P[nmax][imax], where nmax = kmax*jmax

    ATM->PSurf_sp =  (double *)calloc(ATM->kmax * ATM->jmax, sizeof(double));    // sp = zeros(dm1,dm2)
    ATM->HSurf_el =  (double *)calloc(ATM->kmax * ATM->jmax, sizeof(double));    // el = zeros(dm1,dm2)
    ATM->TSurf_skt = (double *)calloc(ATM->kmax * ATM->jmax, sizeof(double));    // skt = zeros(dm1,dm2)
    ATM->Bemis = (double *)calloc(n_channels * ATM->kmax * ATM->jmax, sizeof(double));    // Bemis = zeros(n_channels,dm1,dm2)

    ATM->ST =        (double *)calloc(ATM->kmax * ATM->jmax, sizeof(double));     // lsm = zeros(dm1,dm2); ST = nwpATM.lsm(:)
    ATM->Lat_atmos = (double *)calloc(ATM->kmax * ATM->jmax, sizeof(double));     // lat = zeros(dm1,dm2)    
    ATM->Lon_atmos = (double *)calloc(ATM->kmax * ATM->jmax, sizeof(double));     // lat = zeros(dm1,dm2)
    ATM->Lat_samp = (double *)calloc(ATM->kmax * ATM->jmax, sizeof(double));     // lat = zeros(dm1,dm2)    
    ATM->Lon_samp = (double *)calloc(ATM->kmax * ATM->jmax, sizeof(double));     // lat = zeros(dm1,dm2)
    ATM->SatZen =    (double *)calloc(ATM->kmax * ATM->jmax, sizeof(double));    // SatZen = zeros(dm1,dm2)
    ATM->TCW = (double *)calloc(ATM->kmax * ATM->jmax, sizeof(double));

    // Allocate 3d arrays
    mat3d_init_zeros_rttov(ATM->T, ATM->imax, ATM->kmax, ATM->jmax);            // t = zeros(dm1,dm2,dm3)
    mat3d_init_zeros_rttov(ATM->Q, ATM->imax, ATM->kmax, ATM->jmax);            // q = zeros(dm1,dm2,dm3)

    // Allocate 1d arrays
    ATM->Tt = (double *)calloc(ATM->imax * ATM->kmax * ATM->jmax, sizeof(double));     // ADD [SL] 
    ATM->Qt = (double *)calloc(ATM->imax * ATM->kmax * ATM->jmax, sizeof(double));     // ADD [SL]

    // Check arrays

    assert(ATM != NULL); 
    assert(ATM->T2 != NULL);
    assert(ATM->Q2 != NULL);
    assert(ATM->P != NULL);
    assert(ATM->PSurf_sp != NULL);
    assert(ATM->HSurf_el != NULL);
    assert(ATM->TSurf_skt != NULL);
    assert(ATM->ST != NULL);
    assert(ATM->Lat_atmos != NULL);
    assert(ATM->Lon_atmos != NULL);
    assert(ATM->SatZen != NULL);
    assert(ATM->T != NULL);
    assert(ATM->Q != NULL);
    assert(ATM->Tt != NULL);
    assert(ATM->Qt != NULL);
    assert(ATM->Bemis != NULL);
}

// Clear ATM structure, deallocating arrays [SL]
void clear_ATM_rttov(ATM_rttov *ATM)
{
    assert(ATM != NULL);

    // Deallocate arrays

    if (ATM->T2 != NULL) free(ATM->T2);
    if (ATM->Q2 != NULL) free(ATM->Q2);
    if (ATM->P != NULL) free(ATM->P);

    if (ATM->PSurf_sp != NULL) free(ATM->PSurf_sp);
    if (ATM->HSurf_el != NULL) free(ATM->HSurf_el);
    if (ATM->TSurf_skt != NULL) free(ATM->TSurf_skt);

    if (ATM->ST != NULL) free(ATM->ST);
    if (ATM->Lat_atmos != NULL) free(ATM->Lat_atmos);
    if (ATM->Lon_atmos != NULL) free(ATM->Lon_atmos);
    if (ATM->Lat_samp != NULL) free(ATM->Lat_samp);
    if (ATM->Lon_samp != NULL) free(ATM->Lon_samp);
    if (ATM->SatZen != NULL) free(ATM->SatZen);
    if (ATM->TCW != NULL) free(ATM->TCW);
    if (ATM->Bemis != NULL) free(ATM->Bemis);

    mat3d_clear_rttov(ATM->T, ATM->imax);
    mat3d_clear_rttov(ATM->Q, ATM->imax);
    free(ATM->T);
    free(ATM->Q);

    if (ATM->Tt != NULL) free(ATM->Tt);
    if (ATM->Qt != NULL) free(ATM->Qt);
}

//-------------------------------------------------------------
// Read RTTOV Atmos
//-------------------------------------------------------------
void set_rttov_atmos(ATM_rttov *ATM,
    Mat3d *WV_mix, Mat3d *T_d, double Pres[], 
    Matrix *Senszen, Matrix *Psurf, 
    Matrix *Hsurf, Matrix *Lat, Matrix *Lon, Matrix *tcw_in,
    Matrix *Tskin, Matrix *t2, Matrix *q2)
{
    // t2, and q2 may be empty datasets (size1 == 0).
    // If t2 or q2 is empty, the results will be calculated.

    // Constants
    //const double md = 28.966;   // Md, molecular mass of dry air
    //const double mw = 18.016;   // Mv, molecular mass of water
    //const double wv_scale = 1000.0;  // 1e3 scaling factor
    //const double wv_coeff = (md / mw) * wv_scale;  // ADD [SL] (Md/Mw)*1e3
    const double eps_zero1 = 0.0001;                // ADD [SL] 0.0001
    const double min_tsurf = 90.0;

    // Array sizes
    int i = 0;
    int k = 0;
    int j = 0;

    int imax = 0;
    int kmax = 0;
    int jmax = 0;
    int nmax = 0;

    // Vectors
    Vector Tdata;   // vector Tdata[] = T_d[][k][j] matrix slice
    Vector Wdata;   // vector Wdata[] = WVpro[][k][j] matrix slice
    vec_init(&Tdata);
    vec_init(&Wdata);

    // Array 2d
    double *WVsurface = NULL;            // WVsurface, used for q2 assignment

    // Parameters
    bool bc1 = false;                    // flag for c value checks
    bool bc2 = false;

    double wv_surf = 0.0;                // WVsurf
    double t_surf = 0.0;                // ADD [SL]
    double t_mean = 0.0;                // ADD [SL]

    double p_val = 0.0;                    // ADD [SL]

    double psurf_mean = 0.0;            // ADD [SL]
    double hsurf_mean = 0.0;            // ADD [SL]
    double tsurf_mean = 0.0;            // ADD [SL]
    double t2_mean = 0.0;                // ADD [SL]


    // Copy array sizes from ATM, for simplifying the code only

    imax = ATM->imax;
    kmax = ATM->kmax;
    jmax = ATM->jmax;
    nmax = ATM->nmax;

    // Allocate arrays
    vec_alloc(&Tdata, imax);
    vec_alloc(&Wdata, imax);

    // ADD [SL]: explicit init for WVsurface
    WVsurface = (double *) calloc(kmax * jmax, sizeof(double));

    // Fill near surface levels
    for(k = 0; k < kmax; k++)     // 1:dm1
    {
        for(j = 0; j < jmax; j++)     // 1:dm2
        {
            mat3d_copy_to_vector_kj_rttov(Tdata.vals, k, j, T_d, imax, kmax, jmax);         // Tdata = squeeze(T_d(iq,jq,:))
            mat3d_copy_to_vector_kj_rttov(Wdata.vals, k, j, WV_mix, imax, kmax, jmax);     // WVdata = squeeze(WVpro(iq,jq,:))

            // corresponds to c = find(Tdata<0)
            bc1 = vect_is_all_negative_rttov(Tdata.vals, imax);

            if(bc1)    // numel(c)==dm3
            {
                // CHECK: All Tdata values < 0
                mat3d_set_kj_to_val(T_d, eps_zero1, k, j);// T_d(iq,jq,:) = 0.0001
                ATM->T2[k*jmax+j] = eps_zero1;         // t2(iq,jq) = 0.0001
                ATM->Q2[k*jmax+j] = eps_zero1;        // q2(iq,jq) = 0.0001
            }
            else
            {
                // CHECK: Non-negative Tdata value is found

                bc2 = vect_is_negative_found_rttov(Tdata.vals, imax); // checking for at least one Tdata[i] < 0 value

                if(bc2)    // ~isempty(c) & Tdata(c)<0
                {
                    // Repeat and avoid the negative tempr
                    vect_reset_negative_rttov(Tdata.vals, imax);                         // Tdata(c) = Tdata(c(1)-1)

                    // Set the surface WV
                    wv_surf = vect_get_last_positive_rttov(Wdata.vals, imax);            // iWVsurf = WVdata(WVdata >0) 
                    // WVsurf = iWVsurf(end)

                    t_surf = vect_get_last_positive_rttov(Tdata.vals, imax);                // iTSurf = Tdata(Tdata>0)
                    ATM->T2[k*jmax+j] = t_surf;                                        // t2(iq,jq) = iTSurf(end)
                    WVsurface[k*jmax+j] = wv_surf;                                    // WVsurface(iq,jq) = WVsurf    

                    mat3d_copy_from_vector_kj(Tdata.vals, k, j, T_d); // T_d(iq,jq,:) = Tdata
                }
                else
                {
                    // CHECK: All values are zero or positive

                    wv_surf = Wdata.vals[imax-1];            // WVsurf  = WVdata(end)
                    WVsurface[k*jmax+j] = wv_surf;        // WVsurface(iq,jq) = WVsurf

                    t_surf = Tdata.vals[imax-1];                // t2(iq,jq) = Tdata(end)
                    ATM->T2[k*jmax+j] = t_surf;
                }
            }
        }
    }

    vec_clear(&Tdata);
    vec_clear(&Wdata);

    mat3d_set_zero_to_val(WV_mix, 1.0);        // WVpro(WVpro==0)= 1

    mat3d_copy_mat(ATM->T, T_d, imax, kmax, jmax);                 // t = T_d;
    t_mean = mat3d_get_mean_more_eq_val_rttov(ATM->T, eps_zero1, imax, kmax, jmax);    // t(t<0| t==0.0001)= mean(t(t>0.0001))
    // [BF] mean cannot be less than hard minimum for the dataset.
    if (t_mean < hard_t_min)
    {
        WARNING("Retrieved_Temperature_Profile -- mean reset from %.2f to %.2f\n",
                t_mean, hard_t_min);
        t_mean = hard_t_min;
    }


    // q = WVpro // top to surface
    mat3d_copy_mat(ATM->Q, WV_mix, imax, kmax, jmax);
    // q2 = WVsurface
    if (q2->size1 != 0)
    {
        // q2 is supplied
        mat2d_copy_rttov(ATM->Q2, q2->vals, kmax, jmax);
    }
    else
    {
        // q2 is calculated
        mat2d_copy_rttov(ATM->Q2, WVsurface, kmax, jmax);
    }
    // q2(q2==0)= 0.0001
    mat2d_set_zero_to_val_rttov(ATM->Q2, eps_zero1, kmax, jmax);

    // Set size parameters // ADD COMMENT [SL]

    ATM->nprof = kmax * jmax;    // nprof = dm1*dm2

    ATM->nlev = imax;    // NOTE: REPLACE [SL] 
    // FROM: nlev(1:dm1,1:dm2) = dm3 
    //       sATM1(..., nlev, ...)
    //       nlev = nwpATM.nlev(1) 
    // TO: nlev as a single number, 
    //     where nlev = dm3 == imax

    // Pressure is in one dimension, make it three, by repmat along dim3
    //
    // NOTE: REPLACE [SL] 
    // FROM: p = flipdim(Pres,1)'
    //       p = repmat(p,[dm1,1,dm2])
    //       p = permute(p,[1,3,2])
    //       sATM1(..., p, ...)
    //       P = reshape(nwpATM.p,dm1*dm2,[],1)
    // TO: P[nmax][imax], where nmax = kmax*jmax
    //     P[][i] = Pres[imax - i - 1] reverse order of values,
    //     where each i-th column contains the same pressure value

    // [BF] For Merra, Pres is already flipdim.
    for(i = 0; i < imax; i++)
    {
        // [BF] p_val = Pres[imax-i-1];
        p_val = Pres[i];
        mat2d_set_j_column_to_val_rttov(ATM->P, p_val, i, nmax, imax);
    }

    // Getting rid of negative values

    // PSurf(PSurf<0)= mean(PSurf(PSurf>0))
    mat2d_copy_mat(ATM->PSurf_sp, Psurf, kmax, jmax);
    // sp = PSurf
    psurf_mean = mat2d_get_mean_more_zero_rttov(ATM->PSurf_sp, kmax, jmax); 
    mat2d_set_less_zero_to_val_rttov(ATM->PSurf_sp, psurf_mean, kmax, jmax);

    // HSurf(HSurf<0)= mean(HSurf(HSurf>0))
    mat2d_copy_mat(ATM->HSurf_el, Hsurf, kmax, jmax);
    // el = HSurf
    hsurf_mean = mat2d_get_mean_more_zero_rttov(ATM->HSurf_el, kmax, jmax);
    mat2d_set_less_zero_to_val_rttov(ATM->HSurf_el, hsurf_mean, kmax, jmax);

    // TSurf(TSurf<0)= mean(TSurf(TSurf>0))
    mat2d_copy_mat(ATM->TSurf_skt, Tskin, kmax, jmax);
    // skt = TSurf
    tsurf_mean = mat2d_get_mean_more_zero_rttov(ATM->TSurf_skt, kmax, jmax);
    // [BF] This is hard coded in rttov to 90.0.
    if (tsurf_mean < min_tsurf) tsurf_mean = min_tsurf;
    mat2d_set_less_minval_to_val_rttov(
            ATM->TSurf_skt, min_tsurf, tsurf_mean, kmax, jmax);

    // If actual t2 data is present, replace the calculated values.
    if (t2->size1 != 0)
        mat2d_copy_mat(ATM->T2, t2, kmax, jmax);
    // t2(t2== 0.0001)= mean(t2(t2> 0.0001))
    t2_mean = mat2d_get_mean_more_val_rttov(ATM->T2, eps_zero1, kmax, jmax);
    mat2d_set_eps_to_val_rttov(ATM->T2, eps_zero1, t2_mean, kmax, jmax);

    // Copy Lat, Lon, Senszen arrays
    // Note: Copied 2d matrices are later transposed to match MATLAB 1d output layout

    mat2d_copy_mat(ATM->Lat_atmos, Lat, kmax, jmax);  // Lat_atmos = nwpATM.lat(:)    
    mat2d_copy_mat(ATM->Lon_atmos, Lon, kmax, jmax);  // Lon_atmos = nwpATM.lon(:)
    mat2d_copy_mat(ATM->SatZen, Senszen, kmax, jmax); // SatZen = nwpATM.SatZen(:)

    mat2d_copy_mat(ATM->Lat_samp, Lat, kmax, jmax);   // Lat_atmos = nwpATM.lat(:)    
    mat2d_copy_mat(ATM->Lon_samp, Lon, kmax, jmax);   // Lon_atmos = nwpATM.lon(:)
    mat2d_copy_mat(ATM->TCW, tcw_in, kmax, jmax);

    // [BF] Merra data needs Bemis, but all values are set to "no emis", 
    // i.e., 1e-6
    for (i = 0; i < n_channels * jmax * kmax; i++)
    {
        ATM->Bemis[i] = 1e-6;
    }

    // NOTE [SL]: ST == lsm == array of zeros

    // Transpose matrices
    // 1. Note that matrices are transposed in 
    //    read_rttov (not in write_rttov) function
    // 2. Matching MATLAB layout order for 1d data output
    // 3. Using in-place transpose function
    //
    // Transpose, convert, T, Q 3d matrices defined as double *Mat[N]

    mat3d_transpose_rttov(ATM->Tt, ATM->T, imax, kmax, jmax);
    mat3d_transpose_rttov(ATM->Qt, ATM->Q, imax, kmax, jmax);

    // Transpose 2d matrices, defined as double *Mat
    //
    // UPDATE [SL]: no need to transpose ST, a matrix of zeros

    mat2d_transpose_rttov(ATM->P, nmax, imax);
    mat2d_transpose_rttov(ATM->PSurf_sp,  kmax, jmax);
    mat2d_transpose_rttov(ATM->T2,        kmax, jmax);
    mat2d_transpose_rttov(ATM->Q2,        kmax, jmax);
    mat2d_transpose_rttov(ATM->TSurf_skt, kmax, jmax);
    mat2d_transpose_rttov(ATM->ST,        kmax, jmax); 
    mat2d_transpose_rttov(ATM->HSurf_el,  kmax, jmax);
    mat2d_transpose_rttov(ATM->Lat_atmos, kmax, jmax);
    mat2d_transpose_rttov(ATM->Lon_atmos, kmax, jmax);
    mat2d_transpose_rttov(ATM->SatZen,    kmax, jmax);
    // [BF] Bemis is a constant array of 1e-6, so transposed = nontransposed
    // [BF] do not transpose TCW -- it is not an input to RTTOV
    // Deallocate arrays, without clearing ATM arrays

    if(WVsurface != NULL) free(WVsurface);

    // Note: don't clear ATM arrays
}

//--------------------------------------------------------
// RAD RTTOV structure  [SL]
//--------------------------------------------------------

// Set RAD structure, initializing and setting array values [SL]

void set_RAD_rttov(RAD_rttov *RAD, double *Lat, double *Lon, 
    double *SatZen, int kmax, int jmax)
{
    // Check input

    assert(RAD != NULL);
    assert(Lat != NULL);
    assert(Lon != NULL);
    assert(SatZen != NULL);

    // Set array sizes

    RAD->kmax = kmax;
    RAD->jmax = jmax;
    RAD->nmax = kmax * jmax;

    // Allocate 2d arrays, setting values to zeros

    RAD->Lat =    (double *)calloc(RAD->nmax, sizeof(double));
    RAD->Lon =    (double *)calloc(RAD->nmax, sizeof(double));
    RAD->SatZen = (double *)calloc(RAD->nmax, sizeof(double));
    // [BF] OBE MErra: RAD->QC =     (uint8 *)calloc(RAD->nmax, sizeof(uint8));

    // Check that arrays are allocated

    if(RAD == NULL || RAD->Lat == NULL || 
        RAD->Lon == NULL || RAD->SatZen == NULL)
    {
        ERROR("cannot allocate RAD arrays", "set_RAD_rttov");
        return;
    }

    // Copy data from input arrays

    mat2d_copy_rttov(RAD->Lat, Lat, RAD->kmax, RAD->jmax);
    mat2d_copy_rttov(RAD->Lon, Lon, RAD->kmax, RAD->jmax);
    mat2d_copy_rttov(RAD->SatZen, SatZen, RAD->kmax, RAD->jmax);
    // [BF] OBE Merra: mat2d_copy_uint8_rttov(RAD->QC, BadPixel, RAD->kmax, RAD->jmax);
}

// Clear RAD structure, deallocating arrays [SL]

void clear_RAD_rttov(RAD_rttov *RAD)
{
    // Check input

    assert(RAD != NULL);

    // Deallocate 2d arrays

    if(RAD->Lat != NULL) free(RAD->Lat);
    if(RAD->Lon != NULL) free(RAD->Lon);
    if(RAD->SatZen != NULL) free(RAD->SatZen);
    //if(RAD->QC != NULL) free(RAD->QC);
}

//--------------------------------------------------------
// RadOut RTTOV structure  [SL]
//--------------------------------------------------------

// Initialize RadOut structure, allocating arrays [SL]
// 1. Array size nprof = ATM->nprof
// 2. Array size nvars = pre-defined number of variables

void init_RadOut_rttov(RadOut_rttov *RadOut, int nprof)
{
    assert(RadOut != NULL);

    // Set array sizes
    RadOut->nprof = nprof;

    // Allocate arrays for each channel
    int ichan;
    RadOut->TB_rttov = (double**)malloc(n_channels * sizeof(double*));
    RadOut->RadTot_rttov = (double**)malloc(n_channels * sizeof(double*));
    RadOut->RadUp_rttov = (double**)malloc(n_channels * sizeof(double*));
    RadOut->RadRefDn_rttov = (double**)malloc(n_channels * sizeof(double*));
    RadOut->Trans_rttov = (double**)malloc(n_channels * sizeof(double*));

    for (ichan = 0; ichan < n_channels; ichan++)
    {
        RadOut->TB_rttov[ichan] = (double*)calloc(RadOut->nprof, sizeof(double));
        RadOut->RadTot_rttov[ichan] = (double*)calloc(RadOut->nprof, sizeof(double));
        RadOut->RadUp_rttov[ichan] = (double*)calloc(RadOut->nprof, sizeof(double));
        RadOut->RadRefDn_rttov[ichan] = (double*)calloc(RadOut->nprof, sizeof(double));
        RadOut->Trans_rttov[ichan] = (double*)calloc(RadOut->nprof, sizeof(double));
    }
}

// Clear RadOut structure, deallocating arrays [SL]

void clear_RadOut_rttov(RadOut_rttov *RadOut)
{
    assert(RadOut != NULL);

    // Deallocate arrays
    int ichan;
    for (ichan = 0; ichan < n_channels; ichan++)
    {
        free(RadOut->TB_rttov[ichan]);
        free(RadOut->RadTot_rttov[ichan]);
        free(RadOut->RadUp_rttov[ichan]);
        free(RadOut->RadRefDn_rttov[ichan]);
        free(RadOut->Trans_rttov[ichan]);
    }
    free(RadOut->TB_rttov);
    free(RadOut->RadTot_rttov);
    free(RadOut->RadUp_rttov);
    free(RadOut->RadRefDn_rttov);
    free(RadOut->Trans_rttov);
}

//--------------------------------------------------------
// RTM RTTOV structure  [SL]
//--------------------------------------------------------

// Initialize RTM structure, allocating arrays [SL]

void init_RTM_rttov(RTM_rttov *RTM, int kmax, int jmax)
{
    assert(RTM != NULL);

    // Set array sizes

    RTM->kmax = kmax;
    RTM->jmax = jmax;
    RTM->nmax = kmax * jmax;

    RTM->RadUp_rt = (double**)malloc(n_channels * sizeof(double*));
    assert(RTM->RadUp_rt != NULL);
    RTM->RadRefDn_rt = (double**)malloc(n_channels * sizeof(double*));
    assert(RTM->RadRefDn_rt != NULL);
    RTM->Trans_rt = (double**)malloc(n_channels * sizeof(double*));
    assert(RTM->Trans_rt != NULL);
    
    int i;
    for (i = 0; i < n_channels; i++)
    {
        RTM->RadUp_rt[i] = (double *)calloc(RTM->nmax, sizeof(double));
        assert(RTM->RadUp_rt[i] != NULL);
        RTM->RadRefDn_rt[i] = (double *)calloc(RTM->nmax, sizeof(double));
        assert(RTM->RadRefDn_rt[i] != NULL);
        RTM->Trans_rt[i] = (double *)calloc(RTM->nmax, sizeof(double));
        assert(RTM->Trans_rt[i] != NULL);
    }

    RTM->PWV = (double *)calloc(RTM->nmax, sizeof(double));
}

// Clear RTM structure, deallocating arrays [SL]

void clear_RTM_rttov(RTM_rttov *RTM)
{
    assert(RTM != NULL);

    int i;
    for (i = 0; i < n_channels; i++)
    {
        free(RTM->RadUp_rt[i]);
        free(RTM->RadRefDn_rt[i]);
        free(RTM->Trans_rt[i]);
    }
    free(RTM->RadUp_rt);
    free(RTM->RadRefDn_rt);
    free(RTM->Trans_rt);

    free(RTM->PWV);
}

//-------------------------------------------------------------------
// Write RTTOV Atmos, writing binary profile input to RTTOV [SL]
//
// 1. Includes writing RTTOV input binary file
// 2. Corresponding MATLAB functions: write_atoms_rttov
// 3. Input arrays and sizes are from ATM data structure
//    - Assuming that matrices are already transposed
//      in read_rttov function, matching MATLAB layout
// 4. Input WVS (wvs_case) is defined in MATLAB: runRTTOVunixbin
//    - Case 1: WVS = 0
//            - Scaling Lon_atmos 
//          - Not scaling Q, Qt, Q2
//    - Case 2: WVS = 1
//          - Not scaling Lon Atmos
//          - Scaling Q, Qt, Q2  (scaling both 3d and 1d arrays Q and Qt ) 
//-------------------------------------------------------------------

void write_rttov_atmos(const char *fpath_rttov_profile, ATM_rttov *ATM, 
        int wvs_case)
{
    // Constants
    const double lon_coeff = 360.0;        // ADD [SL] Lon_atoms modulus coefficient
    const double wvs_coeff = 0.7;        // ADD [SL] WVS Q and Qsurf scaling coefficient

    // Check inputs
    assert(fpath_rttov_profile != NULL);
    assert(fpath_rttov_profile[0] != '\0');
    assert(ATM != NULL);
    assert(wvs_case == 0 || wvs_case == 1);

    // Scale ATM arrays according WVS case [SL]
    //
    // Note: Assuming that write_rttov with WVS = 0 is called before WVS = 1

    if(wvs_case == 0)
    {
        // Lon to 0-360
        //
        // Note: scaling Lon array for WVS case = 0 only
        // i.e. only for the 1st write RTTOV profile run

        vect_mod_rttov(ATM->Lon_atmos, lon_coeff, ATM->nmax);    // Lon_atmos=mod(Lon_atmos,360) 
    }
    else // WVS ==1
    {
        // Scale Q and Qsurf according to WVS settings
        //
        // Note: scaling all of Scaling Q, Qt, Q2, 3d, 1d, and 2d matrices

        mat3d_scale_rttov(ATM->Q, wvs_coeff, ATM->imax, ATM->kmax, ATM->jmax);    // Q = Q*0.7 // UPDATE [SL] scale Qt only
        vect_scale_rttov(ATM->Qt, wvs_coeff, ATM->imax * ATM->nmax);            // ADD [SL] scaling Qt in addition to scaling Q
        mat2d_scale_rttov(ATM->Q2, wvs_coeff, ATM->kmax, ATM->jmax);            // Qsurf = Qsurf*0.7
    }

    // Write RTTOV input binary profile file

    write_rttov_profile(fpath_rttov_profile, ATM);

    // Note: don't clear ATM arrays
}

//-------------------------------------------------------------------
// Run RTTOV shell script with binary profile as input to RTTOV [SL]
//-------------------------------------------------------------------

void run_rttov_script(const char *script, const char *exe, const char *coef)
{
    fflush(stdout);
    fflush(stderr);

    int sresult = 0;    // returned value of system call

    // CHECK: UPDATE [SL]
    //
    // Set environment variables for the GFORTRAN
    //
    // MATLAB:
    // setenv('GFORTRAN_STDIN_UNIT', '5');
    // setenv('GFORTRAN_STDOUT_UNIT', '6');
    // setenv('GFORTRAN_STDERR_UNIT', '0');

    // CHECK: UPDATE [SL]
    //
    // Make the sh file executable
    //
    // MATLAB: 
    // !chmod +x \\LAKEALMANOR\latyshev\rttov112\rttov_test\run_mod07_fwd.sh

    // Run RTTOV shell script
    //
    // MATLAB: system('~/rttov112/rttov_test/run_mod07_fwd.sh');
    // NOTE: 
    // 1. Not checking sresult value, since it might be system-dependant
    // 2. Not assigning value returned by the system results in compiler warning

    char script_command[PATH_MAX];
    strcpy(script_command, script);
    strcat(script_command, " ");
    strcat(script_command, exe); // arg 1 = name of executable 
    strcat(script_command, " ");
    strcat(script_command, coef); // arg 2 = coefficient file name

    sresult = system(script_command);

    if (sresult != 0)
    {
        ERREXIT(20, "%s exited with error code %d", script_command, sresult);
    }

    // CHECK: UPDATE [SL]
    //
    // Back to default environment variables
    //
    // MATLAB:
    // setenv('GFORTRAN_STDIN_UNIT', '-1');
    // setenv('GFORTRAN_STDOUT_UNIT', '-1');
    // setenv('GFORTRAN_STDERR_UNIT', '-1');
}

//-------------------------------------------
// Read and interpret RTTOV data output [SL]
//-------------------------------------------

void read_interp_rttov(const char* fpath_rttov_data, RTM_rttov *RTM, 
        RAD_rttov *RAD, ATM_rttov *ATM)
{
    FILE* fp = NULL;        // rad_out.dat data file
    RadOut_rttov RadOut;    // structure for cropped rad_out.dat data

    int ichan;              // channel index
    int n = 0;              // number of lines read from the text file
    int sresult = 0;        // result returned by fscanf
    int channel_no;         // Variable to read channel numbers from radout
    double radout_value;    // Variable to read double files from the radout file
#ifdef DEBUG_RTTOV
    int dbg_linno = 0;
#endif

    // Initialize RadOut structure

    init_RadOut_rttov(&RadOut, ATM->nprof);

    // Open data file
    fp = fopen(fpath_rttov_data, "r"); // r for reading text data

    if (fp == NULL)
        ERREXIT(21, "%s cannot open rad_out.dat", "read_interp_rttov");

    // Read Rad_out vector from text (table) data file
    // Text File LineFormat: 
    //  Chan1 ... ChanN TB1 ... TBN RadUp1 ... RadupN etc.

    n = 0;    // reset line number // Note: n is in the range [0, nprof]

    while (!feof(fp) && n < ATM->nprof) // CHECK: EOF on the system
    {
#ifdef DEBUG_RTTOV
        dbg_linno++;
#endif
        // Read a single line of table text file
        // UPDATE: CHECK: value of sresult (for error checking)

        // Discard the channel numbers at the start of each line.
        for (ichan = 0; ichan < n_channels; ichan++)
        {
            sresult = fscanf(fp, "%d", &channel_no);
#ifdef DEBUG_RTTOV
            if (channel_no != (band[ichan] + 1))
                printf("rad_out.dat line %d expecting channel %d got %d\n",
                       dbg_linno, band[ichan]+1, channel_no);
#endif
        }
        // Could have hit EOF
        if (feof(fp)) break;
        
        // Read and store TB data
        for (ichan = 0; ichan < n_channels; ichan++)
        {
            sresult = fscanf(fp, "%lf", &radout_value);
            if (sresult < 1)
                ERREXIT(22, "rad_out.dat TB input failure on line %d channel %d", 
                        n, ichan);
            RadOut.TB_rttov[ichan][n] = radout_value;
#ifdef DEBUG_RTTOV
            if (radout_value < 150.0 || radout_value > 350.0)
                printf("rad_out.dat line %d channel=%d TB=%.1f\n",
                       dbg_linno, band[ichan]+1, radout_value);
#endif
        }

        // Read and store RadTot data
        for (ichan = 0; ichan < n_channels; ichan++)
        {
            sresult = fscanf(fp, "%lf", &radout_value);
            if (sresult < 1)
                ERREXIT(22, "rad_out.dat RadTot input failure on line %d channel %d", 
                        n, ichan);
            RadOut.RadTot_rttov[ichan][n] = radout_value;
#ifdef DEBUG_RTTOV
            if (radout_value < 0.5 || radout_value > 15.0)
                printf("rad_out.dat line %d channel=%d RadTot=%.1f\n",
                       dbg_linno, band[ichan]+1, radout_value);
#endif
        }

        // Read and store RadUp data
        for (ichan = 0; ichan < n_channels; ichan++)
        {
            sresult = fscanf(fp, "%lf", &radout_value);
            if (sresult < 1)
                ERREXIT(22, "rad_out.dat RadUp input failure on line %d channel %d", 
                        n, ichan);
            RadOut.RadUp_rttov[ichan][n] = radout_value;
#ifdef DEBUG_RTTOV
            if (radout_value < 0.5 || radout_value > 15.0)
                printf("rad_out.dat line %d channel=%d RadUp=%.1f\n",
                       dbg_linno, band[ichan]+1, radout_value);
            if (ichan == 1 && n == 1250)
            {
                printf("RadUp(1,0,1250)=%.4f\n", radout_value);
            }
#endif
        }

        // Read and store RadRefDn data
        for (ichan = 0; ichan < n_channels; ichan++)
        {
            sresult = fscanf(fp, "%lf", &radout_value);
            if (sresult < 1)
                ERREXIT(22, "rad_out.dat RadRefDn input failure on line %d channel %d", 
                        n, ichan);
            RadOut.RadRefDn_rttov[ichan][n] = radout_value;
#ifdef DEBUG_RTTOV
            if (radout_value < 0.5 || radout_value > 15.0)
                printf("rad_out.dat line %d channel=%d RadRefDn=%.1f\n",
                       dbg_linno, band[ichan]+1, radout_value);
#endif
        }

        // Read and store Trans data
        for (ichan = 0; ichan < n_channels; ichan++)
        {
            sresult = fscanf(fp, "%lf", &radout_value);
            if (sresult < 1)
                ERREXIT(22, "rad_out.dat Trans input failure on line %d channel %d", 
                        n, ichan);
            RadOut.Trans_rttov[ichan][n] = radout_value;
#ifdef DEBUG_RTTOV
            if (radout_value < 0.85 || radout_value > 1.0)
                printf("rad_out.dat line %d channel=%d RadTRans=%.1f\n",
                       dbg_linno, band[ichan]+1, radout_value);
#endif
        }

        n++; // Increment line number

    }

    // Close data file
    fclose(fp);

    if (n != ATM->nprof)
        WARNING("%s expected %d lines, contained %d lines",
              fpath_rttov_data, ATM->nprof, n);

    // ================ RTTOV RADCON LUT ==============
    // Get number of lines
    char flut_name[PATH_MAX];
    strcpy(flut_name, OSP_dir);
    strcat(flut_name, rttov_radcon_filename);
    FILE *flut = fopen(flut_name, "r");
    if (!flut)
    {
        perror(flut_name);
        ERREXIT(23, "RttovRadconLUT open error on file ", flut_name);
    }
    
    int nlines_radconlut = 0;
    int ch;
    while (EOF != (ch=getc(flut)))
    {
        if (ch == '\n') nlines_radconlut++;
    }
    assert(nlines_radconlut > 1);
    rewind(flut);

    // Create space for the LUTdata
    double R_micron[MAX_CHANNELS][nlines_radconlut];
    double R_cm[MAX_CHANNELS][nlines_radconlut];
    double lut1;
    
    // Read the LUT data
    int fin;
    for (fin = 0; fin < nlines_radconlut; fin++)
    {
        sresult = fscanf(flut, "%lf", &lut1); // Discard the 1st value
        // For LUT file, need to read all channels.
        for (ichan = 0; ichan < MAX_CHANNELS; ichan++)
        {
            sresult = fscanf(flut, "%lf", &R_cm[ichan][fin]);
            if (sresult < 1)
                ERREXIT(24, "cm field input error in %s at channel %d on line %d",
                    rttov_radcon_filename, ichan, fin);
        }
        for (ichan = 0; ichan < MAX_CHANNELS; ichan++)
        {
            sresult = fscanf(flut, "%lf", &R_micron[ichan][fin]);
            if (sresult < 1)
                ERREXIT(24, "micron field input error in %s at channel %d on line %d",
                    rttov_radcon_filename, ichan, fin);
        }
    }
    fclose(flut);

    // [BF] Interpolate RadRefRn and RadUp. 
    int ilut; // Use band to get index into LUT channel (column).
    for (ichan = 0; ichan < n_channels; ichan++)
    {
        ilut = band[ichan];
        for (n = 0; n < ATM->nprof; n++)
        {
            RadOut.RadRefDn_rttov[ichan][n] = 
                interp1d_npts(R_cm[ilut], R_micron[ilut], 
                    RadOut.RadRefDn_rttov[ichan][n], nlines_radconlut);
            RadOut.RadUp_rttov[ichan][n] = 
                interp1d_npts(R_cm[ilut], R_micron[ilut], 
                    RadOut.RadUp_rttov[ichan][n], nlines_radconlut);
        }
    }

    //
    // ================ end LUT ==============

    // Transpose rows/cols of rad_out data
    for (ichan = 0; ichan < n_channels; ichan++)
    {
        mat2d_transpose_rttov(RadOut.RadUp_rttov[ichan], ATM->jmax, ATM->kmax);
        mat2d_transpose_rttov(RadOut.RadRefDn_rttov[ichan], ATM->jmax, ATM->kmax);
        mat2d_transpose_rttov(RadOut.Trans_rttov[ichan], ATM->jmax, ATM->kmax);
    }

    // Use interp2 to map multiple datasets to the granule grid.
    const int NDS = n_channels * 3 + 1;
    double *inds[NDS];
    double *outds[NDS];
    int i_nds = 0;
    for (ichan = 0; ichan < n_channels; ichan++)
    {
        inds[i_nds] = RadOut.RadUp_rttov[ichan];
        outds[i_nds] = RTM->RadUp_rt[ichan];
        i_nds++;

        inds[i_nds] = RadOut.RadRefDn_rttov[ichan];
        outds[i_nds] = RTM->RadRefDn_rt[ichan];
        i_nds++;
        
        inds[i_nds] = RadOut.Trans_rttov[ichan];
        outds[i_nds] = RTM->Trans_rt[ichan];
        i_nds++;
    }

    inds[i_nds] = ATM->TCW;
    outds[i_nds] = RTM->PWV;
    i_nds++;

    multi_interp2(ATM->Lat_samp, ATM->Lon_samp, ATM->kmax, ATM->jmax, 
        RAD->Lat, RAD->Lon, RAD->nmax, inds, outds, i_nds);

// @TODO add some DEBUG to look at the input and output values.

    // Done with cropped data
    clear_RadOut_rttov(&RadOut);
}

// Write single integer as float to binary file
// 1. Assuming that file is already open for writing
// 2. Converting integer value to float

void write_val_as_float_to_binary_rttov(FILE *fp, int val_int)
{
    size_t w_count = 0;                // number of written elements
    float val = (float) val_int;    // convert integer to float

    if(fp == NULL)
    {
        ERROR("write_val_as_float_to_binary_file", NONE);
        return;
    }

    // Write value to binary file
    w_count = fwrite(&val, sizeof(float), 1, fp); // Note: 1 is for single value

    if(w_count != 1)
    {
        ERROR("writing value to file in %s", "write_val_as_float_to_binary_file");
        return;
    }
}

// Write array to file fp in binary format as float (single values)
// 1. Assuming that file is already open for writing
// 2. Converting double to float array
// 3. Input imax = number of vector elements

void write_vector_as_float_to_binary_rttov(FILE *fp, double *vect_d, int imax)
{
    int i = 0;
    float *vect = NULL;
    size_t imax_count = (size_t) imax;		// number of elements to write
    size_t w_count = 0;						// number of written elements

    assert(fp != NULL);
    assert(vect_d != NULL);

    // Allocate float array, setting values to zeros
    vect = (float *) calloc(imax, sizeof(float)); 
    assert(vect != NULL);

    // Copy, converting double-to-float values
    for(i = 0; i < imax; i++)
    {
        vect[i] = (float) vect_d[i];	// copy values
    }

    // Write vector to binary file
    w_count = fwrite(vect, sizeof(float), imax_count, fp);

    if(w_count != imax_count)
        ERROR("wrong number of values for prof_in.bin: expected %d, found %d ",
                imax_count, w_count);

    // Deallocate float array
    free(vect);
}

// Write RTTOV input binary profile file

void write_rttov_profile(const char *fpath, ATM_rttov *ATM)
{
    // File pointer 
    FILE *fp = NULL;

    assert(fpath != NULL);
    assert(ATM != NULL);

    // Open file
    fp = fopen(fpath, "wb"); // wb for writing to binary

    if(fp == NULL)
        ERREXIT(32, "cannot open output file: %s", fpath);

    // Write array n sizes as single values

    write_val_as_float_to_binary_rttov(fp, ATM->nprof);
    write_val_as_float_to_binary_rttov(fp, ATM->nlev);

    // Write array values in binary format

    write_vector_as_float_to_binary_rttov(fp, ATM->P,         ATM->nmax * ATM->imax);
    write_vector_as_float_to_binary_rttov(fp, ATM->Tt,        ATM->nmax * ATM->imax);
    write_vector_as_float_to_binary_rttov(fp, ATM->Qt,        ATM->nmax * ATM->imax);
    write_vector_as_float_to_binary_rttov(fp, ATM->PSurf_sp,  ATM->nmax);
    write_vector_as_float_to_binary_rttov(fp, ATM->T2,        ATM->nmax);
    write_vector_as_float_to_binary_rttov(fp, ATM->Q2,        ATM->nmax);
    write_vector_as_float_to_binary_rttov(fp, ATM->TSurf_skt, ATM->nmax);
    write_vector_as_float_to_binary_rttov(fp, ATM->ST,        ATM->nmax);
    write_vector_as_float_to_binary_rttov(fp, ATM->HSurf_el,  ATM->nmax);
    write_vector_as_float_to_binary_rttov(fp, ATM->SatZen,    ATM->nmax);
    write_vector_as_float_to_binary_rttov(fp, ATM->Lat_atmos, ATM->nmax);
    write_vector_as_float_to_binary_rttov(fp, ATM->Lon_atmos, ATM->nmax);
    write_vector_as_float_to_binary_rttov(fp, ATM->Bemis,     ATM->nmax * n_channels);

    // Close file
    fclose(fp);
}
 
/**
 * @brief Planck function
 * 
 * @param[in] emax      Emax parameter
 * @param[in] surfradin surfradi{1:N}{x,y}
 * @param[in] skyradin  skyr{1:N}{x,y}
 * @param[in] it        number of iterations
 * @param[out] e
 * @param[out] Tnem
 * @param[out] Re
 * @param[out] kiter
 * @return QA
 */
int NEM_planck(
    // inputs:
    double emax, double surfradin[], double skyradin[], int it,
    // outputs:
    double *e, double *Tnem, double *Re, int *kiter)
    /* return value is QA */
{
#ifdef DEBUG_PIXEL
    if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
    {
        printf("DEBUG(%d,%d) NEM_planck emax=%.4f it=%d\n"
                "        surfradin=(", iline, ipixel, emax, it);
        int idb;
        for (idb = 0; idb < n_channels; idb++)
        {
            if (idb > 0) printf(", ");
            printf("%.4f", surfradin[idb]);
        }
        printf(")\n        skyradin =(");
        for (idb = 0; idb < n_channels; idb++)
        {
            if (idb > 0) printf(", ");
            printf("%.4f", skyradin[idb]);
        }
        printf(")\n");
    }
#endif
    double R[MAX_CHANNELS], Rold[MAX_CHANNELS], T[MAX_CHANNELS], Te[MAX_CHANNELS], B, diff[MAX_CHANNELS];
    int iband, dcon, ddiv,k; 

    assert(nlut_lines > 1); // LUT table must be preloaded with at least 2 entries.
    /// @TODO: make a function that returns a reference and loads LUT on first request.

    *Tnem = nan1;
    for(iband=0; iband<n_channels; iband++) { 
        R[iband] = nan1;
        T[iband] = nan1;
        if (is_nan(surfradin[iband])) 
        {
#ifdef DEBUG_PIXEL
            if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
            {
                printf("        QA==1: surfradin[%d] is nan\n", iband);
            }
#endif
            return 1;
        }
        if (is_nan(skyradin[iband])) 
        {
#ifdef DEBUG_PIXEL
            if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
            {
                printf("        QA==1: skyradin[%d] is nan\n", iband);
            }
#endif
            return 1;
        }
        R[iband] = surfradin[iband]-(1-emax)*skyradin[iband];
        // Replaced with LUT
        // wave5 = pow(wave[iband],5);
        // if(((c1*emax)/(PI*R[iband]*wave5)+1) < 0) {
#ifdef ERREXIT_MSG
        //     INFO("img arg: iband = %d, c1 = %g, emax = %g, R = %g, wave5 = %g\n",
        //         iband,c1,emax,R[iband],wave5);
#endif
        // return(-1);
        // }
        //T[iband] = (log((c1*emax)/(PI*R[iband]*wave5)+1)); 
        //if(fabs(T[iband]) > 0) T[iband] = (c2/wave[iband])/T[iband]; 
        //else T[iband] = nan1;
        
        // LUT for T
        T[iband] = interp1d_npts(lut[band[iband]+1], lut[0], R[iband]/emax, nlut_lines);
        
        /* Tnem is the maximum brightness temperature */
        if(T[iband] > *Tnem) *Tnem = T[iband]; 
    }  
    for(iband=0; iband<n_channels; iband++) { 
        // wave5 = pow(wave[iband],5);
        //  Compute blackbody radiance using Tnem
        //B = c1/(wave5*PI*(pow(E,(c2/(wave[iband]*(*Tnem))))-1));
        B = interp1d_npts(lut[0], lut[band[iband]+1], *Tnem, nlut_lines);
        //  Compute emissivity
        e[iband] = R[iband]/B;
    }

    /* Iteratively correct for downwelling radiance */
    for(iband=0; iband<n_channels; iband++) Rold[iband] = R[iband];
    for(k=1; k<=it; k++) {
        *kiter = k;
#ifdef DEBUG_PIXEL
        if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
        {
            printf("        Iteration %d\n", k);
        }
#endif
        if(k == it) 
        {
#ifdef DEBUG_PIXEL
            if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
            {
                printf("        QA==1: Planck max of %d reached.\n", it);
            }
#endif
            return(1);
        }
        *Tnem = nan1;
        for(iband=0; iband<n_channels; iband++) { 
            //Re[iband] = nan1; 
            //Te[iband] = nan1;
            Re[iband] = surfradin[iband]-(1-e[iband])*skyradin[iband];
            //wave5 = pow(wave[iband],5);
            //if(((c1*e[iband])/(PI*Re[iband]*wave5)+1) < 0) return(-1);
            //Te[iband] = (log((c1*e[iband])/(PI*Re[iband]*wave5)+1));
            Te[iband] = interp1d_npts(lut[band[iband]+1], lut[0], Re[iband]/e[iband], nlut_lines);
            // add error return if invalid Te in any band
            //if(fabs(Te[iband]) > 0) Te[iband] = (c2/wave[iband])/Te[iband]; 
            if(is_nan(Te[iband])) 
            {
                // With change to LUT, this cannot happen.
                WARNING("Te[%d] = nan",iband);
                return(1);
            }
            if(Te[iband] > *Tnem) *Tnem = Te[iband]; 
            diff[iband] = Re[iband]-Rold[iband];
#ifdef DEBUG_PIXEL
            if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
            {
                printf("        iband=%d Re=%.4f Te=%.4f diff=%.4f\n",
                        iband, Re[iband], Te[iband], diff[iband]);
            }
#endif
        }
        //  Compute the delta of corrected radiances - For convergence to occur,
        //  it must occur in all bands
        dcon = 0; ddiv = 0;
        for(iband=0; iband<n_channels; iband++) { 
            if(fabs(diff[iband]) < 0.05) dcon++;
            if(diff[iband] > 0.05) ddiv ++;
        }
#ifdef DEBUG_PIXEL
        if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
        {
            printf("        dcon=%d ddiv=%d\n", dcon, ddiv);
        }
#endif
        // must have more than 2 iterations 
        if(dcon == n_channels && k > 2) 
        {
#ifdef DEBUG_PIXEL
            if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
            {
                printf("        QA==0: Tnem=%.4f e=(%.4f",
                        *Tnem, e[0]);
                int idb;
                for (idb = 1; idb < n_channels; idb++)
                    printf(", %.4f", e[idb]);
                printf(")\n");
            }
#endif
            return(0);
        }

        if(ddiv == n_channels && k > 2)
        {
#ifdef DEBUG_PIXEL
            if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
            {
                printf("        QA==1: diverging results\n");
            }
#endif
            return(1); /* must have more than 2 iterations */
        }

        for(iband=0; iband<n_channels; iband++) { 
            //wave5 = pow(wave[iband],5);
            //B = c1/(wave5*PI*(pow(E,(c2/(wave[iband]*(*Tnem))))-1)); 
            B = interp1d_npts(lut[0], lut[band[iband]+1], *Tnem, nlut_lines);
            e[iband] = Re[iband]/B;
#ifdef DEBUG
            if(e[iband] < 0) {
                WARNING("EMIS < zero: e[%d]=%lf Re=%lf B=%lf Tnem = %lf",
                        iband, e[iband], Re[iband], B, *Tnem);
            }
#endif
            Rold[iband] = Re[iband];
        }
    }
#ifdef DEBUG_PIXEL
    if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
    {
        printf("        QA==1: End of loop without convergence\n");
    }
#endif
    return 1;
}

char *getRuntimeParameter(const char *name)
{
    ConfigElement *elem = findGroupElement(run_params, "RuntimeParameters", name);
    if (elem != NULL)
    {
        return getValue(elem, 0);
    }
    return NULL;
}

double getRuntimeParameter_f64(const char *name, double default_value)
{
    char *paramstr = getRuntimeParameter(name);
    if (paramstr != NULL)
        return atof(paramstr);
    return default_value;
}

int getRuntimeParameter_int(const char *name, int default_value)
{
    char *paramstr = getRuntimeParameter(name);
    if (paramstr != NULL)
        return atoi(paramstr);
    return default_value;
}

void apply_tes_algorithm(Matrix *Ts, Mat3d *emisf, MatUint16 *QC,
        Mat3d *surfradi, Mat3d *skyradi, MatUint8 *gp_water,
        Mat3d *t1r, Mat3d *Y)
{
    // Get size of the swath.
    int n_lines = surfradi->size2;
    int n_pixels = surfradi->size3;

    // Initialize the output arrays.
    mat_alloc(Ts, n_lines, n_pixels);
    mat3d_alloc(emisf, n_channels, n_lines, n_pixels);
    mat_uint16_alloc(QC, n_lines, n_pixels);

    int bmax;
    int lutindx;
    double surfradin[MAX_CHANNELS], skyradin[MAX_CHANNELS]; //, emis[MAX_CHANNELS];
    //double betamin, betamax, MMD;
    double Tf, Reff[MAX_CHANNELS];
    //double beta[MAX_CHANNELS], emin, bm;
    //double wave5;
    double bm2, beta2[MAX_CHANNELS], beta2min, beta2max, emis2[MAX_CHANNELS];
    double MMD2, emin2;
    int nbm2;

    int i;
    double emax;
    double *co;
    int16 qc_val;
    int nbare = 0;
    int nveg = 0;

    // Planck function outputs
    int kiter;
    double ef[MAX_CHANNELS];
    double Tnemf;
    int QA;
    bool ef_neg = false;

#ifdef DEBUG_TES
    long bad_emis_count = 0;
    long bad_temp_count = 0;
#endif

    double emis_val;
    for (iline = 0; iline < n_lines; iline++)
    {
        for(ipixel=0; ipixel < n_pixels; ipixel++)
        {
            qc_val = 0;
            // If all surfradi values are zero, set both bits in QC.
            bool non_zero = false;
            bool has_zero = false;
            for (i = 0; i < n_channels; i++)
            {
                if (fabs(mat3d_get(surfradi, i, iline, ipixel)) > epsilon)
                {
                    non_zero = true;
                }
                else
                {
                    has_zero = true;
                }
            }
            if (!non_zero) 
            {
                mat_uint16_set(QC, iline, ipixel, bits[0]|bits[1]);
                continue;
            }
            if(is_nan(mat3d_get(surfradi, 0, iline, ipixel))) {
                mat_uint16_set(QC, iline, ipixel, bits[1]);
                continue;
            }
            if (has_zero)
            {
                continue;
            }

            for(i=0; i<n_channels; i++) {
                surfradin[i] = mat3d_get(surfradi, i, iline, ipixel);
                skyradin[i] = mat3d_get(skyradi, i, iline, ipixel);
            }

            if(mat_uint8_get(gp_water, iline, ipixel) == 0) 
            {
                // Use VEG coefficients for non-gray.
                emax = emax_veg;
                co = co_veg;
                nveg++;
            }
            else 
            {
                // Use BARE coefficients for gray pixels.
                emax = emax_bare;
                co = co_bare;
                nbare++;
            }

            // Call Planck function
            QA = NEM_planck(emax, surfradin, skyradin, it,
                ef, &Tnemf, Reff, &kiter);

            if(QA == 1) {
                mat_uint16_set(QC, iline, ipixel,
                        qc_val | (bits[0]+bits[1])); 
#ifdef DEBUG_PIXEL
                if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
                {
                    printf("NEM_planck QA==1 emax=%.4f surfradin=%.4f skyradin=%.4f "
                           "it=%d ef=(%.4f, %.4f, %.4f) Tnemf=%.4f Reff=%.4f kiter=%d\n",
                           emax, surfradin[lst_band_index], skyradin[lst_band_index], it,
                           ef[0], ef[1], ef[2], Tnemf, *Reff, kiter);
                }
#endif
                continue;
            }
            ef_neg = false;
            for (i = 0; i < n_channels; i++)
            {
                if (ef[i] < 0.0)
                {
                    ef_neg = true;
                    break;
                }
            }
            if (ef_neg)
            {
                mat_uint16_set(QC, iline, ipixel, qc_val);
                continue;
            }

            //if(kiter >= 7) QC[ipixel] = QC[ipixel]; // note: no effect
            if(kiter == 6) qc_val |= bits[6]; 
            if(kiter == 5) qc_val |= bits[7]; 
            if(kiter < 5) qc_val |= (bits[6]+bits[7]);

            //betamax = nan1; betamin = -nan1;
            //bm = ef[0];
            //for(i=1; i<n_channels; i++) {
            //    if(ef[i] > bm) bm = ef[i];
            //}
            // Compute mean of bands 2,3,4
            beta2max = nan1; beta2min = -nan1;
            bm2 = 0;
            nbm2 = 0;
            for(i = band_index[BAND_2]; i <= band_index[BAND_4]; i++) {
                if(!is_nan(ef[i])) {
                    bm2 += ef[i];
                    nbm2++;
                }
            }
            if (nbm2 > 0) bm2 /= nbm2;

            for (i=0; i<n_channels; i++) 
            {
                //beta[i] = ef[i]/bm;

                // back to ef/avg(NEM value, but check for nan)
                if(nbm2 > 0) {
                    beta2[i] = ef[i]/bm2;       
                } else {
                    WARNING("invalid bm2 in finish module "
                          "ef = %lf %lf %lf %lf %lf nbm2 = %d",
                          ef[0],ef[1],ef[2],ef[3],ef[4],nbm2);
                }
                //if(beta[i] > betamax) betamax = beta[i];
                //if(beta[i] < betamin) betamin = beta[i];

                if(beta2[i] > beta2max) beta2max = beta2[i];
                if(beta2[i] < beta2min) beta2min = beta2[i];
            }
            //MMD = betamax-betamin;
            MMD2 = beta2max-beta2min;

            //emin = co[0]-co[1]*pow(MMD,co[2]);
            emin2 = co[0]-co[1]*pow(MMD2,co[2]);

            // Compute emissivity for each channel and for wide band.
            for(i=0; i<n_channels; i++)
            {
                //emis[i] = beta[i]*(emin/betamin);
                emis2[i] = beta2[i]*(emin2/beta2min);

                if(MMD2 < 0 || nbm2 == 0) 
                    emis_val = 0.0;
                else 
                {
                    emis_val = beta2[i]*(emin2/beta2min);
                }
                mat3d_set(emisf, i, iline, ipixel, emis_val);
#ifdef DEBUG_PIXEL
                if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
                {
                    printf("EMIS:   kiter=%d iband=%d [%d][%d]\n",
                            kiter, i, iline, ipixel);
                    printf("        emin2=%.4f beta2min=%.4f emisf=%.4f\n",
                            emin2, beta2min, emis_val);
                    printf("        beta2[%d]=%.4f beta2max=%.4f beta2min=%.4f\n",
                            i, beta2[i], beta2max, beta2min);
                }
#endif
#ifdef DEBUG_TES
                if (emis_val > 1.5 || emis_val < 0.0)
                {
                    bad_emis_count++;
                    if (bad_emis_count == 10)
                        WARNING("etc.", NONE);
                    else if (bad_emis_count < 10)
                        WARNING("emisf=%.4f is out of bounds at [%d][%d][%d]",
                            emis_val, i, iline, ipixel);
                }
#endif
            }
            // Use the clearest channel for LST calculation
            bmax = lst_band_index;
            lutindx = band[lst_band_index] + 1;
#ifdef DEBUG_PIXEL
            if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
            {
                printf("*DEBUG* at %d %d emis2[%d] = %.4f\n", 
                    iline, ipixel, bmax, emis2[bmax]);
            }
#endif
            double Reffc = Reff[bmax] / mat3d_get(emisf, bmax, iline, ipixel);
            assert(nlut_lines > 1);
            Tf = interp1d_npts(lut[lutindx], lut[0], Reffc, nlut_lines);
#ifdef DEBUG_PIXEL
            if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
            {
                printf("INTERP: bmax=%d lut[0]=%.1f...%.1f lut[%d]=%.1f...%.1f "
                        "nlut_lines=%d]\n",
                        bmax, lut[0][0], lut[0][nlut_lines-1], 
                        lutindx, lut[lutindx][0], lut[lutindx][nlut_lines-1],nlut_lines);
                printf("        Reffc=%.1f\n", Reffc);
                printf("LST:    Tf[%d][%d]=%.1f Reff[%d]=%.4f "
                        "emisf[%d]=%.4f Reffc=%.4f\n", 
                        iline, ipixel, Tf, bmax, Reff[bmax], bmax, 
                        mat3d_get(emisf, bmax, iline, ipixel), Reffc);
            }
#endif
            if(emis2[bmax] < 0.0) {
                mat_set(Ts, iline, ipixel, 0.0);
                for(i=0; i<n_channels; i++) 
                {
                    mat3d_set(emisf, i, iline, ipixel, 0.0);
                }
#ifdef DEBUG_PIXEL
                if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
                {
                    printf("Ts[%d][%d]=0.0 due to emis2[%d]=%.4f\n", 
                            iline, ipixel, bmax, emis2[bmax]);
                }
#endif
            } 
            else 
            {
                mat_set(Ts, iline, ipixel, Tf);
#ifdef DEBUG_TES
                if (Tf > 360.0 || 
                    (Tf > 0.0 && Tf < hard_t_min))
                {
                    bad_temp_count++;
                    if (bad_temp_count == 10)
                        WARNING("etc.", NONE);
                    else if (bad_temp_count < 10)
                        WARNING("Ts=%.4f is out of bounds at [%d][%d]",
                                 Tf, iline, ipixel);
                }
#endif
            }

            // Set quality flag based on band 4 and 5 emis data, 
            // and t1r for band 2.
            if( (mat3d_get(emisf, band_index[BAND_4], iline, ipixel) < 0.95 && 
                 mat3d_get(emisf, band_index[BAND_5], iline, ipixel) < 0.95)
                || mat3d_get(t1r, band_index[BAND_2], iline, ipixel) < 0.4) {
                    qc_val &= 0xFFFC;  // Clear bits 0, 1
                    qc_val |= bits[0]; // Set bit 0
            }
            // ECOSTRESS uses band 2 here
            double skyr1_over_y1 = 
                mat3d_get(skyradi, band_index[BAND_2], iline, ipixel) / 
                mat3d_get(Y, band_index[BAND_2], iline, ipixel);
            if (skyr1_over_y1 >= 0.3) 
            {
                // QC[ipixel] = QC[ipixel]; -- do nothing
            }
            else if(skyr1_over_y1 > 0.2) 
                qc_val |= bits[8];
            else if(skyr1_over_y1 > 0.1) 
                qc_val |= bits[9];
            else
                qc_val |= (bits[8]+bits[9]);
            if(MMD2 >= 0.15) 
            {
                // QC[ipixel] = QC[ipixel]; -- do nothing
            }
            else if(MMD2 > 0.1) qc_val |= bits[10];
            else if(MMD2 > 0.03) qc_val |= bits[11];
            else
                qc_val |= (bits[10]+bits[11]);
            mat_uint16_set(QC, iline, ipixel, qc_val);
#ifdef DEBUG_PIXEL
            if (iline == DEBUG_LINE && ipixel == DEBUG_PIXEL)
            {
                printf("QC[%d][%d]=%04X\n", 
                        iline, ipixel, qc_val);
            }
#endif
        }
    }
#ifdef DEBUG_TES
    if (bad_emis_count != 0)
        WARNING("Found emisf out of bounds %ld times.", bad_emis_count);
    if (bad_temp_count != 0)
        WARNING("Found Temp out of bounds %ld times.", bad_temp_count);
#endif
}
