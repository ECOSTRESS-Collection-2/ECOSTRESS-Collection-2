// See tes_util.h for credits and documentation.
// @copyright (c) Copyright 2016, Jet Propulsion Laboratories, Pasadena, CA
#include "error.h"
#include "fileio.h"
#include "lste_lib.h"
#include "tes_util.h"
#include <assert.h>
#include <stdio.h>
#include <hdf.h>
#include <grib_api.h>

#ifndef PATH_MAX
#define PATH_MAX 1024
#endif

#define COMMENT
//#define TESTING
//#define DEBUG
#ifdef DEBUG
#include "timer.h"
#endif
//#define DEBUG_NCEP
#ifdef DEBUG_NCEP
#define DEBUG_NCEP_LVL 24
#define DEBUG_NCEP_ROW 124
#define DEBUG_NCEP_COL 71
static bool debug_ncep_on = false;
#endif

static const double real_nan = REAL_NAN;

// for NCEP -- Pressure (hPa or mb)
#define MAX_NCEP_LVLS 100
static int npres = 0;
static double pres[MAX_NCEP_LVLS] = 
       {10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 
        300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 
        800.0, 850.0, 900.0, 925.0, 950.0, 975.0, 1000.0};

static int lastday[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static void check_leap_year(int year)
{
    static bool checked = false;
    if (checked) return;
    if (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0))
    {
        // leap year
        lastday[2] = 29;
    }
    checked = true;
}

// Remember NWP filename
char nwp_name[FILEIO_STR_SIZE] = "";
char *get_nwp_name()
{
    return nwp_name;
}

// [BF] Need to copy from 4D to 3D, 3D to 2D, and 1D to 1D.
// These functions pass m, a structure from matrix,
// as the target, and v, a pointer to an aray of double values.
void copy3d(Mat3d *m, double *v)
{
    int i;
    for (i = 0; i < m->size1 * m->size2 * m->size3; i++)
    {
        m->vals[i] = v[i];
    }
}

void copy2d(Matrix *m, double *v)
{
    int i;
    for (i = 0; i < m->size1 * m->size2; i++)
    {
        m->vals[i] = v[i];
    }
}

void copy1d(Vector *m, double *v)
{
    int i;
    for (i = 0; i < m->size; i++)
    {
        m->vals[i] = v[i];
    }
}

void constrict_cloud(MatUint8 *cloud, int n_neighbors)
{
    // TODO -- if n_neighbors > 2, corners will always be elliminated.
    int count;
    // Look at each pixel.
    int r, c;
    int nr, nc;
    int ix = 0;
    int nix;
    for (r = 0; r < cloud->size1; r++)
    {
        for (c = 0; c < cloud->size2; c++, ix++)
        {
            // No need to count neighbors when cloud is zero.
            if (cloud->vals[ix] == 0)
                continue;
            // count non-zero neighbors
            count = 0;
            for (nr = r-1; nr <= r+1; nr++)
            {
                // Skip when off the edges
                if (nr < 0 || nr >= cloud->size1)
                    continue;
                for (nc = c - 1; nc <= c + 1; nc++)
                {
                    // Skip when off the edges
                    if (nc < 0 || nc >= cloud->size2)
                        continue;
                    // Skip current location (self not a neighbor)
                    if (nr == r && nc == c)
                        continue;
                    // Get neighbor locatio.
                    nix = nr * cloud->size2 + nc;
                    // Count non-zero neighbor.
                    if (cloud->vals[nix] != 0)
                        count++;
                }
            }
            // If pixel does not have enough neighbors, change it to zero.
            if (count <= n_neighbors)
                cloud->vals[ix] = 0;
        }
    }
}

void extend_cloud(MatUint8 *cloud_in, MatUint8 *cextend, int cloud_extend, 
    uint8 extended_value)
{
    mat_uint8_clear(cextend);
    mat_uint8_alloc(cextend, cloud_in->size1, cloud_in->size2);
    // Extend Cloudy pixels in cloud_in into the cextend matrix
    int i, j, ii, jj;
    uint8 cext, cm_val, c_in;
    for (i = 0; i < cloud_in->size1; i++)
    {
        for (j = 0; j < cloud_in->size2; j++)
        {
            cm_val = mat_uint8_get(cloud_in, i, j);
            if (cm_val != 0)
            {
                // Found a Cloudy pixel. Extend it into the output.
                for (ii = i - cloud_extend; ii <= i + cloud_extend; ii++)
                {
                    if (ii < 0) continue;
                    if (ii >= cloud_in->size1) break;
                    for (jj = j - cloud_extend; jj <= j + cloud_extend; jj++)
                    {
                        if (jj < 0) continue;
                        if (jj >= cloud_in->size2) break;
                        c_in = mat_uint8_get(cloud_in, ii, jj);
                        // Keep original cloud setting
                        if (c_in != 0)
                            cext = c_in;
                        else
                            cext = extended_value;
                        mat_uint8_set(cextend, ii, jj, cext);
                    }
                }
            }
        }
    }
}

// --------------------------------------------------------------
// TesDate: Valid data ranges (constants)
// --------------------------------------------------------------

const int kTesDateYearRange[2] =    { 0, 3000 };
const int kTesDateMonthRange[2] =   { 1, 12 };
const int kTesDateDayRange[2] =     { 1, 31 };
const int kTesDateDoyRange[2] =     { 1, 366 };    // day of year
const int kTesDateHourRange[2] =    { 0, 23 };    // hour
const int kTesDateMinuteRange[2] =  { 0, 59 };    // minute

// --------------------------------------------------------------
// TesDate: check valid data ranges
// --------------------------------------------------------------

bool tes_check_date(const TesDate *date)
{
    assert(date != NULL);
    if(kTesDateYearRange[0]   <= date->year   && date->year   <=  kTesDateYearRange[1] &&
       kTesDateMonthRange[0]  <= date->month  && date->month  <=  kTesDateMonthRange[1] &&
       kTesDateDayRange[0]    <= date->day    && date->day    <=  kTesDateDayRange[1] &&
       kTesDateHourRange[0]   <= date->hour   && date->hour   <=  kTesDateHourRange[1] &&
       kTesDateMinuteRange[0] <= date->minute && date->minute <=  kTesDateMinuteRange[1])
    {
        return true;    // values are within specified range
    }
    else
    {
        return false;    // values are outside of valid range
    }
}

// --------------------------------------------------------------
// TesDate: Copy date, date = date_src
// --------------------------------------------------------------

void tes_copy_date(TesDate *date, const TesDate *date_src)
{
    assert(date != NULL && date_src != NULL);
    date->year = date_src->year;
    date->month = date_src->month;
    date->day = date_src->day;
    date->hour = date_src->hour;
    date->minute = date_src->minute;
}

// --------------------------------------------------------------
// TesDate: Set date, initializing values
// --------------------------------------------------------------

int tes_set_date(TesDate *date, int _year, int _month, int _day, 
                int _hour, int _minute)
{
    assert(date != NULL);
        
    // Set date values
    date->year = _year;
    date->month = _month;
    date->day = _day;
    date->hour = _hour;
    date->minute = _minute;
    if (!tes_check_date(date))
        return -1;
    return 0;
}

void tes_add_minutes(TesDate *date, int mins)
{
    date->minute += mins;
    if (date->minute >= 60)
    {
        date->hour += date->minute / 60;
        date->minute %= 60;
        if (date->hour > 23)
        {
            check_leap_year(date->year);
            date->day += date->hour / 24;
            date->hour %= 24;
            if (date->day > lastday[date->month])
            {
                date->day -= lastday[date->month];
                date->month++;
                if (date->month > 12)
                {
                    date->month -= 12;
                    date->year++;
                }
            }
        }
    }
}

// --------------------------------------------------------------
// TesDoyDate: Copy date, date = date_src
// --------------------------------------------------------------

void tes_copy_doy_date(TesDoyDate *date, const TesDoyDate *date_src)
{
    assert(date != NULL && date_src != NULL);
    date->year = date_src->year;
    date->doy = date_src->doy;
    date->hour = date_src->hour;
    date->minute = date_src->minute;
}

// --------------------------------------------------------------
// TesDate: Set date, initializing values
// --------------------------------------------------------------

int tes_set_doy_date(TesDoyDate *date, int year, int doy, int hour, int minute)
{
    assert(date != NULL);

    // Set date values
    date->year = year;
    date->doy = doy;
    date->hour = hour;
    date->minute = minute;
    // TBD check doy date
    return 0;
}

// --------------------------------------------------------------
// Matrix 2d double: Set values corresponding to uint8 mask == 1 to NaN
// --------------------------------------------------------------

void tes_mat_set_mask1_to_nan(Matrix *mat, const MatUint8 *mask)
{
    if (mat != NULL && mat->vals != NULL)
    {
        // Reset values to NaN
        size_t i = 0;
        size_t imax = mat->size1 * mat->size2;    // total size

        for(i = 0; i < imax; i++)
        {
            if(mask->vals[i] == 1) mat->vals[i] = REAL_NAN; // double NaN
        }
    }
}


// --------------------------------------------------------------
// Interpolate 3D matrix column by level to replace NaNs
// --------------------------------------------------------------
void interp_column(Vector *xvec, Mat3d *ymat, int zmax, int y, int x)
{
    // zmax is the number of levels
    double xvals[zmax];
    double yvals[zmax];
    int    nnans = 0; // Number of NaNs
    int    nnums = 0; // Number of non-NaNs.
    int    naniz[zmax];
    int    i, z;
    double yval;
    
    // Split the values up into indexes for NaNs and arrays of non-NaNs.
    for (z = 0; z < zmax; z++)
    {
        yval = mat3d_get(ymat, z, y, x);
        if (is_nan(yval))
        {
            naniz[nnans++] = z; // Save index to be interpolated
        }
        else
        {
            xvals[nnums] = vec_get(xvec, z);
            yvals[nnums] = yval;
            nnums++;
        }
    }
    
    // If no NaNs, nothing to do.
    if (nnans < 1) return;
    
    // If not min 2 non-NaNs, can't interpolate.
    assert(nnums > 1);
    
    // Replace NaNs with interpolated non-NaN values.
    for (i = 0; i < nnans; i++)
    {
        z = naniz[i];
        // Interpolate the replacement value 
        yval = interp1d_npts(xvals, yvals, vec_get(xvec, z), nnums);
        mat3d_set(ymat, z, y, x, yval);
    }
}

void init_rad(RAD *rad)
{
    mat_init(&rad->Lat);
    mat_init(&rad->Lon);
    mat_init(&rad->El);
    mat_init(&rad->Satzen);
    mat_init(&rad->Solzen);
    mat_uint8_init(&rad->Water);
    mat_init(&rad->Rad[0]);
    mat_init(&rad->Rad[1]);
    mat_init(&rad->Rad[2]);
    mat_init(&rad->Rad[3]);
    mat_init(&rad->Rad[4]);
    mat_init(&rad->SWIR);
    mat_int32_init(&rad->DataQ[0]);
    mat_int32_init(&rad->DataQ[1]);
    mat_int32_init(&rad->DataQ[2]);
    mat_int32_init(&rad->DataQ[3]);
    mat_int32_init(&rad->DataQ[4]);
    strcpy(rad->OrbitCorrectionPerformed, "unknown");
    rad->OrbitCorrectionPerformed[TES_STR_SIZE - 1] = '\0';
    strcpy(rad->FieldOfViewObstruction, "unknown");
    rad->FieldOfViewObstruction[TES_STR_SIZE - 1] = '\0';
}

void clear_rad(RAD *rad)
{
    mat_clear(&rad->Lat);
    mat_clear(&rad->Lon);
    mat_clear(&rad->El);
    mat_clear(&rad->Satzen);
    mat_clear(&rad->Solzen);
    mat_uint8_clear(&rad->Water);
    mat_clear(&rad->Rad[0]);
    mat_clear(&rad->Rad[1]);
    mat_clear(&rad->Rad[2]);
    mat_clear(&rad->Rad[3]);
    mat_clear(&rad->Rad[4]);
    mat_clear(&rad->SWIR);
    mat_int32_clear(&rad->DataQ[0]);
    mat_int32_clear(&rad->DataQ[1]);
    mat_int32_clear(&rad->DataQ[2]);
    mat_int32_clear(&rad->DataQ[3]);
    mat_int32_clear(&rad->DataQ[4]);
}

void tes_read_geo_data(RAD *rad, const char *GEO_filename)
{
    int stat;
    Matrix land_fraction;
    mat_init(&land_fraction);
    // Open the hdf file
    hid_t geo_id = open_hdf5(GEO_filename, H5F_ACC_RDONLY);
    if (geo_id == FAIL)
        ERREXIT(25, "Could not open file %s", GEO_filename);
    // Load the data
    stat = hdf5read_mat(geo_id, "Geolocation/latitude", &rad->Lat);
    if (stat < 0)
        ERREXIT(26, "Unable to read Latitude from %s", GEO_filename);
    stat = hdf5read_mat(geo_id, "Geolocation/longitude", &rad->Lon);
    if (stat < 0)
        ERREXIT(26, "Unable to read Longitude from %s", GEO_filename);
    stat = hdf5read_mat(geo_id, "Geolocation/height", &rad->El);
    if (stat < 0)
        ERREXIT(26, "Unable to read Height from %s", GEO_filename);
    stat = hdf5read_mat(geo_id, "Geolocation/view_zenith", &rad->Satzen);
    if (stat < 0)
        ERREXIT(26, "Unable to read view_zenith from %s", GEO_filename);
    stat = hdf5read_mat(geo_id, "Geolocation/solar_zenith", &rad->Solzen);
    if (stat < 0)
        ERREXIT(26, "Unable to read solar_zenith from %s", GEO_filename);
    stat = hdf5read_mat(geo_id, "Geolocation/land_fraction", &land_fraction);
    if (stat < 0)
        ERREXIT(26, "Unable to read land_fraction from %s", GEO_filename);

    // Get the product metadata "/L1GEOMetadata/OrbitCorrectionPerformed"
    hid_t gid, dsetid, dtypeid,strtype;
    char *valchar[1]; // HDF5 returns a point to a string, and this data should be a string.
    gid = H5Gopen2(geo_id, "L1GEOMetadata", H5P_DEFAULT);
    if (gid < 0)
        WARNING("Unable to open group L1GEOMetadata from %s", GEO_filename);
    else
    {
        dsetid = H5Dopen2(gid, "OrbitCorrectionPerformed", H5P_DEFAULT);
        if (dsetid < 0)
            WARNING("Unable to open /L1GEOMetadata/OrbitCorrectionPerformed "
                    "from %s", GEO_filename);
        else
        {
            dtypeid = H5Dget_type(dsetid);
            strtype = H5Tcopy(dtypeid);
            stat = H5Dread(dsetid, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, valchar);
            if (stat < 0)
                WARNING("Unable to read /L1GEOMetadata/OrbitCorrectionPerformed "
                        "from %s", GEO_filename);
            else
                strncpy(rad->OrbitCorrectionPerformed, valchar[0], TES_STR_SIZE - 1);
            H5Dclose(dsetid);
        }
        H5Gclose(gid);
    }

    // Get the product metadata "/StandardMetadata/FieldOfViewObstruction"
    gid = H5Gopen2(geo_id, "StandardMetadata", H5P_DEFAULT);
    if (gid < 0)
        WARNING("Unable to open group StandardMetadata from %s", GEO_filename);
    else
    {
        dsetid = H5Dopen2(gid, "FieldOfViewObstruction", H5P_DEFAULT);
        if (dsetid < 0)
            WARNING("Unable to open /StandardMetadata/FieldOfViewObstruction "
                    "from %s", GEO_filename);
        else
        {
            dtypeid = H5Dget_type(dsetid);
            strtype = H5Tcopy(dtypeid);
            stat = H5Dread(dsetid, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, valchar);
            if (stat < 0)
                WARNING("Unable to read /StandardMetadata/FieldOfViewObstruction "
                        "from %s", GEO_filename);
            else
                strncpy(rad->FieldOfViewObstruction, valchar[0], TES_STR_SIZE - 1);
            H5Dclose(dsetid);
        }
    H5Gclose(gid);
    }

    close_hdf5(geo_id);

    // Create Water (a.k.a., oceanpix) from land_fraction
    mat_uint8_alloc(&rad->Water, land_fraction.size1, land_fraction.size2);

    // Perform scaling
    int i;
    const int DSIZE = rad->El.size1 * rad->El.size2;
    for (i = 0; i < DSIZE; i++)
    {
        rad->El.vals[i] /= 1000.0; // convert Height to km
        /* Remove scaling for zenith angles. Jan 2018
        rad->Satzen.vals[i] *= 0.01;
        rad->Solzen.vals[i] *= 0.01;
        */
        if (land_fraction.vals[i] < 50.0) // If less than 50% land,
            rad->Water.vals[i] = 1;       // then mark it as water.
    }
}

int get_rad_nbands(const char *RAD_filename)
{
    int i;
    int stat;
    int nbands = 0;
    Vector bsa;
    vec_init(&bsa);

    hid_t h5fid = open_hdf5(RAD_filename, H5F_ACC_RDONLY);
    if (h5fid == FAIL)
        ERREXIT(27, "Could not open file %s", RAD_filename);
    // Read BandSpecification metadata
    stat = hdf5read_vec(h5fid, "/L1B_RADMetadata/BandSpecification", &bsa);
    if (stat >= 0)
    {
        // Count bands that are present. Skip [0] (SWIR band).
        for (i = 1; i < 6; i++)
        {
            if (bsa.vals[i] > 0.1)
                nbands++;
        }
        INFO("BandSpecification indicates %d bands for L2.", nbands);
    }
    else
    {
        nbands = 5;
        INFO("Missing BandSpecification indicates 5 bands for L2.", NONE);
    }
    close_hdf5(h5fid);
    vec_clear(&bsa);
    return nbands;
}

void tes_read_rad_data(RAD *rad, const char *RAD_filename,
         int n_bands, int *band_numbers)
{
    int b;
    int band_no;
    char name_buf[128];
    int stat;
    // Open the hdf file
    hid_t rad_id = open_hdf5(RAD_filename, H5F_ACC_RDONLY);
    if (rad_id == FAIL)
        ERREXIT(27, "Could not open file %s", RAD_filename);
    // Load the data
    for (b = 0; b < n_bands; b++)
    {
        band_no = band_numbers[b] + 1;
        sprintf(name_buf, "/Radiance/radiance_%d", band_no);
        stat = hdf5read_mat(rad_id, name_buf, &rad->Rad[b]);
        if (stat < 0)
            ERREXIT(28, "Unable to read %s from %s", 
                    name_buf, RAD_filename);
        // Replace radiance fill values -9999, -9998, and -9997 with NaN
        // [1.0.4] All negative RAD values should be NAN.
        int ir;
        for (ir = 0; ir < rad->Rad[b].size1 * rad->Rad[b].size2; ir++)
        {
            if (rad->Rad[b].vals[ir] < -9990.0)
                rad->Rad[b].vals[ir] = REAL_NAN;
        }
        sprintf(name_buf, "/Radiance/data_quality_%d", band_no);
        stat = hdf5read_mat_int32(rad_id, name_buf, &rad->DataQ[b]);
        if (stat < 0)
            ERREXIT(28, "Unable to read %s from %s", 
                    name_buf, RAD_filename);
    }
    rad->nchannels = n_bands;
/* [BF] Removed SWIR band processing
    stat = hdf5read_mat(rad_id, "/SWIR/swir_dn", &rad->SWIR);
    if (stat < 0)
        ERREXIT(28, "Unable to read swir_dn from %s", RAD_filename);
    hdf5_scale_and_offset_mat(rad_id, "SWIR/swir_dn",  &rad->SWIR);
*/
    close_hdf5(rad_id);
}

// --------------------------------------------------------------
// init - constructor
// --------------------------------------------------------------

void tes_atm_init(TesATM *atm)
{
    vec_init(&(atm->lat));
    vec_init(&(atm->lon));
    vec_init(&(atm->lev));
    mat3d_init(&atm->t);
    mat3d_init(&atm->q);
    mat_init(&atm->sp);
    mat_init(&atm->skt);
    mat_init(&atm->t2);
    mat_init(&atm->q2);
    mat_init(&atm->tcw);
}

// --------------------------------------------------------------
// ATM: Clear, deallocate -- Destructor
// --------------------------------------------------------------

void tes_atm_clear(TesATM *atm)
{
    assert(atm != NULL);
    
    vec_clear(&(atm->lat));
    vec_clear(&(atm->lon));
    vec_clear(&(atm->lev));
    mat3d_clear(&(atm->t));
    mat3d_clear(&(atm->q));
    mat_clear(&(atm->sp));
    mat_clear(&(atm->skt));
    mat_clear(&(atm->t2));
    mat_clear(&(atm->q2));
    mat_clear(&(atm->tcw));
}

void crop_atm(TesATM *outATM, TesATM *inATM, CropSize *crop)
{
    tes_atm_init(outATM); // Init output struct.

    // Init index settings
    int iLatStart = 0;
    int iLatEnd = inATM->lat.size;
    int iLonStart = 0;
    int iLonEnd = inATM->lon.size;
    // Determine limited Lat range
    while (iLatStart < iLatEnd
        && (inATM->lat.vals[iLatStart] < crop->minLat 
            || inATM->lat.vals[iLatStart] > crop->maxLat))
        iLatStart++;
    while (iLatEnd > iLatStart
        && (inATM->lat.vals[iLatEnd - 1] < crop->minLat 
            || inATM->lat.vals[iLatEnd - 1] > crop->maxLat))
        iLatEnd--;
    // Determine limited Lon range
    while (iLonStart < iLonEnd
        && (inATM->lon.vals[iLonStart] < crop->minLon 
            || inATM->lon.vals[iLonStart] > crop->maxLon))
        iLonStart++;
    while (iLonEnd > iLonStart
        && (inATM->lon.vals[iLonEnd - 1] < crop->minLon 
            || inATM->lon.vals[iLonEnd - 1] > crop->maxLon))
        iLonEnd--;
    // Get the sizes.
    int iLatSize = iLatEnd - iLatStart;
    int iLonSize = iLonEnd - iLonStart;
    assert(iLatSize > 0 && iLonSize > 0);
    // Allocate the output struct
    vec_alloc(&(outATM->lev), inATM->lev.size);
    vec_alloc(&(outATM->lat), iLatSize);
    vec_alloc(&(outATM->lon), iLonSize);
    mat3d_alloc(&(outATM->t), inATM->lev.size, iLatSize, iLonSize);
    mat3d_alloc(&(outATM->q), inATM->lev.size, iLatSize, iLonSize);
    mat_alloc(&(outATM->sp), iLatSize, iLonSize);
    // Allocate optional parts.
    if (inATM->skt.size1 > 0)
        mat_alloc(&(outATM->skt), iLatSize, iLonSize);
    if (inATM->t2.size1 > 0)
        mat_alloc(&(outATM->t2), iLatSize, iLonSize);
    if (inATM->q2.size1 > 0)
        mat_alloc(&(outATM->q2), iLatSize, iLonSize);
    if (inATM->tcw.size1 > 0)
        mat_alloc(&(outATM->tcw), iLatSize, iLonSize);
    // Copy cropped part of data
    int ilat, ilon, ilev;
    for (ilat = 0; ilat < iLatSize; ilat++)
    {
        for (ilon = 0; ilon < iLonSize; ilon++)
        {
            for (ilev = 0; ilev < inATM->lev.size; ilev++)
            {
                mat3d_set(&(outATM->t), ilev, ilat, ilon,
                    mat3d_get(&(inATM->t), ilev, ilat+iLatStart, ilon+iLonStart));
                mat3d_set(&(outATM->q), ilev, ilat, ilon,
                    mat3d_get(&(inATM->q), ilev, ilat+iLatStart, ilon+iLonStart));
            }
            mat_set(&(outATM->sp), ilat, ilon,
                mat_get(&(inATM->sp), ilat+iLatStart, ilon+iLonStart));
            if (inATM->skt.size1 > 0)
                mat_set(&(outATM->skt), ilat, ilon,
                    mat_get(&(inATM->skt), ilat+iLatStart, ilon+iLonStart));
            if (inATM->t2.size1 > 0)
                mat_set(&(outATM->t2), ilat, ilon,
                    mat_get(&(inATM->t2), ilat+iLatStart, ilon+iLonStart));
            if (inATM->q2.size1 > 0)
                mat_set(&(outATM->q2), ilat, ilon,
                    mat_get(&(inATM->q2), ilat+iLatStart, ilon+iLonStart));
            if (inATM->tcw.size1 > 0)
                mat_set(&(outATM->tcw), ilat, ilon,
                    mat_get(&(inATM->tcw), ilat+iLatStart, ilon+iLonStart));
        }
    }

    // Copy vectors
    for (ilat = 0; ilat < iLatSize; ilat++)
    {
        outATM->lat.vals[ilat] = inATM->lat.vals[ilat+iLatStart];
    }
    for (ilon = 0; ilon < iLonSize; ilon++)
    {
        outATM->lon.vals[ilon] = inATM->lon.vals[ilon+iLonStart];
    }
    for (ilev = 0; ilev < inATM->lev.size; ilev++)
    {
        outATM->lev.vals[ilev] = inATM->lev.vals[ilev];
    }

#ifdef DEBUG
    printf("crop_atm: from %lux%lu to %dx%d from %d,%d to %d,%d\n",
        inATM->lat.size, inATM->lon.size, iLatSize, iLonSize, 
        iLatStart, iLonStart, iLatEnd, iLonEnd);
#endif
}

void interp_atmos(TesATM *atm0, TesATM *atm1, TesATM *atm, int hh, int mm, int period_hrs)
{
    // Interpolate from base hour: 00, 06, 12, or 18
    // to the current hour + minute/60.
    double dhh = (hh % period_hrs) + mm / 60.0;
    double ri = dhh / (double)period_hrs;
    
    // Allocate atm output structure
    mat3d_alloc(&atm->t, atm0->t.size1, atm0->t.size2, atm0->t.size3);
    mat3d_alloc(&atm->q, atm0->q.size1, atm0->q.size2, atm0->q.size3);
    mat_alloc(&atm->sp, atm0->sp.size1, atm0->sp.size2);
    vec_alloc(&atm->lat, atm0->lat.size);
    vec_alloc(&atm->lon, atm0->lon.size);
    vec_alloc(&atm->lev, atm0->lev.size);
    // Some datasets are optional.
    if (atm0->skt.size1 > 0)
        mat_alloc(&atm->skt, atm0->skt.size1, atm0->skt.size2);
    if (atm0->t2.size1 > 0)
        mat_alloc(&atm->t2, atm0->t2.size1, atm0->t2.size2);
    if (atm0->q2.size1 > 0)
        mat_alloc(&atm->q2, atm0->q2.size1, atm0->q2.size2);
    if (atm0->tcw.size1 > 0)
        mat_alloc(&atm->tcw, atm0->tcw.size1, atm0->tcw.size2);

    int i, j, k;
    int ioff, joff;
    // Interpolate 3D data
    for (i = 0; i < atm->t.size1; i++)
    {
        ioff = i * atm->t.size2 * atm->t.size3;
        for (j = 0; j < atm->t.size2; j++)
        {
            joff = j * atm->t.size3;
            for (k = 0; k < atm->t.size3; k++)
            {
                atm->t.vals[ioff + joff + k] = 
                    atm0->t.vals[ioff + joff + k] +
                    ri * 
                    (atm1->t.vals[ioff + joff + k] - 
                     atm0->t.vals[ioff + joff + k]);
                atm->q.vals[ioff + joff + k] = 
                    atm0->q.vals[ioff + joff + k] +
                    ri * 
                    (atm1->q.vals[ioff + joff + k] - 
                     atm0->q.vals[ioff + joff + k]);
            }
        }
    }
    // Interpolate 2D data
    for (i = 0; i < atm->sp.size1; i++)
    {
        ioff = i * atm->sp.size2;
        for (j = 0; j < atm->sp.size2; j++)
        {
            atm->sp.vals[ioff + j] =
                atm0->sp.vals[ioff + j] +
                ri * 
                (atm1->sp.vals[ioff + j] - 
                 atm0->sp.vals[ioff + j]);
            // The following are not always supplied:
            if (atm->skt.size1 > 0)
                atm->skt.vals[ioff + j] =
                    atm0->skt.vals[ioff + j] +
                    ri * 
                    (atm1->skt.vals[ioff + j] - 
                     atm0->skt.vals[ioff + j]);
            if (atm->t2.size1 > 0)
                atm->t2.vals[ioff + j] =
                    atm0->t2.vals[ioff + j] +
                    ri * 
                    (atm1->t2.vals[ioff + j] - 
                     atm0->t2.vals[ioff + j]);
            if (atm->q2.size1 > 0)
                atm->q2.vals[ioff + j] =
                    atm0->q2.vals[ioff + j] +
                    ri * 
                    (atm1->q2.vals[ioff + j] - 
                     atm0->q2.vals[ioff + j]);
            if (atm->tcw.size1 > 0)
                atm->tcw.vals[ioff + j] =
                    atm0->tcw.vals[ioff + j] +
                    ri * 
                    (atm1->tcw.vals[ioff + j] - 
                     atm0->tcw.vals[ioff + j]);
        }
    }
    // Non-interpolated items
    copy1d(&atm->lat, atm0->lat.vals);
    copy1d(&atm->lon, atm0->lon.vals);
    copy1d(&atm->lev, atm0->lev.vals);
}

// Read one data file in NetCDF format.
void read_atmos_netcdf(
        const char* netcdffilename,
        Vector *height, Vector *lat, Vector *lon, 
        Mat4d *t, Mat4d *q, Mat3d *sp
        )
{    
    int ncid;
    ncid = open_netcdf(netcdffilename, NC_NOWRITE);
    if (ncid == -1)
        ERREXIT(31, "Unable to open input file %s", netcdffilename);

    int stat;
    // Read datasets
    stat = netcdfread_vec(ncid, "lev", height);
    if (stat < 0) ERREXIT(31, "Could not read dataset 'lev' from %s", netcdffilename);
    stat = netcdfread_vec(ncid, "lat", lat);
    if (stat < 0) ERREXIT(31, "Could not read dataset 'lat' from %s", netcdffilename);
    stat = netcdfread_vec(ncid, "lon", lon);
    if (stat < 0) ERREXIT(31, "Could not read dataset 'lon' from %s", netcdffilename);
    
    stat = netcdfread_mat4d(ncid, "T", t);
    if (stat < 0) ERREXIT(31, "Could not read dataset 'T' from %s", netcdffilename);
    stat = netcdfread_mat4d(ncid, "QV", q);
    if (stat < 0) ERREXIT(31, "Could not read dataset 'QV' from %s", netcdffilename);
    stat = netcdfread_mat3d(ncid, "PS", sp);
    if (stat < 0) ERREXIT(31, "Could not read dataset 'PS' from %s", netcdffilename);
    
    // Close NetCDF file
    close_netcdf(ncid);
}

// ========================================================
// GEOS Processing
// ========================================================
/// @brief Converts a year, month, day, hour to a GEOS filename
///
/// @param[in]  year        Year (4 digits)
/// @param[in]  month       Month (1 to 12)
/// @param[in]  day         Day of month (1 to 31)
/// @param[in]  hour        Hour of day (0 to 23)
/// @param[in]  geos_dir    Merra collection top directory
/// @param[out] geosname    Name of file [PATH_MAX]
/// @return 	0=Success   -1=Failure
int geos_filename(int year, int month, int day, int hour, 
        const char *geos_dir, char *geosname)
{
    // Each file is for a three hour period. Mod 3 to get hour of file.
    int h3 = hour - (hour % 3);
    // Execute a linux find command to look for the GEOS file.
    char searchcmd[PATH_MAX];
    sprintf(searchcmd, 
            "find -L %s -name 'GEOS*3d_asm_Np*%4d%02d%02d_%02d00*.nc4'", 
            geos_dir, year, month, day, h3);
    FILE *pipe = popen(searchcmd, "r");
    if (pipe == NULL)
        ERREXIT(29, "Could not execute command: %s", searchcmd);
    char *res = fgets(geosname, PATH_MAX-1, pipe);
    pclose(pipe);
    if (res == NULL)
    {
        WARNING("Unable to find GEOS file using:\n    %s", searchcmd);
        return -1;
    }
    // Remove newline & other junk from end of line
    int ends = strlen(geosname) - 1;
    while (ends > 0 && geosname[ends] < ' ')
        geosname[ends--] = '\0';

#ifdef COMMENT
    printf("Found GEOS file: %s\n", geosname);
#endif

    // Indicate success.
    return 0;
}

void tes_read_interp_atmos_geos(TesATM *atm, const TesDate *tdate, 
        const char *geos_dir, CropSize *crop)
{
    // Prepare the atm structs for before and after.
    TesATM atm0, atm1;
    tes_atm_init(&atm0);
    tes_atm_init(&atm1);
    Mat4d T4d, Q4d;
    mat4d_init(&T4d);
    mat4d_init(&Q4d);
    Mat3d SP3d;
    mat3d_init(&SP3d);

    // Get the date/time that is mod three hours before tdate.
    TesDate before;
    tes_copy_date(&before, tdate);
    before.hour -= before.hour % 3;
    before.minute = 0;

    // Get the date/time three hours later.
    TesDate after;
    tes_copy_date(&after, &before);
    tes_add_minutes(&after, 3*60);

    // Find the GEOS before file.
    char geos_name1[TES_STR_SIZE];
    if (geos_filename(before.year, before.month, before.day, before.hour, 
                      geos_dir, geos_name1) < 0)
        ERREXIT(31, "Could not find GEOS file 1 under %s", geos_dir);
    read_atmos_netcdf(geos_name1, &atm0.lev, &atm0.lat, &atm0.lon, 
        &T4d, &Q4d, &SP3d);
    // Convert 4D data with lead dimension of 1 to 3D
    atm0.t.size1 = T4d.size2;
    atm0.t.size2 = T4d.size3;
    atm0.t.size3 = T4d.size4;
    atm0.t.vals = T4d.vals;
    atm0.q.size1 = Q4d.size2;
    atm0.q.size2 = Q4d.size3;
    atm0.q.size3 = Q4d.size4;
    atm0.q.vals = Q4d.vals;
    // Likewise, 3D to 2D.
    atm0.sp.size1 = SP3d.size2;
    atm0.sp.size2 = SP3d.size3;
    atm0.sp.vals = SP3d.vals;

    strcpy(nwp_name, geos_name1);
    strcat(nwp_name, ", ");

    // Find the GEOS after file.
    char geos_name2[TES_STR_SIZE];
    if (geos_filename(after.year, after.month, after.day, after.hour, 
                      geos_dir, geos_name2) < 0)
        ERREXIT(31, "Could not find GEOS file 2 under %s", geos_dir);
    strcat(nwp_name, geos_name2);
    // Reuse temporary matrixes.
    mat4d_init(&T4d);
    mat4d_init(&Q4d);
    mat3d_init(&SP3d);

    read_atmos_netcdf(geos_name2, &atm1.lev, &atm1.lat, &atm1.lon, 
        &T4d, &Q4d, &SP3d);
    atm1.t.size1 = T4d.size2;
    atm1.t.size2 = T4d.size3;
    atm1.t.size3 = T4d.size4;
    atm1.t.vals = T4d.vals;
    atm1.q.size1 = Q4d.size2;
    atm1.q.size2 = Q4d.size3;
    atm1.q.size3 = Q4d.size4;
    atm1.q.vals = Q4d.vals;
    atm1.sp.size1 = SP3d.size2;
    atm1.sp.size2 = SP3d.size3;
    atm1.sp.vals = SP3d.vals;

    // Remove NaNs from 3D datasets
    int x, y;
    for (y = 0; y < atm0.t.size2; y++)
    {
        for (x = 0; x < atm0.t.size3; x++)
        {
            interp_column(&atm0.lev, &atm0.t, atm0.t.size1, y, x);
            interp_column(&atm0.lev, &atm0.q, atm0.q.size1, y, x);
            interp_column(&atm1.lev, &atm1.t, atm1.t.size1, y, x);
            interp_column(&atm1.lev, &atm1.q, atm1.q.size1, y, x);
        }
    }

    // Interpolate between the before and after data.
    TesATM crop_atm0;
    tes_atm_init(&crop_atm0);
    TesATM crop_atm1;
    tes_atm_init(&crop_atm1);
    // Crop the atms.
    crop_atm(&crop_atm0, &atm0, crop);
    crop_atm(&crop_atm1, &atm1, crop);
    int interp_hh = tdate->hour % 3; // Each file is 3 hours.
    interp_atmos(&crop_atm0, &crop_atm1, atm, interp_hh, tdate->minute, 3);

#ifdef TESTING
    TesATM interp_atm;
    tes_atm_init(&interp_atm);
    interp_atmos(&atm0, &atm1, &interp_atm, interp_hh, tdate->minute, 3);
    mat3d_flipdim(&interp_atm.t);
    mat3d_flipdim(&interp_atm.q);
    hid_t dbg_id = create_hdf5("GEOS_debug.h5");
    if (dbg_id == FAIL)
        ERROR("Could not create %s", "GEOS_debug.h5");
    herr_t dbg_stat  = hdf5write_mat3d(dbg_id, "NWP_T", &interp_atm.t);
    dbg_stat += hdf5write_mat3d(dbg_id, "NWP_QV", &interp_atm.q);
    dbg_stat += hdf5write_mat(dbg_id, "NWP_SP", &interp_atm.sp);
    close_hdf5(dbg_id);
    tes_atm_clear(&interp_atm);
#endif

    // Tidy up temp structs
    tes_atm_clear(&atm0);
    tes_atm_clear(&atm1);
    tes_atm_clear(&crop_atm0);
    tes_atm_clear(&crop_atm1);

    // Adjust data
    const double Rair = 287.0;
    const double Rwater = 461.5;
    const double wv_coeff = 1.0 / (1e-6 * (Rair/Rwater));

    int i;    
    for (i = 0; 
         i < atm->q.size1 * atm->q.size2 * atm->q.size3; 
         i++)
    {
        if (!is_nan(atm->q.vals[i]))
        {
            atm->q.vals[i] *= wv_coeff;
            if (atm->q.vals[i] > 1e5) 
                atm->q.vals[i] = real_nan;
            else if (atm->q.vals[i] < 0.000001)
                atm->q.vals[i] = 0.000001;
        }
        if (!is_nan(atm->t.vals[i]))
        {
            if (atm->t.vals[i] < hard_t_min)
                atm->t.vals[i] = hard_t_min;
            else if (atm->t.vals[i] > hard_t_max)
                atm->t.vals[i] = hard_t_max;
        }
    }

    for (i = 0; i < atm->sp.size1 * atm->sp.size2; i++)
    {
        if (!is_nan(atm->sp.vals[i]))
        {
            atm->sp.vals[i] /= 100.0;
            if (atm->sp.vals[i] < hard_p_min)
                atm->sp.vals[i] = hard_p_min;
            else if (atm->sp.vals[i] > hard_p_max)
                atm->sp.vals[i] = hard_p_max;
        }
    }

    // Need to flip dims for T and QV.
    mat3d_flipdim(&atm->t);
    mat3d_flipdim(&atm->q);
}

// ========================================================
// MERRA Processing
// ========================================================
/// @brief Converts a doy (day-of-year) date into a MERRA filename
///
/// Part of the problem is that Merra data uses various version numbers,
/// like MERRA2_300_ or MERRA2_400_. This function creates a matching
/// string, then searches for the actual filename on disk.
///
/// @param[in]      year        Year (4 digits)
/// @param[in]      doy         Day of year (1 to 366)
/// @param[in]      merradir    Merra collection top directory
/// @param[out]     merraname   Name of file [PATH_MAX]
/// @return 	0=Success   -1=Failure
//
//  MERRA filenames are in the form:
//    MERRA2_300.inst6_3d_ana_Np.yyyymmdd.nc4
//  where yyyy is the year, mm the month, and dd the day, e.g., 
//    MERRA2_300.inst6_3d_ana_Np.20040803.nc4

int merra_doy_filename(int year, int doy, const char *merradir, char *merraname)
{
    assert(year > 1900);
    assert(doy > 0);
    assert(doy <= 366);
    check_leap_year(year);
    int mm = 1;
    int dd = doy;
    while (dd > lastday[mm])
    {
        dd -= lastday[mm];
        mm++;
    }

    char searchcmd[PATH_MAX];
    sprintf(searchcmd, "find -L %s -name '*MERRA2*inst6*Np*%4d%02d%02d.nc4'", 
            merradir, year, mm, dd);
    FILE *pipe = popen(searchcmd, "r");
    if (pipe == NULL)
        ERREXIT(29, "Could not execute command: %s", searchcmd);
    char *res = fgets(merraname, PATH_MAX-1, pipe);
    pclose(pipe);
    if (res == NULL)
    {
        WARNING("Unable to get output of: %s", searchcmd);
        return -1;
    }
    // Remove newline & other junk from end of line
    int ends = strlen(merraname) - 1;
    while (ends > 0 && merraname[ends] < ' ')
        merraname[ends--] = '\0';

    char *check_merra = strstr(merraname, "MERRA2");
    if (check_merra == NULL)
    {
        WARNING("Could not find MERRA filename in %s", merraname);
        return -1;
    }
    return 0;
}

/// @brief Creates the MERRA filename for the date
///
/// Part of the problem is that Merra data uses various version numbers,
/// like MERRA2_300_ or MERRA2_400_. This function creates a match
/// string, then searches for the actual filename on disk.
///
/// @param[in]      year        Year (4 digits)
/// @param[in]      month       Month (1 to 12)
/// @param[in]      day         Day of month (1 to 31)
/// @param[in]      merradir    Merra collection top directory
/// @param[out]     merraname   Name of file [PATH_MAX]
/// @return 	0=Success   -1=Failure
//
//  MERRA filenames are in the form:
//    MERRA2_300.inst6_3d_ana_Np.yyyymmdd.nc4
//  where yyyy is the year, mm the month, and dd the day, e.g., 
//    MERRA2_300.inst6_3d_ana_Np.20040803.nc4

int merra_filename(int year, int month, int day, 
        const char *merradir, char *merraname)
{
    assert(year > 1900);
    assert(month >= 1);
    assert(month <= 12);
    assert(day > 0);
    assert(day <= 31);
    
    char searchcmd[PATH_MAX];
    sprintf(searchcmd, "find -L %s -name '*MERRA2*inst6*Np*%4d%02d%02d*.nc4'", 
            merradir, year, month, day);
    FILE *pipe = popen(searchcmd, "r");
    if (pipe == NULL)
        ERREXIT(29, "Could not execute command: %s", searchcmd);
    char *res = fgets(merraname, PATH_MAX-1, pipe);
    pclose(pipe);
    if (res == NULL)
    {
        WARNING("Unable to get output of: %s", searchcmd);
        return -1;
    }
    // Remove newline & other junk from end of line
    int ends = strlen(merraname) - 1;
    while (ends > 0 && merraname[ends] < ' ')
        merraname[ends--] = '\0';

    char *check_merra = strstr(merraname, "MERRA2");
    if (check_merra == NULL)
    {
        WARNING("Failed to get a valid MERRA2 file from the command: %s", 
                searchcmd);
        WARNING("Command returned output: %s", merraname);
        return -1;
    }
    return 0;
}

/* -- HDF4 files are not currently supported.

// Read one MERRA data file in HDF4 format.
void read_atmos_hdf(
        const char* hdffilename,
        Vector *height, Vector *lat, Vector *lon, 
        Mat4d *t, Mat4d *q, Mat3d *sp
        )
{    
    int stat = 0;

    // Open HDF file in read-only mode
    // Selecting merra1 or merra2 input file
    hdf4id_t sd_fid = 0;    // CHECK: change default to HDF FAIL or -1
    
    sd_fid = SDstart(hdffilename, DFACC_RDONLY ); // open merra1 file
  
    // Check that file is open
    if(sd_fid == FAIL)
        ERROR("SDstart cannot open HDF file %s", hdffilename);
 
    // Read HDF 1d array datasets
    // Pressure = double(hdfread(file_MERRA,'Height'))
    stat = hdfread_vec(sd_fid, "Height", height);
    if (stat < 0)
        ERROR("Could not read Height from %s", hdffilename);
    // Lat_atmos = double(hdfread(file_MERRA,'YDim'))
    stat = hdfread_vec(sd_fid, "YDim", lat);
    if (stat < 0)
        ERROR("Could not read YDim from %s", hdffilename);
    // Lon_atmos = double(hdfread(file_MERRA,'XDim'))
    stat = hdfread_vec(sd_fid, "XDim", lon);
    if (stat < 0)
        ERROR("Could not read XDim from %s", hdffilename);
    
    // Read HDF 4d float32 array datasets
    // Temperature = double(hdfread(file_MERRA,'T'))
    stat = hdfread_mat4d(sd_fid, "T", t);
    if (stat < 0)
        ERROR("Could not read T from %s", hdffilename);
    // H2O_MR = double(hdfread(file_MERRA,'QV'))
    stat = hdfread_mat4d(sd_fid, "QV", q);
    if (stat < 0)
        ERROR("Could not read QV from %s", hdffilename);
    
    // Read HDF 3d float32 array dataset
    // Surface_Pressure = double(hdfread(file_MERRA,'PS'))
    stat = hdfread_mat3d(sd_fid, "PS", sp);
    if (stat < 0)
        ERROR("Could not read PS from %s", hdffilename);
    
    // Close HDF file
    SDend(sd_fid);
}
*/

// --------------------------------------------------------------
// Read ATM data for two time periods and interpolate between them.
// --------------------------------------------------------------
    
// Read HDF (or NetCDF) file and datasets
// 1. Corresponding to MATLAB function: read_atmos_merra
// 2. Assuming that input ATM matrices are empty
// 3. merra_filename exists
// 4. ihour = 0 for 0000, 1 for 0600, 2 for 1200, or 3 for 1800.
// If ihour == 3, the second Merra filename will be used to read
// the next day, and the second dataset will use ihour 0. 
//

void tes_read_interp_atmos_merra(TesATM *atm, const TesDate *tdate, 
     const char *merra_dir)
{
    assert(atm != NULL);
    assert(tdate != NULL);

    int i, x, y;
    int stat;
    // Get hour and minute of data from tdate data.
    int hh = tdate->hour;
    int mm = tdate->minute;
    int ihour = hh / 6; // Get the index for the quarter of the day

    // Prepare the atm structs for before and after.
    TesATM atm0, atm1;
    tes_atm_init(&atm0);
    tes_atm_init(&atm1);

    // Vectors and matrices read from HDF arrays
    Vector Pressure;        vec_init(&Pressure);
    Vector Lat;             vec_init(&Lat);
    Vector Lon;             vec_init(&Lon);
    Mat4d Temperature;      mat4d_init(&Temperature);
    Mat4d H2O_MR;           mat4d_init(&H2O_MR);
    Mat3d Surface_Pressure; mat3d_init(&Surface_Pressure);
    
    // Get Merra filename
    // Get hour and minute of data from tdate data.
    char pathname[FILEIO_STR_SIZE];
    stat = merra_filename(tdate->year, tdate->month, tdate->day, 
            merra_dir, pathname);
    if (stat < 0)
    {
        ERREXIT(30, "MERRA file not found -- %s", pathname);
    }
#ifdef COMMENT
    INFO("MERRA2 file 1: %s", pathname);
#endif

    char *merra_root = strstr(pathname, "MERRA");
    if (merra_root == NULL)
        merra_root = pathname;
    strcpy(nwp_name, merra_root);

    read_atmos_netcdf(pathname, 
        &Pressure, &Lat, &Lon, &Temperature, &H2O_MR, &Surface_Pressure);

    // Scale and perform clamps, conversions, etc.
    // Convert water vapor:
    // % Q_MR to Q_PPMV
    // %                                                   MW(Dry Air)
    // %           ppmv(MOL) = 1000 . Mixing_Ratio(MOL) . -------------
    // %                                                   MW(MOL)
    // %       where MW(Dry Air) = Average molecular weight of dry air (28.966).
    // %             MW(MOL)     = Molecular weight of the gas in question (H2O =
    // %             18.015).
    // % q = q*1000; % kg/kg to g/kg
    // % q = 1000*q*(28.966/18.015); % g/kg to PPMV
    // % 
    // % or
    // % wv_mmr = 1.e-6 * wv_ppmv_layer * (Rair / Rwater)
    // % wv_mmr in kg/kg, Rair = 287.0, Rwater = 461.5
    const double Rair = 287.0;
    const double Rwater = 461.5;
    const double wv_coeff = 1.0 / (1e-6 * (Rair/Rwater));
    
    // 4D
    for (i = 0; 
         i < H2O_MR.size1 * H2O_MR.size2 * H2O_MR.size3 * H2O_MR.size4; 
         i++)
    {
        if (!is_nan(H2O_MR.vals[i]))
        {
            H2O_MR.vals[i] *= wv_coeff;
            if (H2O_MR.vals[i] > 1e5) 
                H2O_MR.vals[i] = real_nan;
            else if (H2O_MR.vals[i] < 0.000001)
                H2O_MR.vals[i] = 0.000001;
        }
        if (!is_nan(Temperature.vals[i]) && Temperature.vals[i] < 150.0)
        {
#ifdef DEBUG
            fprintf(stderr, "T[%d] clamped to 150.0\n", i);
#endif
            Temperature.vals[i] = 150.0;
        }
    }
    // 3D
    for (i = 0; i < Surface_Pressure.size1 * Surface_Pressure.size2 
         * Surface_Pressure.size3; i++)
    {
        if (!is_nan(Surface_Pressure.vals[i]))
            Surface_Pressure.vals[i] /= 100.0;
    }
    
    // Copy data to the first ATM 
    int i3doff = ihour * Surface_Pressure.size2 * Surface_Pressure.size3;
    int i4doff = ihour * Temperature.size2 * Temperature.size3 
        * Temperature.size4;
    mat3d_alloc(&atm0.t, Temperature.size2, Temperature.size3,
        Temperature.size4);
    copy3d(&atm0.t, &Temperature.vals[i4doff]);
    mat3d_alloc(&atm0.q, H2O_MR.size2, H2O_MR.size3,
        H2O_MR.size4);
    copy3d(&atm0.q, &H2O_MR.vals[i4doff]);
    mat_alloc(&atm0.sp, Surface_Pressure.size2, Surface_Pressure.size3);
    copy2d(&atm0.sp, &Surface_Pressure.vals[i3doff]);
    vec_alloc(&atm0.lat, Lat.size);
    copy1d(&atm0.lat, Lat.vals);
    vec_alloc(&atm0.lon, Lon.size);
    copy1d(&atm0.lon, Lon.vals);
    vec_alloc(&atm0.lev, Pressure.size);
    copy1d(&atm0.lev, Pressure.vals);

    ihour = (ihour + 1) % 4;
    if (ihour == 0)
    {
	    check_leap_year(tdate->year);
        // Next quarter day is in the next day. Read the next MERRA file.
        int nextday_year = tdate->year;
        int nextday_month = tdate->month;
        int nextday_day = tdate->day + 1;
        if (nextday_day > lastday[tdate->month])
        {
            // Is this the last day of year?
            if (tdate->month == 12)
            {
                // Happy New Year!
                nextday_year++;
                nextday_month = 1;
            }
            else
            {
                nextday_month++;
            }
            nextday_day = 1;
        }
        stat = merra_filename(nextday_year, nextday_month, nextday_day, 
                merra_dir, pathname);
        if (stat < 0)
        {
            ERREXIT(30, "Could not find second Merra file for %d.%d.%d", 
                nextday_year, nextday_month, nextday_day);
        }
#ifdef COMMENT
        INFO("MERRA2 file 2: %s", pathname);
#endif
    }

    read_atmos_netcdf(pathname, 
        &Pressure, &Lat, &Lon, &Temperature, &H2O_MR, &Surface_Pressure);

    // Repeat scale and clamping for second input set
    // TBD make a subroutine that reads the dataset and put all this 
    // duplicate code into it so, if there is a problem, it gets fixed in
    // all places when all places = one place. Also, the clamps and
    // conversions can be applied to the end data that is 1/4th of
    // the 3D and 4D datasets.
    for (i = 0; 
         i < H2O_MR.size1 * H2O_MR.size2 * H2O_MR.size3 * H2O_MR.size4; 
         i++)
    {
        if (!is_nan(H2O_MR.vals[i]))
        {
            H2O_MR.vals[i] *= wv_coeff;
            if (H2O_MR.vals[i] > 1e5) 
                H2O_MR.vals[i] = REAL_NAN;
            else if (H2O_MR.vals[i] < 0.000001)
                H2O_MR.vals[i] = 0.000001;
        }

        if (!is_nan(Temperature.vals[i]) && Temperature.vals[i] < 150.0)
        {
#ifdef DEBUG
            fprintf(stderr, "T2 clamped to 150 at %d\n", i);
#endif
            Temperature.vals[i] = 150.0;
        }
    }

    for (i = 0; i < Surface_Pressure.size1 * Surface_Pressure.size2 
         * Surface_Pressure.size3; i++)
    {
        if (!is_nan(Surface_Pressure.vals[i]))
            Surface_Pressure.vals[i] /= 100.0;
    }

    // Copy data to the second ATM 
    i3doff = ihour * Surface_Pressure.size2 * Surface_Pressure.size3;
    i4doff = ihour * Temperature.size2 * Temperature.size3 
        * Temperature.size4;
    mat3d_alloc(&atm1.t, Temperature.size2, Temperature.size3,
        Temperature.size4);
    copy3d(&atm1.t, &Temperature.vals[i4doff]);
    mat3d_alloc(&atm1.q, H2O_MR.size2, H2O_MR.size3,
        H2O_MR.size4);
    copy3d(&atm1.q, &H2O_MR.vals[i4doff]);
    mat_alloc(&atm1.sp, Surface_Pressure.size2, Surface_Pressure.size3);
    copy2d(&atm1.sp, &Surface_Pressure.vals[i3doff]);
    vec_alloc(&atm1.lat, Lat.size);
    copy1d(&atm1.lat, Lat.vals);
    vec_alloc(&atm1.lon, Lon.size);
    copy1d(&atm1.lon, Lon.vals);
    vec_alloc(&atm1.lev, Pressure.size);
    copy1d(&atm1.lev, Pressure.vals);

#ifdef DEBUG
    time_t startinterp = getusecs();
#endif 

    // Remove NaNs from 3D datasets
    for (y = 0; y < atm0.t.size2; y++)
    {
        for (x = 0; x < atm0.t.size3; x++)
        {
            interp_column(&Pressure, &atm0.t, atm0.t.size1, y, x);
            interp_column(&Pressure, &atm0.q, atm0.q.size1, y, x);
            interp_column(&Pressure, &atm1.t, atm1.t.size1, y, x);
            interp_column(&Pressure, &atm1.q, atm1.q.size1, y, x);
        }
    }

#ifdef DEBUG
    time_t elapsedinterp = elapsed(startinterp);
    fprintf(stderr, "column interpolations took %.1f secs.\n", 
            toseconds(elapsedinterp));
#endif
    
    // Interpolate between datasets to the given hour
    interp_atmos(&atm0, &atm1, atm, hh, mm, 6);

    // Need to flip dims for T and QV.
    mat3d_flipdim(&atm->t);
    mat3d_flipdim(&atm->q);
}


// ========================================================
// NCEP Processing
// ========================================================
//calculate_saturation_vapor_pressure_liquid
//    case 'Murphy&Koop2005'
//        if ~all(123 < T & T < 332)
//            warning('calculate_saturation_vapor_pressure_liquid:outRange', ...
//                'Temperature out of range [123-332] K.');
//        end
//        temp = 54.842763 - 6763.22 ./ T - 4.210 .* log(T) + 0.000367 .* T ...
//          + tanh( 0.0415 * (T - 218.8) ) ...
//          .* (53.878 - 1331.22 ./ T - 9.44523 .* log(T) + 0.014025 .* T);
//        es = exp(temp);
//        % D. M. MURPHY and T. KOOP
//        % Review of the vapour pressures of ice and supercooled water for
//        % atmospheric applications
//        % Q. J. R. Meteorol. Soc. (2005), 131, pp. 1539-1565 
//        % doi: 10.1256/qj.04.94
//        % <http://dx.doi.org/10.1256/qj.04.94>
//        % 
//        % "Widely used expressions for water vapour (Goff and Gratch 1946; Hyland and
// Wexler 1983) are being applied outside the ranges of data used by the original authors
// for their fits. This work may be the first time that data on the molar heat capacity of
// supercooled water [i.e., at temperatures below its freezing temperature] have been used
// to constrain its vapour pressure.
static bool have_warned = false;
double calculate_saturation_vapor_pressure_liquid(double T)
{
    if (T <= 123.0 || T >= 332.0)
    {
        if (!have_warned)
        {
            WARNING("calculate_saturation_vapor_pressure_liquid -- "
                    "temperature %.1f out of range [123-332] K", T);
            have_warned = true;
        }
    }

    double temp = 54.842763 - 6763.22 / T - 4.210 * log(T) + 0.000367 * T 
        + tanh( 0.0415 * (T - 218.8) )
        * (53.878 - 1331.22 / T - 9.44523 * log(T) + 0.014025 * T);
#ifdef DEBUG_NCEP
    if (debug_ncep_on)
    {
        printf("            temp=%g exp(temp)=%g\n",
                temp, exp(temp));
    }
#endif
    return exp(temp);
}

// Convert relative humidity to mixing level using temperatture (K) and 
// pressure (Pa).
//    case 'relative humidity'
//        RH = in;
//        es = calculate_saturation_vapor_pressure_liquid(T, method_es);
//        e = RH./100 .* es;
//    case 'mixing ratio'
//        Pd = P - e;
//        r = (e./Pd) .* c;
//        out = r;
double convert_relative_humidity(double p, double t, double rh)
{
    const double M_wet = 18.0152;  // kg / kmol,  molar mass of water vapor
    const double M_dry = 28.9644;  // kg / kmol,  molar mass of dry air
    const double c = M_wet / M_dry;

    double es = calculate_saturation_vapor_pressure_liquid(t);
    double e  = rh / 100.0 * es;
    double pd = p - e;
#ifdef DEBUG_NCEP
    if (debug_ncep_on)
    {
        printf("            p=%g t=%g rh=%g c=%g\n",
            p, t, rh, c);
        printf("            es=%g e=%g pd=%g ppwv=%g\n",
            es, e, pd, (e / pd) * c);
    }
#endif
    return (e / pd) * c;
}

/// @brief Determines the NCEP filename for the date.
///
/// @param[in]      year        Year (4 digits)
/// @param[in]      month       Month (1 to 12)
/// @param[in]      day         Day of month (1 to 31)
/// @param[in]      hour        Hour: one of 00, 06, 12, or 18
/// @param[in]      ncepdir     NCEP collection top directory
/// @param[out]     ncepname    Name of file [PATH_MAX]
/// @return 	0=Success   -1=Failure
//
//  NCEP filenames are in the form:
//    gdas1.PGrbF00.yymmdd.hhz
//  where yy is the year, mm the month, dd the day, and hh the hour, e.g., 
//    gdas1.PGrbF00.120811.18z

int ncep_filename(int year, int month, int day, int hour,
        const char *ncepdir, char *ncepname)
{
    year = year % 100;
    assert(month >= 1);
    assert(month <= 12);
    assert(day > 0);
    assert(day <= 31);
    assert(hour == 0 || hour == 6 || hour == 12 || hour == 18);

    char searchcmd[PATH_MAX];
    sprintf(searchcmd, "find -L %s -name 'gdas1*%02d%02d%02d.%02dz'", 
            ncepdir, year, month, day, hour);
    FILE *pipe = popen(searchcmd, "r");
    if (pipe == NULL)
        ERREXIT(29, "Could not execute command: %s", searchcmd);
    char *res = fgets(ncepname, PATH_MAX-1, pipe);
    pclose(pipe);
    if (res == NULL)
    {
        WARNING("Unable to get output of: %s", searchcmd);
        return -1;
    }
    // Remove newline & other junk from end of line
    int ends = strlen(ncepname) - 1;
    while (ends > 0 && ncepname[ends] < ' ')
        ncepname[ends--] = '\0';

    char *check_ncep = strstr(ncepname, "gdas1");
    if (check_ncep == NULL)
    {
        WARNING("Failed to get a valid NCEP filename from the command: %s", 
                searchcmd);
        WARNING("Command returned output: %s", ncepname);
        return -1;
    }
    return 0;
}

// NCEP grib1 grids are ordered by Latitude and Longitude:
//     Latitude  (rows) 90 ... 0, -1 ... -90
//     Longitude (cols) 0 ... 179, -180 ... -1. 
//
// RTTOV needs the order to be:
//     Latitude -90 ... -1, 0 ... 90
//     -180 ... -1, 0 ... 179. 
// To get this, the columns of the first and last half of each row
// are swapped. The rows are inverted.
// 
void copy_and_swap_lon_halfs(
        double *pncep, double *prttov, int maxlat, int maxlon)
{
    int inlat = maxlat - 1;
    int inlon = maxlon / 2; // start at middle of lon
    int outlat = 0;
    int outlon = 0;
    int count = 0;
    int inpos, outpos;
    while (count < maxlat * maxlon)
    {
        inpos = inlat * maxlon + inlon;
        outpos = outlat * maxlon + outlon;
        prttov[outpos] = pncep[inpos];

        inlon++;
        if (inlon >= maxlon) // wrap around
            inlon = 0;
        else if (inlon == maxlon / 2) // End at middle
            inlat--;

        outlon++;
        if (outlon >= maxlon)
        {
            outlat++;
            outlon = 0;
        }

        count++;
    }
}

// Find or add Pressure level. If overflow, return -1;
// Else, return pressure level index.
int get_pressure_level(double p)
{
    int i;
    for (i = 0; i < npres; i++)
    {
        if (dequal(p, pres[i]))
            return i;
    }
    // New pressure level
    if (npres >= MAX_NCEP_LVLS)
    {
        WARNING("Exceeded max of %d NCEP pressure levels at P=%.1f",
            MAX_NCEP_LVLS, p);
            return -1;
    }
    i = npres;
    pres[npres++] = p;
    return i;
}

// This function helps read_atmos_ncep initialize matrixes.
void init_ncep(
    grib_handle *h, size_t *pNi, size_t *pNj, TesATM *atm, 
        Matrix *gvals, size_t *nvals)
{
    // Get coordinates.
    GRIB_CHECK(grib_get_long(h, "Ni", (long*)pNi), 0);
    GRIB_CHECK(grib_get_long(h, "Nj", (long*)pNj), 0);
    mat_alloc(gvals, *pNi, *pNj);
    mat3d_alloc(&atm->t, MAX_NCEP_LVLS, *pNj, *pNi);
    mat3d_alloc(&atm->q, MAX_NCEP_LVLS, *pNj, *pNi);
    mat3d_fill(&atm->q, 0.000001);
    *nvals = *pNi * *pNj;
}

// This function reads a grib NCEP file and interprets the fields
// into the TesATM struct: Time, lat, lon, lev, p, t, q, sp.
// It also produces TCW from the pwat data.
void read_atmos_ncep(const char* ncepfilename, TesATM *atm)
{
    int err;
    size_t Ni, Nj;
    size_t nvals;
    Matrix grib_values; mat_init(&grib_values);
    long level;
    bool init = false;
    size_t slen = 1024;
    char shortName[1024];
    char typeOfLevel[1024];
    int tlvl = 0;
    int qlvl = 0;
    int desti;

    FILE* gribf = fopen(ncepfilename, "r");
    if (gribf == NULL)
        ERREXIT(35, "Could not open grib file: %s", ncepfilename);
    grib_handle* h = NULL;
    while((h = grib_handle_new_from_file(0, gribf, &err)) != NULL)
    {
        // Found a message. Get the short name.
        GRIB_CHECK(grib_get_length(h, "shortName", &slen), 0);
        grib_get_string(h, "shortName", shortName, &slen);
        if (strcmp(shortName, "t") == 0)
        {
            // Has t matrix been allocated?
            if (!init)
            {
                init_ncep(h, &Ni, &Nj, atm, &grib_values, &nvals);
                init = true;
            }

            // Found Temperature. Want level type isobaricInhPa.
            GRIB_CHECK(grib_get_length(h, "typeOfLevel", &slen), 0);
            grib_get_string(h, "typeOfLevel", typeOfLevel, &slen);
            if (strcmp(typeOfLevel, "isobaricInhPa") == 0)
            {
                // Get the pressure level.
                GRIB_CHECK(grib_get_long(h, "level", &level), 0);
                tlvl = get_pressure_level((double)level);
                if (tlvl >= 0)
                {
                    // Index to start of current t level
                    desti = tlvl * nvals;
                    // Get grib level ata.
                    GRIB_CHECK(grib_get_double_array(
                        h, "values", grib_values.vals, &nvals),0);
                    copy_and_swap_lon_halfs(
                        grib_values.vals, &atm->t.vals[desti], Nj, Ni);
                }
            }
            else if (strcmp(typeOfLevel, "surface") == 0)
            {
                // This is the skin temperature.
                mat_alloc(&atm->skt, Nj, Ni);
                GRIB_CHECK(grib_get_double_array(
                    h, "values", grib_values.vals, &nvals),0);
                copy_and_swap_lon_halfs(
                        grib_values.vals, atm->skt.vals, Nj, Ni);
            }
        }
        else if (strcmp(shortName, "r") == 0)
        {
            // Found Relative Humidity. Want level type isobaricInhPa.
            GRIB_CHECK(grib_get_length(h, "typeOfLevel", &slen), 0);
            grib_get_string(h, "typeOfLevel", typeOfLevel, &slen);
            if (strcmp(typeOfLevel, "isobaricInhPa") == 0)
            {
                // Want this Humidity. Has matrix been allocated?
                if (!init)
                {
                    init_ncep(h, &Ni, &Nj, atm, &grib_values, &nvals);
                    init = true;
                }

                // Get the pressure level.
                GRIB_CHECK(grib_get_long(h, "level", &level), 0);
                qlvl = get_pressure_level((double)level);
                if (qlvl >= 0)
                {
                    // Index to start of current t level
                    desti = qlvl * nvals;
                    // Get grib level ata.
                    GRIB_CHECK(grib_get_double_array(
                        h, "values", grib_values.vals, &nvals),0);
                    copy_and_swap_lon_halfs(
                            grib_values.vals, &atm->q.vals[desti], Nj, Ni);
                }
            }
        }
        else if (strcmp(shortName, "sp") == 0)
        {
            // Found Surface pressure
            if (!init)
            {
                init_ncep(h, &Ni, &Nj, atm, &grib_values, &nvals);
                init = true;
            }
            mat_alloc(&atm->sp, Nj, Ni);

            // Get the surface pressure.
            GRIB_CHECK(grib_get_double_array(
                h, "values", grib_values.vals, &nvals),0);
            copy_and_swap_lon_halfs(
                    grib_values.vals, atm->sp.vals, Nj, Ni);
        }
        else if (strcmp(shortName, "2t") == 0)
        {
            // Found 2 meter temperature
            if (!init)
            {
                init_ncep(h, &Ni, &Nj, atm, &grib_values, &nvals);
                init = true;
            }
            mat_alloc(&atm->t2, Nj, Ni);

            // Get the surface pressure.
            GRIB_CHECK(grib_get_double_array(
                h, "values", grib_values.vals, &nvals),0);
            copy_and_swap_lon_halfs(
                    grib_values.vals, atm->t2.vals, Nj, Ni);
        }
        else if (strcmp(shortName, "pwat") == 0)
        {
            // Found Precipitable water
            if (!init)
            {
                init_ncep(h, &Ni, &Nj, atm, &grib_values, &nvals);
                init = true;
            }
            mat_alloc(&atm->tcw, Nj, Ni);

            // Get the surface pressure.
            GRIB_CHECK(grib_get_double_array(
                h, "values", grib_values.vals, &nvals),0);
            copy_and_swap_lon_halfs(
                    grib_values.vals, atm->tcw.vals, Nj, Ni);
        }

        // Release the handle to the current message.
        grib_handle_delete(h);
    }
    fclose(gribf);
    mat_clear(&grib_values);

    // "trim" 3D arrays
    atm->t.size1 = npres;
    atm->q.size1 = npres;

    // Convert values:
    int x, y, z;
    double p, q;
    mat_alloc(&atm->q2, Nj, Ni);
    for (y = 0; y < Nj; y++)
    {
        for (x = 0; x < Ni; x++)
        {
            // 3D
            for (z = 0; z < npres; z++)
            {
                q = mat3d_get(&atm->q, z, y, x);
#ifdef DEBUG_NCEP
                if (z == DEBUG_NCEP_LVL &&
                    y == DEBUG_NCEP_ROW &&
                    x == DEBUG_NCEP_COL) 
                {
                    debug_ncep_on = true;
                    printf("NCEP_DEBUG: z=%d y=%d x=%d Q=%.1f T=%.1f\n",
                        z, y, x, q, mat3d_get(&atm->t, z, y, x));
                }
#endif
                q = convert_relative_humidity(
                    pres[z]* 100.0, mat3d_get(&atm->t, z, y, x), q);
                if (q <= 0.000001) q = 0.000001; //% avoid zero
                mat3d_set(&atm->q, z, y, x, q/(1e-6*(287.0/461.5)));
                //mat3d_set(&atm->q, z, y, x, q);
#ifdef DEBUG_NCEP
                if (debug_ncep_on)
                {
                    printf("            Q=%g\n", 
                            q/(1e-6*(287.0/461.5)));
                }
                debug_ncep_on = false;
#endif
            }

            // Find the pressure that is above the surface pressure.
            p = mat_get(&atm->sp, y, x);
            for (z = npres - 1; z > 0; z--)
            {
                if (pres[z] < p)
                {
                    // Found pressure level less than surface.
                    break;
                }
            }

            mat_set(&atm->q2, y, x, mat3d_get(&atm->q, z, y, x));

            // sp = sp/100; % Pa to hPa
            mat_set(&atm->sp, y, x, mat_get(&atm->sp, y, x)/100.0);
            if (atm->tcw.size1 > 0)
                mat_set(&atm->tcw, y, x, mat_get(&atm->tcw, y, x)/10.0);
        }
    }

    // Finally, fill out the static data for vectors.
    vec_alloc(&atm->lat, Nj);
    double latval = -90.0;
    for (y = 0; y < Nj; y++)
    {
        atm->lat.vals[y] = latval;
        latval += 1.0;
    }

    vec_alloc(&atm->lon, Ni);
    double lonval = -180.0;
    for (x = 0; x < Ni; x++)
    {
        atm->lon.vals[x] = lonval;
        lonval += 1.0;
    }

    vec_alloc(&atm->lev, npres);
    for (z = 0; z < npres; z++)
        atm->lev.vals[z] = pres[z];
}

void tes_read_interp_atmos_ncep(TesATM *atm, 
        const TesDate *tdate, const char *ncep_dir)
{
    int stat;

    // Prepare the atm structs and tcw matrix for before and after.
    TesATM atm1, atm2;
    tes_atm_init(&atm1);
    tes_atm_init(&atm2);

    // Determine the first filename for NCEP data.
    int hh = tdate->hour;
    int mm = tdate->minute;
    int ihour = hh / 6; // Get the index for the quarter of the day
    int ncephour = ihour * 6;
    int ncepday = tdate->day;
    int ncepmonth = tdate->month;
    int ncepyear = tdate->year;
    check_leap_year(ncepyear);

    char pathname[FILEIO_STR_SIZE];
    stat = ncep_filename(ncepyear, ncepmonth, ncepday, ncephour, 
            ncep_dir, pathname);
    if (stat < 0)
    {
        ERREXIT(37, "NCEP file 1 not found -- %s", pathname);
    }
    strcpy(nwp_name, pathname);

    #ifdef COMMENT
    INFO("NCEP file 1: %s", pathname);
#endif

    // Input the first NCEP file into atm1, tcw1.
    read_atmos_ncep(pathname, &atm1);

    // Determine the second NCEP filename.
    ncephour += 6;
    if (ncephour > 18)
    {
        ncephour = 0;
        ncepday++;
        if (ncepday > lastday[ncepmonth])
        {
            ncepday = 1;
            ncepmonth++;
            if (ncepmonth > 12)
            {
                ncepmonth = 1;
                ncepyear++;
            }
        }
    }
    stat = ncep_filename(ncepyear, ncepmonth, ncepday, ncephour, 
            ncep_dir, pathname);
    if (stat < 0)
    {
        ERREXIT(37, "NCEP file 2 not found -- %s", pathname);
    }
#ifdef COMMENT
    INFO("NCEP file 2: %s", pathname);
#endif

    // Input the second NCEP file into atm2, tcw2.
    read_atmos_ncep(pathname, &atm2);

    // Interpolate between datasets to the given hour/minute
    interp_atmos(&atm1, &atm2, atm, hh, mm, 6);

    // Free local structures.
    tes_atm_clear(&atm1);
    tes_atm_clear(&atm2);

#ifdef DEBUG_NCEP
    hid_t dbgf = create_hdf5("output/DebugNCEP_interp.h5");
    stat =  hdf5write_vec(dbgf, "P", &atm->lev);
    stat += hdf5write_vec(dbgf, "Lat", &atm->lat);
    stat += hdf5write_vec(dbgf, "Lon", &atm->lon);
    stat += hdf5write_mat3d(dbgf, "T", &atm->t);
    stat += hdf5write_mat3d(dbgf, "Q", &atm->q);
    stat += hdf5write_mat(dbgf, "SP", &atm->sp);
    if (atm->skt.size1 > 0)
       stat += hdf5write_mat(dbgf, "SKT", &atm->skt);
    if (atm->t2.size1 > 0)
       stat += hdf5write_mat(dbgf, "T2", &atm->t2);
    if (atm->q2.size1 > 0)
       stat += hdf5write_mat(dbgf, "Q2", &atm->q2);
    if (atm->tcw.size1 > 0)
       stat += hdf5write_mat(dbgf, "TCW", &atm->tcw);
    close_hdf5(dbgf);
    if (stat < 0)
        WARNING("%d errors writing datasets to DebugNCEP_interp.h5", -stat);
#endif
}
