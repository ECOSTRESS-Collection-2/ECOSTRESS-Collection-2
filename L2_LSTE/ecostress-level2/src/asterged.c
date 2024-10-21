// asterged.c
// See asterged.h for credits and documentation

#include <assert.h>
#include <dirent.h>
#include <math.h>
#include <mfhdf.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "asterged.h"
#include "bool.h"
#include "fileio.h"
#include "interps.h"
#include "lste_lib.h"
#include "maptogrid.h"

//#define DEBUG_ASTER
//#define DEBUG_EMIS_NANS
//#define ENABLE_TIMING
#ifdef ENABLE_TIMING
#include "timer.h"
#endif

// Parameters
//#define USE_INTERP_ID
#ifdef USE_INTERP_ID
double r2 = 30.0;
#endif


void ag_init(AsterGed *ag)
{
    ag->sensor = 0;
    ag->startLat = 0;
    ag->startLon = 0;
    ag->endLat = 0;
    ag->endLon = 0;
    ag->nLats = 0;
    ag->nLons = 0;
    mat_init(&ag->emis);
    mat_init(&ag->ndvi);
    ag->nFiles = 0;
    ag->nDefault = 0;
    ag->verbose = 0;
}

void ag_clear(AsterGed *ag)
{
    mat_clear(&ag->emis);
    mat_clear(&ag->ndvi);
}

double ag_get_emis(AsterGed *ag, double lat, double lon)
{
    size_t ilat = (size_t)((lat - ag->startLat) * 100);
    size_t ilon = (size_t)((lon - ag->startLon) * 100);
    if (ag->verbose > 8)
        printf("get(%.2f, %.2f) gets cell %lu, %lu\n",
                lat, lon, ilat, ilon);
    return mat_get(&ag->emis, ilat, ilon);
}

double ag_get_ndvi(AsterGed *ag, double lat, double lon)
{
    if (!ag->read_ndvi)
        return ASTER_NAN;
    size_t ilat = (size_t)((lat - ag->startLat) * 100);
    size_t ilon = (size_t)((lon - ag->startLon) * 100);
    return mat_get(&ag->ndvi, ilat, ilon);
}

int ag_load(AsterGed *ag, const char *asterdir, double minLat, double minLon,
            double maxLat, double maxLon)
{
    const size_t Band11 = 1;
    const size_t Band12 = 2;
    const size_t Band13 = 3;
    const size_t Band14 = 4;

    char aster_filename[FILEIO_MAX_FILENAME];
    if (ag->verbose > 0)
    {
        printf("ag_load: \n    dir=%s minLat=%.2f minLon=%.2f "
        "maxLat=%.2f maxLon=%.2f\n",
        asterdir, minLat, minLon, maxLat, maxLon);
    }

    ag->startLat = floor(minLat);
    ag->startLon = floor(minLon);
    ag->endLat = ceil(maxLat);
    ag->endLon = ceil(maxLon);
    ag->nLats = (int)((ag->endLat - ag->startLat) * 100);
    ag->nLons = (int)((ag->endLon - ag->startLon) * 100);
    if (ag->verbose > 1)
    {
        printf("    Start Lat/Lon = %.2f,%.2f\n"
               "    End Lat/Lon   = %.2f,%.2f\n"
               "    nLats=%lu nLons=%lu\n",
                ag->startLat, ag->startLon,
                ag->endLat, ag->endLon, ag->nLats, ag->nLons);
    }
    // Allocate the arrays
    mat_alloc(&ag->emis, ag->nLats, ag->nLons);
    if (ag->read_ndvi)
    {
        mat_alloc(&ag->ndvi, ag->nLats, ag->nLons);
    }

    // Loop to fill in the array 100 x 100 from each file that 
    // covers the area.

    // Have to read the file one degree higher
    int latChunk = (int)ag->startLat + 1;
    int lonChunk = (int)ag->startLon;
    int iLat = 0;
    int iLon = 0;
    char latnum[18], lonnum[8];
    Mat3dInt16 echunk;
    mat3d_int16_init(&echunk);
    MatInt16 nchunk;
    mat_int16_init(&nchunk);
    int x, y, i, j;
    double ndvival = 0.0;
    double e11 = 0.0;
    double e12 = 0.0;
    double e13 = 0.0;
    double e14 = 0.0;
    double ev1 = 0.0;
    double ev2 = 0.0;
    double ev3 = 0.0;
    double emin = 0.0;

    while (iLat < ag->nLats)
    {
        while (iLon < ag->nLons)
        {
            j = iLat;
            i = iLon;
            // Form the filename for the current chunk
            if (latChunk < 0)
                sprintf(latnum, "-%02d", -latChunk);
            else
                sprintf(latnum, "%02d", latChunk);
            if (lonChunk < 0)
                sprintf(lonnum, "-%03d", -lonChunk);
            else
                sprintf(lonnum, "%03d", lonChunk);
            sprintf(aster_filename, "%s/asterged.v003.%s.%s.0010.h5", 
                asterdir, latnum, lonnum);
            if (ag->verbose > 2)
                printf("    File: %s\n", aster_filename);

            // Determine if the file exists
            if (access(aster_filename, F_OK) == 0)
            {
                // File exists
                if (access(aster_filename, R_OK) != 0)
                {
                    // User does not have read access
                    fflush(stdout);
                    fprintf(stderr, "error -- No read access to file %s\n",
                        aster_filename);
                    return -1;
                }

                // Open the file
                hid_t agfid = open_hdf5(aster_filename, H5F_ACC_RDONLY);
                if (agfid == -1)
                {
                    // Error opening file
                    fflush(stdout);
                    fprintf(stderr, "error -- open failed for file %s\n",
                        aster_filename);
                    return -1;
                }
                
                // Read the data
                if (hdf5read_mat3d_int16(agfid, "/Emissivity/Mean", &echunk) < 0)
                {
                    close_hdf5(agfid);
                    return -1;
                }

                if (ag->verbose > 3)
                    printf("    Loaded /Emissivity/Mean %lu x %lu x %lu: "
                           "%d %d %d %d ...\n",
                            echunk.size1, echunk.size2, echunk.size3,
                            echunk.vals[0], echunk.vals[1],
                            echunk.vals[2], echunk.vals[3]);
                if (ag->read_ndvi)
                {
                    if (hdf5read_mat_int16(agfid, "/NDVI/Mean", &nchunk) < 0)
                    {
                        close_hdf5(agfid);
                        return -1;
                    }

                    if (ag->verbose > 3)
                        printf("    Loaded /NDVI/Mean %lu x %lu: "
                               "%d %d %d %d ...\n",
                                nchunk.size1, nchunk.size2,
                                nchunk.vals[0], nchunk.vals[1],
                                nchunk.vals[2], nchunk.vals[3]);
                }
                close_hdf5(agfid);
                ag->nFiles++;

                // Copy and convert chunk data into overall matrices.
                // Invert Latitude from high-low to low-high.
                for (y = echunk.size2 - 1; y >= 0; y--)
                {
                    for (x = 0; x < echunk.size3; x++)
                    {
                        if (ag->read_ndvi)
                        {
                            ndvival = mat_int16_get(&nchunk, y, x) 
                                      * ASTER_GED_NDVI_SCALE;
                            if (ndvival < ASTER_GED_NDVI_MIN)
                                ndvival = ASTER_NAN;
                            mat_set(&ag->ndvi, j, i, ndvival);
                            if (ag->verbose > 9)
                            {
                                printf("        set ndvi[%d,%d]=%.4f\n",
                                        j, i, ndvival);
                            }
                        }
                        if (ag->sensor == ECOSTRESS)
                        {
                            e12 =  mat3d_int16_get(&echunk, Band12, y, x)
                                      * ASTER_GED_EMIS_SCALE;
                        }
                        else
                        {
                            e11 = mat3d_int16_get(&echunk, Band11, y, x)
                                      * ASTER_GED_EMIS_SCALE;
                            if (e11 < ASTER_GED_EMIS_MIN)
                                e11 = ASTER_NAN;
                            e13 = mat3d_int16_get(&echunk, Band13, y, x)
                                      * ASTER_GED_EMIS_SCALE;
                            if (e13 < ASTER_GED_EMIS_MIN)
                                e13 = ASTER_NAN;
                            e14 = mat3d_int16_get(&echunk, Band14, y, x)
                                      * ASTER_GED_EMIS_SCALE;
                            if (e14 < ASTER_GED_EMIS_MIN)
                                e14 = ASTER_NAN;
                        }
                        if (e11 == ASTER_NAN || 
                            e13 == ASTER_NAN || 
                            e14 == ASTER_NAN)
                        {
                            ev1 = ASTER_NAN;
                            ev2 = ASTER_NAN;
                            ev3 = ASTER_NAN;
                        }
                        else if (ag->sensor == MODIS)
                        {
                            ev1 = 0.9878 * e11 + 0.0125;
                            ev2 = 0.4050 * e13 + 0.5945 * e14 + 0.0013;
                            ev3 = 0.1175 * e13 + 0.3572 * e14 + 0.5168;
                        }
                        else if (ag->sensor == VIIRS)
                        {
                            ev1 = 0.9902 * e11 + 0.0097;
                            ev2 = 0.8805 * e13 + 0.1375 * e14 - 0.0180;
                            ev3 = -0.3272 * e13 + 1.1799 * e14 + 0.1466;
                        }
                        else if (ag->sensor == ECOSTRESS)
                        {
                            ev3 = 0.9809 * e12 + 0.0176;
                        }
                        else
                        {
                            fprintf(stderr, "Invalid SENSOR type passed to ASTER_GED processing.\n");
                            return -1;
                        }
                        if (ag->sensor == ECOSTRESS)
                            emin = ev3;
                        else
                            emin = MIN(MIN(ev1, ev2), ev3);
                        mat_set(&ag->emis, j, i, emin);
                        if (ag->verbose > 9)
                        {
                            printf("        set emis[%d,%d]=%.4f\n",
                                    j, i, emin);
                        }
                        i++;
                    }
                    j++;
                    i = iLon;
                }

            }
            else
            {
                if (ag->verbose > 3)
                    printf("    Not found: used default values "
                           "emis=%.2f, NDVI=%.2f\n",
                           ASTER_GED_EMIS_DEFAULT,
                           ASTER_GED_NDVI_DEFAULT);
                ag->nDefault++;
                // The file does not exist -- apply default values for the chunk
                for (y = 0; y < 100; y++)
                {
                    for (x = 0; x < 100; x++)
                    {
                        mat_set(&ag->emis, j, i, ASTER_GED_EMIS_DEFAULT);
                        if (ag->read_ndvi)
                            mat_set(&ag->ndvi, j, i, ASTER_GED_NDVI_DEFAULT);
                        if (ag->verbose > 9)
                        {
                            printf("        set emis[%d,%d]=%.4f\n",
                                    j, i, ASTER_GED_EMIS_DEFAULT);
                            if (ag->read_ndvi)
                                printf("        set ndvi[%d,%d]=%.4f\n",
                                    j, i, ASTER_GED_NDVI_DEFAULT);
                        }
                        i++;
                    }
                    j++;
                    i = iLon;
                }
            }

            // Move to the next 100 x 100 chunk
            iLon += 100;
            lonChunk++;
        }
        // Move to the chunks for the next latitude  level
        iLon = 0;
        iLat += 100;
        lonChunk = (int)ag->startLon;
        latChunk++;
    }

    if (ag->verbose > 0)
    {
        printf("ag_load: Loaded %d files and used %d default chunks.\n",
                ag->nFiles, ag->nDefault);
    }
    if (ag->verbose > 4)
    {
        printf("         emis=%.4f %.4f %.4f %.4f ...\n",
                ag->emis.vals[0], ag->emis.vals[1],
                ag->emis.vals[2], ag->emis.vals[3]);
        if (ag->read_ndvi)
            printf("         ndvi=%.4f %.4f %.4f %.4f ...\n",
                ag->ndvi.vals[0], ag->ndvi.vals[1],
                ag->ndvi.vals[2], ag->ndvi.vals[3]);
    }
    mat3d_int16_clear(&echunk);
    mat_int16_clear(&nchunk);
    return 0; // Success
}

int read_aster_ged(const char *ASTER_dir, Matrix *gran_lat, 
    Matrix *gran_lon, Matrix *emis_aster, MatUint8 *oceanpix, 
    MatUint8 *snow_obs, Matrix *ndvi_obs, bool adjust_aster,
    Sensor sensor_type)
{
#ifdef ENABLE_TIMING
    time_t start_init = getusecs();
#endif
    // Main structure for Aster GED data input
    AsterGed ag;

    ag_init(&ag); // call constructor
    ag.sensor = sensor_type; 
    ag.read_ndvi = (ndvi_obs != NULL);

#ifdef DEBUG_ASTER
    ag.verbose = 9; // debug message level, default = 0
#endif

    // Check ASTER_dir. It must exist, have read and search privileges,
    // and also contain some files. Since ASTER is a fluid data set, 
    // not every lat/lon has a file. When a file is not found, the
    // ASTER load function assumes it is entirely ocean and supplies
    // default ocean values. If an empty directory was given, the loader
    // could not tell that it was wrong. It would generate ASTER values
    // using the default for Ocean for the entire granule. By checking for a 
    // few files, the risk of not detecting a bad directory name is reduced.
    if (access(ASTER_dir, F_OK) != 0)
    {
        fprintf(stderr, "ASTER_DIR does not exist -- %s\n", ASTER_dir);
        return -1;
    }

    if (access(ASTER_dir, R_OK | X_OK) != 0)
    {
        fprintf(stderr, "ASTER_DIR lacks read or search privileges -- %s\n", ASTER_dir);
        return -1;
    }

    // Don't read the entire directory because it would tke a long time, but
    // make sure there are at least 4 file sin the directory or it might not be 
    // ASTER data.
    DIR *astdir = opendir(ASTER_dir);
    struct dirent *pent;
    int nents = 0;
    const int AT_LEAST = 32;
    while (nents < AT_LEAST && (pent = readdir(astdir)) != NULL)
    {
        if (!strncmp(pent->d_name, "asterged", 8))
            nents++;
    }
    closedir(astdir);
    if (nents < AT_LEAST)
    {
        fprintf(stderr, "ASTER_DIR does not contain enough asterged files -- %s\n", ASTER_dir);
        return -1;
    }

    // Determine the ASTER limits to load.
    double minLat, maxLat;
    minmax2d(gran_lat->vals, gran_lat->size1, gran_lat->size2, -9999.0, 
        &minLat, &maxLat);
    double minLon, maxLon;
    minmax2d(gran_lon->vals, gran_lon->size1, gran_lon->size2, -9999.0, 
        &minLon, &maxLon);

    if (is_nan(minLat) || is_nan(minLon))
    {
        // This condition can happen is all the Latitude or Longitude data
        // is set to fill values.
        fprintf(stderr, "ASTER_GED failed because there are no valid Latitude/Longitude points in the granule.\n");
        return -1;
    }

#ifdef ENABLE_TIMING
    time_t elapsed_init = elapsed(start_init);
    printf("ASTERGED: init                 %9.3f secs\n",
        toseconds(elapsed_init));
#endif

#ifdef DEBUG_ASTER
    printf("\nSwath coordinates:\nLatitude:\n");
    int ie = gran_lat->size1 - 1;
    int je = gran_lat->size2 - 1;
    printf("%.2f %.2f %.2f\t...\t%.2f %.2f %.2f\n",
       mat_get(gran_lat, 0, 0), mat_get(gran_lat, 0, 1), mat_get(gran_lat, 0, 2),
       mat_get(gran_lat, 0, je-2), mat_get(gran_lat, 0, je-1), mat_get(gran_lat, 0, je) );
    printf("%.2f %.2f %.2f\t...\t%.2f %.2f %.2f\n",
       mat_get(gran_lat, 1, 0), mat_get(gran_lat, 1, 1), mat_get(gran_lat, 1, 2),
       mat_get(gran_lat, 1, je-2), mat_get(gran_lat, 1, je-1), mat_get(gran_lat, 1, je) );
    printf("%.2f %.2f %.2f\t...\t%.2f %.2f %.2f\n",
       mat_get(gran_lat, 2, 0), mat_get(gran_lat, 2, 1), mat_get(gran_lat, 2, 2),
       mat_get(gran_lat, 2, je-2), mat_get(gran_lat, 2, je-1), mat_get(gran_lat, 2, je) );
    printf("%.2f %.2f %.2f\t...\t%.2f %.2f %.2f\n",
       mat_get(gran_lat, ie-2, 0), mat_get(gran_lat, ie-2, 1), mat_get(gran_lat, ie-2, 2),
       mat_get(gran_lat, ie-2, je-2), mat_get(gran_lat, ie-2, je-1), mat_get(gran_lat, ie-2, je) );
    printf("%.2f %.2f %.2f\t...\t%.2f %.2f %.2f\n",
       mat_get(gran_lat, ie-1, 0), mat_get(gran_lat, ie-1, 1), mat_get(gran_lat, ie-1, 2),
       mat_get(gran_lat, ie-1, je-2), mat_get(gran_lat, ie-1, je-1), mat_get(gran_lat, ie-1, je) );
    printf("%.2f %.2f %.2f\t...\t%.2f %.2f %.2f\n",
       mat_get(gran_lat, ie, 0), mat_get(gran_lat, ie, 1), mat_get(gran_lat, ie, 2),
       mat_get(gran_lat, ie, je-2), mat_get(gran_lat, ie, je-1), mat_get(gran_lat, ie, je) );
    printf("\nLongitude:\n");
    printf("%.2f %.2f %.2f\t...\t%.2f %.2f %.2f\n",
       mat_get(gran_lon, 0, 0), mat_get(gran_lon, 0, 1), mat_get(gran_lon, 0, 2),
       mat_get(gran_lon, 0, je-2), mat_get(gran_lon, 0, je-1), mat_get(gran_lon, 0, je) );
    printf("%.2f %.2f %.2f\t...\t%.2f %.2f %.2f\n",
       mat_get(gran_lon, 1, 0), mat_get(gran_lon, 1, 1), mat_get(gran_lon, 1, 2),
       mat_get(gran_lon, 1, je-2), mat_get(gran_lon, 1, je-1), mat_get(gran_lon, 1, je) );
    printf("%.2f %.2f %.2f\t...\t%.2f %.2f %.2f\n",
       mat_get(gran_lon, 2, 0), mat_get(gran_lon, 2, 1), mat_get(gran_lon, 2, 2),
       mat_get(gran_lon, 2, je-2), mat_get(gran_lon, 2, je-1), mat_get(gran_lon, 2, je) );
    printf("%.2f %.2f %.2f\t...\t%.2f %.2f %.2f\n",
       mat_get(gran_lon, ie-2, 0), mat_get(gran_lon, ie-2, 1), mat_get(gran_lon, ie-2, 2),
       mat_get(gran_lon, ie-2, je-2), mat_get(gran_lon, ie-2, je-1), mat_get(gran_lon, ie-2, je) );
    printf("%.2f %.2f %.2f\t...\t%.2f %.2f %.2f\n",
       mat_get(gran_lon, ie-1, 0), mat_get(gran_lon, ie-1, 1), mat_get(gran_lon, ie-1, 2),
       mat_get(gran_lon, ie-1, je-2), mat_get(gran_lon, ie-1, je-1), mat_get(gran_lon, ie-1, je) );
    printf("%.2f %.2f %.2f\t...\t%.2f %.2f %.2f\n",
       mat_get(gran_lon, ie, 0), mat_get(gran_lon, ie, 1), mat_get(gran_lon, ie, 2),
       mat_get(gran_lon, ie, je-2), mat_get(gran_lon, ie, je-1), mat_get(gran_lon, ie, je) );
    printf("\n");
#endif

    // Load Aster GED data
#ifdef ENABLE_TIMING
    time_t start_agload = getusecs();
#endif

    if (ag_load(&ag, ASTER_dir, minLat, minLon, maxLat, maxLon) < 0)
    {
        return -1;
    }

#ifdef ENABLE_TIMING
    time_t elapsed_agload = elapsed(start_agload);
    printf("ASTERGED: ag_load              %9.3f secs\n",
        toseconds(elapsed_agload));

    time_t start_swmap = getusecs();
#endif

    // Create Lat/Lon map for ASTER data
    Matrix asterLat, asterLon;
    mat_init(&asterLat);
    mat_init(&asterLon);
    meshgrat(&asterLat, &asterLon, ag.startLat, ag.startLon, 
             ag.endLat, ag.endLon, 0.01, 0.01);

    // Map the ASTER data onto the swath
    int nsets = 1;
    Matrix griddedEmis;
    mat_init(&griddedEmis);
    mat_alloc(&griddedEmis, gran_lat->size1, gran_lat->size2);
    Matrix griddedNdvi;
    mat_init(&griddedNdvi);
    if (ag.read_ndvi)
    {
        mat_alloc(&griddedNdvi, gran_lat->size1, gran_lat->size2);
        nsets = 2;
    }
    double *insets[] = {ag.emis.vals, ag.ndvi.vals};
    double *outsets[] = {griddedEmis.vals, griddedNdvi.vals};
    multi_interp2(asterLat.vals, asterLon.vals, asterLat.size1, asterLat.size2,
        gran_lat->vals, gran_lon->vals, gran_lat->size1*gran_lat->size2,
        insets, outsets, nsets);

#ifdef ENABLE_TIMING
    time_t elapsed_swmap = elapsed(start_swmap);
    printf("ASTERGED: map to swath         %9.3f secs\n",
        toseconds(elapsed_swmap));
#endif

#ifdef ENABLE_TIMING
    time_t start_fill = getusecs();
#endif

    Matrix emodfill;
    mat_init(&emodfill);
    mat_copy(&emodfill, &griddedEmis);
    
    Matrix Nfill;
    mat_init(&Nfill);
    if (ag.read_ndvi)
        mat_copy(&Nfill, &griddedNdvi);

#ifdef USE_INTERP_ID
    // Interpolate to remove NaNs.
    interp_id_columns(emodfill.vals, emodfill.size1, emodfill.size2);
    interp_id_columns(Nfill.vals, Nfill.size1, Nfill.size2);
#else
    // Alternative to interpolate: use default values
    int ifill;
    for (ifill = 0; ifill < emodfill.size1*emodfill.size2; ifill++)
    {
        if (is_nan(emodfill.vals[ifill]))
            emodfill.vals[ifill] = ASTER_GED_EMIS_NOMINAL;
        if (ag.read_ndvi)
        {
            if (is_nan(Nfill.vals[ifill]))
                Nfill.vals[ifill] = ASTER_GED_NDVI_NOMINAL;
        }
    }
#endif

#ifdef ENABLE_TIMING
    time_t elapsed_fill = elapsed(start_fill);
    printf("ASTERGED: fill NaNs            %9.3f secs\n",
        toseconds(elapsed_fill));

    time_t start_obs = getusecs();
#endif

    // ASTER adjustment for observation phenology and snow
    mat_alloc(emis_aster, emodfill.size1, emodfill.size2);
    
    int i, j;
    int nocean = 0;
    int nsnow = 0;
    double ndvi_val;
    const double mintest = 0.15;
    for (i = 0; i < emodfill.size1; i++)
    {
        for (j = 0; j < emodfill.size2; j++)
        {
            if (ag.read_ndvi)
            {
                // Limit NDVI to minimum of zero.
                if (mat_get(&Nfill, i, j) < 0.0)
                    mat_set(&Nfill, i, j, 0.0);
            }

            // Set emis to 0.99 over the ocean
            if (mat_uint8_get(oceanpix, i, j) != 0)
            {
                mat_set(&emodfill, i, j, ASTER_GED_EMIS_OCEAN);
                mat_set(emis_aster, i, j, ASTER_GED_EMIS_OCEAN);
                if (ag.read_ndvi)
                    mat_set(&Nfill, i, j, ASTER_GED_NDVI_OCEAN);
                nocean++;
            }
            else if (adjust_aster && ag.read_ndvi)
            {
                // Adjust EMIS
                double fv_swath  = 
                    1.0 - ((1.0 - mat_get(&Nfill, i, j))/(1.0 - mintest));
                if (fv_swath < 0.0) fv_swath = 0.0;

                ndvi_val = mat_get(ndvi_obs, i, j);
                if (is_nan(ndvi_val))
                    ndvi_val = mat_get(&Nfill, i, j);
                double fv_obs = 
                    1.0 - ((1.0 - ndvi_val)/(1.0 - mintest));
                if (fv_obs < 0.0) fv_obs = 0.0;
#ifdef DEBUG_EMIS_NANS
                if (is_nan(fv_obs))
                {
                    printf("fv_obs = 1 - ((1 - ndvi_val(%d,%d)(=%.4f))/(0.85)\n",
                            i, j, ndvi_val);
                }
#endif
                double emis_bare = 
                    (mat_get(&emodfill, i, j) - 0.985 * fv_swath)/
                    (1.0 - fv_swath);
                double emis_aster_val =  0.985 * fv_obs + emis_bare * (1.0 - fv_obs);
                if (mat_uint8_get(snow_obs, i, j) == 1)
                {
                    emis_aster_val = 0.985;  // Medium snow
                    nsnow++;
                }
                else if (emis_aster_val < ASTER_GED_EMIS_MIN)
                {
                    //emis_aster_val = _NaN_; // Don't use this pixel
                    // [BF] Instead, use nominal value
                    emis_aster_val = ASTER_GED_EMIS_NOMINAL;
                }
                else if (emis_aster_val > ASTER_GED_EMIS_MAX)
                {
                    emis_aster_val = ASTER_GED_EMIS_MAX; // Clamp max val
                }
                mat_set(emis_aster, i, j, emis_aster_val);
#ifdef DEBUG_EMIS_NANS
                if (is_nan(emis_aster_val))
                {
                    printf("emis_aster(%d, %d) set to NaN because of calulation:"
                           "fv_swath=%.4f fv_obs=%.4f emis_bare=%.4f\n",
                            i, j, fv_swath, fv_obs, emis_bare);
                }
#endif
            }
            else // no adjustment
            {
                 mat_set(emis_aster, i, j, 
                     mat_get(&emodfill, i, j));
            }
        }
    }

#ifdef DEBUG_ASTER
    printf("ASTER GED: %d ocean pixels, %d snow pixels\n", nocean, nsnow);
#endif

#ifdef ENABLE_TIMING
    time_t elapsed_obs = elapsed(start_obs);
    printf("ASTERGED: adjustments          %9.3f secs\n",
        toseconds(elapsed_obs));
#endif

#ifdef DEBUG_OUTPUTS

#ifdef ENABLE_TIMING
    time_t start_outs = getusecs();
#endif

    // Output the datasets to be available in AsterGed.hdf
    //     MATLAB         aster_ged_4modis.c
    //     =============  ===================
    //     emis_tot       ag.emis              Loaded data from Aster collection
    //     NDVI_tot       ag.ndvi              Loaded data from Aster collection
    //     Late_sub       asterLat             meshgrat Lat xcoiodinates
    //     Lone_sub       asterLon             meshgrat Lat xcoiodinates
    //     emodf          griddedEmis          Aster gridded EMIS
    //     Nf             griddedNdvi          Aster gridded NDVI
    //     emodfill       emodfill             gridded and NaNs interpolated out
    //     Nfill          Nfill                gridded and NaNs interpolated out
    //     emis_sub       n/a                  same as griddedEmis
    //     NDVI_sub       n/a                  same as griddedNdvi
    //     emis_aster     emis_aster           Final adjusted swath emis data

    int stattot = 0;
    int32 out_fid = SDstart("AsterGed.hdf", DFACC_CREATE);
    if (out_fid == FAIL)
    {
        fprintf(stderr, "Failed to create AsterGed.hdf\n");
        return 0;
    }

    //stattot += hdfwrite_mat(out_fid, "emis_tot", &ag.emis);
    //if (ag.ndvi != NULL)
    //    stattot += hdfwrite_mat(out_fid, "NDVI_tot", &ag.ndvi);
    //stattot += hdfwrite_mat(out_fid, "Late_sub", &asterLat);
    //stattot += hdfwrite_mat(out_fid, "Lone_sub", &asterLon);
    //stattot += hdfwrite_mat(out_fid, "emodf", &griddedEmis);
    //stattot += hdfwrite_mat(out_fid, "Nf", &griddedNdvi);
    stattot += hdfwrite_mat(out_fid, "emodfill", &emodfill);
    if (ag.read_ndvi)
        stattot += hdfwrite_mat(out_fid, "Nfill", &Nfill);
    stattot += hdfwrite_mat(out_fid, "emis_aster", emis_aster);
#ifdef DEBUG_ASTER
    stattot += hdfwrite_mat(out_fid, "SwathLatitude", gran_lat);
    stattot += hdfwrite_mat(out_fid, "SwathLongitude", gran_lon);
#endif
    SDend(out_fid);
    if (stattot < 0)
    {
        fprintf(stderr, "Error writing %d datasets to AsterGed.hdf\n",
                -stattot);
    }

#ifdef ENABLE_TIMING
    time_t elapsed_outs = elapsed(start_outs);
    printf("ASTERGED: AsterGed.hdf output  %9.3f secs\n",
        toseconds(elapsed_outs));
#endif
#endif // DEBUG_OUTPUTS

    mat_clear(&emodfill);
    mat_clear(&griddedEmis);
    mat_clear(&griddedNdvi);
    mat_clear(&asterLat);
    mat_clear(&asterLon);

    return 0;
}
