// See details in cloud.h
//
#include "cloud.h"
#include "error.h"
#include "interps.h"
#include "matrix.h"
#include "metadata.h"
#include "tes_util.h"
#include "lste_lib.h"
#include "fileio.h"

#include <string.h>
#include <stdio.h>

void process_cloud(RAD *rad,
                   const char* bt11_lut_file,
                   const char* btdiff_file,
                   double* rad_lut_data[],
                   int rad_lut_nlines,
                   const char* cloud_filename)
{
    const unsigned int BAND_4 = rad->nchannels - 2;
    const unsigned int BAND_5 = BAND_4 + 1;
    const unsigned int COL_1 = 0;
    const unsigned int COL_5 = 4;
    const unsigned int COL_6 = 5;

    // Convert radiances to brightness temperature
    unsigned int nrows = rad->Rad[0].size1;
    unsigned int ncols = rad->Rad[0].size2;
    Matrix TB4; mat_init(&TB4);
    mat_alloc(&TB4, nrows, ncols);
    Matrix TB5; mat_init(&TB5);
    mat_alloc(&TB5, nrows, ncols);
    int i;
    double v;
    for (i = 0; i < nrows * ncols; i++)
    {
        v = interp1d_npts(rad_lut_data[COL_5], rad_lut_data[COL_1],
            rad->Rad[BAND_4].vals[i], rad_lut_nlines);
        TB4.vals[i] = v;
        v = interp1d_npts(rad_lut_data[COL_6], rad_lut_data[COL_1],
            rad->Rad[BAND_5].vals[i], rad_lut_nlines);
        TB5.vals[i] = v;
    }

    // Cloud test 1: Brightness temperature look-up-table (LUT) approach
    hid_t h5fid = open_hdf5(bt11_lut_file, H5F_ACC_RDONLY);
    if (h5fid < 0)
        ERREXIT(14, "Unable to open lookup table from %s", bt11_lut_file);
    int stat = -1;

    Matrix t_lut_lat; mat_init(&t_lut_lat);
    stat = hdf5read_mat(h5fid, "/Geolocation/Latitude", &t_lut_lat);
    if (stat < 0)
        ERREXIT(14, "Error reading /Geolocation/Latitude from %s", bt11_lut_file);
    Matrix lut_lat; mat_init(&lut_lat);
    mat_transpose(&lut_lat, &t_lut_lat);
    mat_clear(&t_lut_lat);

    Matrix t_lut_lon; mat_init(&t_lut_lon);
    stat = hdf5read_mat(h5fid, "/Geolocation/Longitude", &t_lut_lon);
    if (stat < 0)
        ERREXIT(14, "Error reading /Geolocation/Longitude from %s", bt11_lut_file);
    Matrix lut_lon; mat_init(&lut_lon);
    mat_transpose(&lut_lon, &t_lut_lon);
    mat_clear(&t_lut_lon);

    // Get ECOSTRESS observation time info
    // sample filename: ECOSTRESS_L2_CLOUD_10390_027_20200506T002649_0601_01.h5
    //                           0^                    23^  28^
    const char* start_pos = rindex(cloud_filename, '/');
    if (start_pos == NULL)
        start_pos = cloud_filename;
    const char* L2_CLOUD = strstr(start_pos, "L2_CLOUD");
    if (L2_CLOUD == NULL)
        ERREXIT(14, "L2_CLOUD not in filename %s", cloud_filename);
    int hh, mm, mth;
    sscanf(L2_CLOUD+23, "%2d", &mth);
    sscanf(L2_CLOUD+28, "%2d%2d", &hh, &mm);

    double hrfrac = hh + (double)mm/60.0;

    // Find nearest Cloud LUT UTC time
    int Ztime0 = 6 * (hh / 6);
    int Ztime1 = (Ztime0 + 6) % 24;
    char cloudvar1[32];
    char cloudvar2[32];
    sprintf(cloudvar1, "/Data/LUT_cloudBT11_%02d_%02d", Ztime0, mth);
    sprintf(cloudvar2, "/Data/LUT_cloudBT11_%02d_%02d", Ztime1, mth);

    Matrix t_BT1; mat_init(&t_BT1);
    stat = hdf5read_mat(h5fid, cloudvar1, &t_BT1);
    if (stat < 0)
        ERREXIT(14, "Error reading %s from %s", cloudvar1, bt11_lut_file);
    Matrix BT1; mat_init(&BT1);
    mat_transpose(&BT1, &t_BT1);
    mat_clear(&t_BT1);

    Matrix t_BT2; mat_init(&t_BT2);
    stat = hdf5read_mat(h5fid, cloudvar2, &t_BT2);
    if (stat < 0)
        ERREXIT(14, "Error reading %s from %s", cloudvar2, bt11_lut_file);
    Matrix BT2; mat_init(&BT2);
    mat_transpose(&BT2, &t_BT2);
    mat_clear(&t_BT2);
    close_hdf5(h5fid);

    // Grid thresholds onto ECOSTRESS scene
    Matrix BTgrid1; mat_init(&BTgrid1);
    mat_alloc(&BTgrid1, nrows,  ncols);
    Matrix BTgrid2; mat_init(&BTgrid2);
    mat_alloc(&BTgrid2, nrows,  ncols);
    double *inds[2] = {BT1.vals, BT2.vals};
    double *outds[2] = {BTgrid1.vals, BTgrid2.vals};
    multi_interp2(lut_lat.vals, lut_lon.vals, lut_lat.size1, lut_lat.size2,
        rad->Lat.vals, rad->Lon.vals, rad->Lat.size1 * rad->Lat.size2,
        inds, outds, 2);
    double m = (hrfrac - Ztime0) / 6.0;
    Matrix BT; mat_init(&BT);
    mat_alloc(&BT, nrows, ncols);
    for (i = 0; i < nrows * ncols; i++)
    {
        BT.vals[i] = BTgrid1.vals[i] + m * (BTgrid2.vals[i] -  BTgrid1.vals[i]);
    }

    // Detected clouds have lower brightness temperatures than thresholds
    MatUint8 cloud1;
    mat_uint8_init(&cloud1);
    mat_uint8_alloc(&cloud1, nrows, ncols);
    for (i = 0; i < nrows * ncols; i++)
    {
        if (isnan(TB4.vals[i]))
            cloud1.vals[i] = 255;
        else
            cloud1.vals[i] = TB4.vals[i] < BT.vals[i];
    }

    // Cleanup cloud test 1
    mat_clear(&lut_lat);
    mat_clear(&lut_lon);
    mat_clear(&BT1);
    mat_clear(&BT2);
    mat_clear(&BTgrid1);
    mat_clear(&BTgrid2);
    mat_clear(&BT);

    // Cloud test 2

    // Read Cloud BTdiff data
    Vector T11; vec_init(&T11);
    Matrix Tmax; mat_init(&Tmax);
    Vector VA;  vec_init(&VA);
    h5fid = open_hdf5(btdiff_file, H5F_ACC_RDONLY);
    if (h5fid == -1)
    {
        ERREXIT(14, "Unable to open cloud BTdiff file %s", btdiff_file);
    }
    stat = hdf5read_vec(h5fid, "/SDS/T11", &T11);
    if (stat == -1)
    {
        ERREXIT(14, "Unable to read /SDS/T11 from %d", btdiff_file);
    }
    stat = hdf5read_mat(h5fid, "/SDS/Tmax", &Tmax);
    if (stat == -1)
    {
        ERREXIT(14, "Unable to read /SDS/Tmax from %d", btdiff_file);
    }
    stat = hdf5read_vec(h5fid, "/SDS/VA", &VA);
    if (stat == -1)
    {
        ERREXIT(14, "Unable to read /SDS/VA from %d", btdiff_file);
    }
    close_hdf5(h5fid);

    // [1.1.2] Reduce Tmax by 10%
    int iTmax;
    for (iTmax = 0; iTmax < Tmax.size1 * Tmax.size2; iTmax++)
    {
        Tmax.vals[iTmax] *= 0.9;
    }

    // Set up BT diff meshgrid for interpolants
    int btdiff_nrows = T11.size - 1;
    int btdiff_ncols = VA.size;
    Matrix btdiff_X; mat_init(&btdiff_X);
    mat_alloc(&btdiff_X, btdiff_nrows, btdiff_ncols);
    Matrix btdiff_Y; mat_init(&btdiff_Y);
    mat_alloc(&btdiff_Y, btdiff_nrows, btdiff_ncols);
    int btrow, btcol;
    for (btrow = 0; btrow < btdiff_nrows; btrow++)
    {
        for (btcol = 0; btcol < btdiff_ncols; btcol++)
        {
            mat_set(&btdiff_X, btrow, btcol, T11.vals[btrow+1]);
            mat_set(&btdiff_Y, btrow, btcol, VA.vals[btcol]);
        }
    }

    Matrix TB4temp; mat_init(&TB4temp);
    mat_alloc(&TB4temp, nrows, ncols);
    Matrix SZtemp; mat_init(&SZtemp);
    mat_alloc(&SZtemp, nrows, ncols);
    int nquery_points = nrows * ncols;
    for (i = 0; i < nquery_points; i++)
    {
        v = TB4.vals[i];
        if (v < 260.0) v = 260.0;
        TB4temp.vals[i] = v;
        v = rad->Satzen.vals[i];
        if (v < 3.0)
            v = 3.0;
        SZtemp.vals[i] = v;
    }

    // Interpolate Tmax onto granule based on TB4 and SZ
    Matrix btdifft_mat; mat_init(&btdifft_mat);
    mat_alloc(&btdifft_mat, nrows, ncols);
    interp2(btdiff_X.vals, btdiff_Y.vals, Tmax.vals, btdiff_nrows, btdiff_ncols,
            TB4temp.vals, SZtemp.vals, btdifft_mat.vals, nquery_points);

    double Tdiff1, btdifft1;
    MatUint8 cloud2;
    mat_uint8_init(&cloud2);
    mat_uint8_alloc(&cloud2, nrows, ncols);
    for (i = 0; i < nquery_points; i++)
    {
        Tdiff1 = TB4.vals[i] - TB5.vals[i];
        // Increase threshold above 1km elevation to account for
        // elevation dependent anomalies.
        btdifft1 = btdifft_mat.vals[i] + 1.0;
        if (isnan(Tdiff1) || isnan(btdifft1))
            cloud2.vals[i] = 255;
        else
            cloud2.vals[i] = Tdiff1 > btdifft1;
    }

    // Contiguity  - set cloudy pixels with 3 or less cloudy neighbors to clear
    const int NumNeighbors = 3;
    constrict_cloud(&cloud2, NumNeighbors);
    // Repeat a second time.
    constrict_cloud(&cloud2, NumNeighbors);

    // Final cloud
    MatUint8 cloud_final;
    mat_uint8_init(&cloud_final);
    mat_uint8_alloc(&cloud_final, nrows, ncols);
    for (i = 0; i < nquery_points; i++)
    {
        cloud_final.vals[i] = cloud1.vals[i] | cloud2.vals[i];
    }

    // Cleanup intermediate matrices.
    mat_clear(&TB4);
    mat_clear(&TB5);
    mat_clear(&SZtemp);
    mat_clear(&TB4temp);
    mat_clear(&btdiff_X);
    mat_clear(&btdiff_Y);
    mat_clear(&btdifft_mat);

    // Output the Cloud_test_1 dataset.
    hid_t cloudout = create_hdf5(cloud_filename);
    if (cloudout < 0)
        ERREXIT(16, "Could not create output file %s", cloud_filename);

    hid_t sds_group = H5Gcreate2(cloudout, "/SDS",
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (sds_group == -1)
        ERREXIT(17, "Could not create group /SDS in %s", cloud_filename);
    H5Gclose(sds_group);

    stat = hdf5write_mat_uint8(cloudout, "/SDS/Cloud_test_1", &cloud1);
    if (stat == -1)
        ERREXIT(17, "Error writing /SDS/Cloud_test_1 to %s", cloud_filename);

    // Add attributes to the output dataset.
    uint8 minCM = 0;
    uint8 maxCM = 1;
    uint8 fillCM = 255;
    stat = writeDatasetMetadataHdf5(cloudout, "/SDS/Cloud_test_1",
            "Brightness temperature LUT test",
            "1 = cloudy",
            DFNT_UINT8, &minCM, &maxCM,
            DFNT_UINT8, &fillCM, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(17, "Failed to write attributes to dataset "
            "/SDS/Cloud_test_1 in %s",
            cloud_filename);

    // Output Cloud_test_2.
    stat = hdf5write_mat_uint8(cloudout, "/SDS/Cloud_test_2", &cloud2);
    if (stat == -1)
        ERREXIT(17, "Error writing /SDS/Cloud_test_2 to %s", cloud_filename);

    // Add attributes to the output dataset.
    stat = writeDatasetMetadataHdf5(cloudout, "/SDS/Cloud_test_2",
            "Brightness temperature difference test",
            "1 = cloudy",
            DFNT_UINT8, &minCM, &maxCM,
            DFNT_UINT8, &fillCM, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(17, "Failed to write attributes to dataset "
            "/SDS/Cloud_test_2 in %s",
            cloud_filename);

    // Output Cloud_test_final.
    stat = hdf5write_mat_uint8(cloudout, "/SDS/Cloud_test_final", &cloud_final);
    if (stat == -1)
        ERREXIT(17, "Error writing /SDS/Cloud_test_final to %s", cloud_filename);

    // Add attributes to the output dataset.
    stat = writeDatasetMetadataHdf5(cloudout, "/SDS/Cloud_test_final",
            "Final cloud mask",
            "1 = cloudy",
            DFNT_UINT8, &minCM, &maxCM,
            DFNT_UINT8, &fillCM, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(17, "Failed to write attributes to dataset "
            "/SDS/Cloud_test_final in %s",
            cloud_filename);

    close_hdf5(cloudout);

    // Cleanup data matrices.
    mat_uint8_clear(&cloud1);
    mat_uint8_clear(&cloud2);
    mat_uint8_clear(&cloud_final);
}
