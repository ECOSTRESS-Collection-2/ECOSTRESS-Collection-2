// tg_wvs.c
// See tg_wvs.h for credits and documentation

#include "tg_wvs.h"
#include "interps.h"
#include "lste_lib.h"
#include "fileio.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>

#define DEBUG
#ifdef DEBUG
#include "describe.h"
#endif

//#define DEBUG_HDF
#ifdef DEBUG_HDF
#include "error.h"
#endif

extern int debug_lste_lib;

//#define DEBUG_LINE 1200
//#define DEBUG_PIXEL 2740

#define MAX_CHANNELS 5
extern int n_channels; // Get n channels from tes_main.
extern int* band; // Get list of band numbers from tes_main.

// Compare function to return values for descending order.
int descending(const void *v1, const void *v2)
{
    double d1 = *(const double *)v1;
    double d2 = *(const double *)v2;
    if (d1 < d2) return 1;  // return "greater" when less
    if (d1 > d2) return -1; // return "less" when greater
    return 0;
}

int tg_wvs(Mat3d *tg, Matrix *pwv, Matrix *senszen,
           Matrix *emis_aster, Mat3d *lrad, Matrix *solzen, 
           const char *coef_file, double *lut[], int nlut_lines)
{
    int i, j, k;
    int stat = 0;

    // This function might need some extra stack space.
    struct rlimit stack_limits;
    getrlimit(RLIMIT_STACK, &stack_limits);
    stack_limits.rlim_cur = stack_limits.rlim_max;
    setrlimit(RLIMIT_STACK, &stack_limits);

#ifdef DEBUG_PIXEL
    printf("tg_wvs: (%d,%d) pwv=%g senszen=%g solzen=%g emis_aster=%g\n        lrad=[",
            DEBUG_LINE, DEBUG_PIXEL, 
            mat_get(pwv, DEBUG_LINE, DEBUG_PIXEL),
            mat_get(senszen, DEBUG_LINE, DEBUG_PIXEL),
            mat_get(solzen, DEBUG_LINE, DEBUG_PIXEL),
            mat_get(emis_aster, DEBUG_LINE, DEBUG_PIXEL));
    for (i = 0; i < n_channels; i++)
    {
        printf("%g", mat3d_get(lrad, i, DEBUG_LINE, DEBUG_PIXEL));
        if (i < n_channels-1) printf(",");
    }
    printf("]\n        coef_file=%s\n", coef_file);
#endif

    // Each band has three primary coefficients and a set of three 
    // secondary coefficients for each band.
    const int N_COEFFS_PER_BAND = 3 + (MAX_CHANNELS * 3);

    // Make a copy of emis_aster
    Matrix emis_aster_loc;
    mat_init(&emis_aster_loc);
    mat_copy(&emis_aster_loc, emis_aster);
    
    // Prepare output matrix
    mat3d_alloc(tg, n_channels, pwv->size1, pwv->size2);

    // Input the coefficient file.
    hid_t coef_fid = open_hdf5(coef_file, H5F_ACC_RDONLY);
    if (coef_fid == -1)
    {
        perror(coef_file);
        return -1;
    }

    // The file has MAX_CHANNELS coefficient datasets.
    Matrix coeff[MAX_CHANNELS];
    char dsname[32];
    stat = 0;
    for (i = 0; i < MAX_CHANNELS; i++)
    {
        mat_init(&coeff[i]);
        // Make the dataset name
        sprintf(dsname, "SDS/COEFF%d", i+1);
        stat += hdf5read_mat(coef_fid, dsname, &coeff[i]);
    }
    close_hdf5(coef_fid);
    if (stat < 0)
    {
        fprintf(stderr, "%d error(s) reading coefficient datasets from %s\n",
                -stat, coef_file);
        return -1;
    }

#ifdef DEBUG
    printf("tg_wvs: input coefficients %lu x %lu from file %s\n", 
            coeff[0].size1, coeff[0].size2, coef_file);
#endif

    // Determine day or night
    int totday = 0;
    int totpix = solzen->size1 * solzen->size2;
    for (i = 0; i < totpix; i++)
    {
        if (solzen->vals[i] > 85.0)
            totday++;
    }
    
    // The coefficient files have 17 rows. The first row of the data is the
    // day/night index. The first half of the coefficient columns are day(0)
    // and the second half are night(1). Easiest way in C is to add an offset
    // to the index.
    int ncoeffs = coeff[0].size2 / 2;
    int dnoff = 0;
    if ((double)totday / (double)totpix > 0.5)
    {
        // night -- index to second half of coefficients
        dnoff += ncoeffs;
    }

    // Extract the vectors for emisi, viewsi, and cwi
    Matrix coeff_in[MAX_CHANNELS];
    for (i = 0; i < MAX_CHANNELS; i++)
    {
        mat_init(&coeff_in[i]);
        mat_alloc(&coeff_in[i], N_COEFFS_PER_BAND, ncoeffs);
    }

    double emis_val;
    Vector emisi, viewsi, cwi;
    vec_init(&emisi);
    vec_init(&viewsi);
    vec_init(&cwi);
    vec_alloc(&emisi, ncoeffs);
    vec_alloc(&viewsi, ncoeffs);
    vec_alloc(&cwi, ncoeffs);
    int nval_coeffs = 0;
    // Loop over each coefficient
    for (i = 0; i < ncoeffs; i++)
    {
        emis_val = mat_get(&coeff[0], 2, i+dnoff);
        emis_val = round2(emis_val, 0.01);
        if (emis_val < 0.7)
            continue; // skip this column of coefficients
        // Accept column of coefficients
        vec_set(&emisi, nval_coeffs, emis_val);
        vec_set(&viewsi, nval_coeffs, mat_get(&coeff[0], 1, i+dnoff));
        vec_set(&cwi, nval_coeffs, mat_get(&coeff[0], 3, i+dnoff));
        // Loop over the values for each coefficient set. 
        for (j = 0; j < N_COEFFS_PER_BAND; j++)
        {
            // Loop over the MAX_CHANNELS
            for (k = 0; k < MAX_CHANNELS; k++)
            {
                mat_set(&coeff_in[k], j, nval_coeffs, 
                        mat_get(&coeff[k], j+4, i+dnoff));
            }
        }
        nval_coeffs++;
    }

    // Sort the datasets.
    Vector esort, vsort, csort;
    vec_init(&esort);
    vec_init(&vsort);
    vec_init(&csort);
    vec_copy(&esort, &emisi);
    vec_copy(&vsort, &viewsi);
    vec_copy(&csort, &cwi);
    qsort(esort.vals, nval_coeffs, sizeof(double), descending);
    qsort(vsort.vals, nval_coeffs, sizeof(double), descending);
    qsort(csort.vals, nval_coeffs, sizeof(double), descending);

    // Get the unique values
    int ne_uniq = 0;
    int nv_uniq = 0;
    int nc_uniq = 0;
    double laste = MATRIX_NAN;
    double lastv = MATRIX_NAN;
    double lastc = MATRIX_NAN;
    Vector emis, views, cw;
    vec_init(&emis);
    vec_init(&views);
    vec_init(&cw);
    vec_alloc(&emis, nval_coeffs);
    vec_alloc(&views, nval_coeffs);
    vec_alloc(&cw, nval_coeffs);
    double v;
    for (i = 0; i < nval_coeffs; i++)
    {
        v = vec_get(&esort, i);
        if (v != laste)
        {
            vec_set(&emis, ne_uniq++, v);
            laste = v;
        }
        v = vec_get(&vsort, i);
        if (v != lastv)
        {
            vec_set(&views, nv_uniq++, v);
            lastv = v;
        }
        v = vec_get(&csort, i);
        if (v != lastc)
        {
            vec_set(&cw, nc_uniq++, v);
            lastc = v;
        }
    }

#ifdef DEBUG
    printf("tg_wvs: emis[%d]  =", ne_uniq);
    for (k = 0; k < ne_uniq; k++)
        printf(" %.4f", vec_get(&emis, k));
    printf("\n");
    printf("tg_wvs: views[%d] =", nv_uniq);
    for (k = 0; k < nv_uniq; k++)
        printf(" %.4f", vec_get(&views, k));
    printf("\n");
    printf("tg_wvs: cw[%d]    =", nc_uniq);
    for (k = 0; k < nc_uniq; k++)
        printf(" %.4f", vec_get(&cw, k));
    printf("\n");
#endif

    // Setup the meshgrid for interpolants (sample indexes)
    Mat3d x, y, z;
    mat3d_init(&x);
    mat3d_init(&y);
    mat3d_init(&z);

    mat3d_alloc(&x, nc_uniq, nv_uniq, ne_uniq);
    mat3d_alloc(&y, nc_uniq, nv_uniq, ne_uniq);
    mat3d_alloc(&z, nc_uniq, nv_uniq, ne_uniq);
    for (i = 0; i < nc_uniq; i++)
    {
        for (j = 0; j < nv_uniq; j++)
        {
            for (k = 0; k < ne_uniq; k++)
            {
                mat3d_set(&x, i, j, k, vec_get(&views, j));
                mat3d_set(&y, i, j, k, vec_get(&emis, k));
                mat3d_set(&z, i, j, k, vec_get(&cw, i));
            }
        }
    }

    // Create N_COEFFS_PER_BAND coefficient sample bands for each channel
    Mat3d cband[MAX_CHANNELS][N_COEFFS_PER_BAND];
    for (i = 0; i < MAX_CHANNELS; i++)
    {
        for (k = 0; k < N_COEFFS_PER_BAND; k++)
        {
            mat3d_init(&cband[i][k]);
            mat3d_alloc(&cband[i][k], nc_uniq, nv_uniq, ne_uniq);
        }
    }

    int want;
    int found = 0;
    double findx, findy, findz;
    int ilev, irow, icol;
    for (ilev = 0; ilev < nc_uniq; ilev++)
    {
        for (irow = 0; irow < nv_uniq; irow++)
        {
            for (icol = 0; icol < ne_uniq; icol++)
            {
                found = 0;
                findx = mat3d_get(&x, ilev, irow, icol);
                findy = mat3d_get(&y, ilev, irow, icol);
                findz = mat3d_get(&z, ilev, irow, icol);
                for (want = 0; want < nval_coeffs; want++)
                {
                    if (viewsi.vals[want] == findx &&
                        emisi.vals[want] == findy &&
                        cwi.vals[want] == findz)
                    {
                        found = 1;
                        for (i = 0; i < MAX_CHANNELS; i++)
                        {
                            for (k = 0; k < N_COEFFS_PER_BAND; k++)
                            {
                                mat3d_set(&cband[i][k], ilev, irow, icol, 
                                    mat_get(&coeff_in[i], k, want));
                            }
                        }
                        break;
                    }
                }
                if (found < 1)
                {
                    printf("Failed to find entry for:\n");
                    printf("    X=%.4f  Y=%.4f  Z=%.4f\n",
                            findx, findy, findz);
                }
                assert(found > 0); // coordinates should find a match
            }
        }
    }

    // Clamp values in emis_aster to 0.7 to max emis value
    clamp2d(&emis_aster_loc, 0.7, vec_get(&emis, 0));

    // Scale water vapor coefficients by 1/10
    for (i = 0; i < z.size1 * z.size2 * z.size3; i++)
    {
        z.vals[i] /= 10.0;
    }

    // Convert NaNs to zeroes in PWV
    for (i = 0; i < pwv->size1 * pwv->size2; i++)
    {
        if (is_nan(pwv->vals[i]))
            pwv->vals[i] = 0.0;
    }

    // Interpolation of c_wvs for N_COEFFS_PER_BAND coefficients in MAX_CHANNELS bands
    Matrix c_wvs[MAX_CHANNELS][N_COEFFS_PER_BAND];
    for (i = 0; i < MAX_CHANNELS; i++)
    {
        for (j = 0; j < N_COEFFS_PER_BAND; j++)
        {
            mat_init(&c_wvs[i][j]);
            mat_alloc(&c_wvs[i][j], pwv->size1, pwv->size2);
        }
    }

#ifdef DEBUG
    printf("cband: %lu x %lu x %lu",
            cband[0][0].size1, cband[0][0].size2, cband[0][0].size3);
    for (i = 0; i < cband[0][0].size1; i++)
    {
        
        printf("\n %02d:", i);
        for (j = 0; j < cband[0][0].size2; j++)
        {
            for (k = 0; k < cband[0][0].size3; k++)
                printf(" %.4f", mat3d_get(&cband[0][0], i, j, k));
            printf("\n    ");
        }
    }
    printf("\n");
#endif

#ifdef DEBUG_HDF
    int dbg_stat = 0;
    hid_t dbg_id = create_hdf5("TgWVSdebug1.h5");
    if (dbg_id == FAIL)
        ERROR("Could not create TgWVSdebug1.h5", NONE);
    dbg_stat  = hdf5write_mat3d(dbg_id, "x", &x);
    dbg_stat += hdf5write_mat3d(dbg_id, "y", &y);
    dbg_stat += hdf5write_mat3d(dbg_id, "z", &z);
    dbg_stat += hdf5write_mat(dbg_id, "senszen", senszen);
    dbg_stat += hdf5write_mat(dbg_id, "solzen", solzen);
    dbg_stat += hdf5write_mat(dbg_id, "emis_aster", &emis_aster_loc);
    dbg_stat += hdf5write_mat3d(dbg_id, "lrad", lrad);
    dbg_stat += hdf5write_mat(dbg_id, "pwv", pwv);

    close_hdf5(dbg_id);
    if (dbg_stat < 0)
        WARNING("%d error(s) writing datasets to TgWVSdebug1.h5",
            -dbg_stat);
#endif

// Setup the input and output datasets for the tri-linear interpolation.
    int n_dsets = MAX_CHANNELS * N_COEFFS_PER_BAND;
    double *insets[n_dsets];
    double *outsets[n_dsets];
    int n = 0;
    for (i = 0; i < MAX_CHANNELS; i++)
    {
        for (k = 0; k < N_COEFFS_PER_BAND; k++)
        {
            insets[n] = cband[i][k].vals;
            outsets[n] = c_wvs[i][k].vals;
            n++;
        }
    }

    // Interpolate all coefficient datasets
    multi_interp3(x.vals, y.vals, z.vals, nv_uniq, ne_uniq, nc_uniq,
        senszen->vals, emis_aster_loc.vals, pwv->vals, pwv->size1 * pwv->size2,
        (const double**)insets, outsets, n_dsets);

    // Clean up memory not needed any longer.
    vec_clear(&csort);
    vec_clear(&cw);
    vec_clear(&cwi);
    vec_clear(&emis);
    vec_clear(&emisi);
    vec_clear(&esort);
    vec_clear(&views);
    vec_clear(&viewsi);
    vec_clear(&vsort);
    mat_clear(&emis_aster_loc);
    mat3d_clear(&x);
    mat3d_clear(&y);
    mat3d_clear(&z);
    for (i = 0; i < MAX_CHANNELS; i++)
    {
        mat_clear(&coeff[i]);
        mat_clear(&coeff_in[i]);
        for (j = 0; j < N_COEFFS_PER_BAND; j++)
        {
            mat3d_clear(&cband[i][j]);
        }
    }

    // Compute Tg
    // Interpolate brightness temperatures from lookup table
    Mat3d tb;
    mat3d_init(&tb);
    mat3d_alloc(&tb, n_channels, lrad->size2, lrad->size3);
    double interp_val, lrad_val;
    debug_lste_lib = 1;
    for (i = 0; i < n_channels; i++)
    {
        for (j = 0; j < tb.size2; j++)
        {
            for (k = 0; k < tb.size3; k++)
            {
                // This statement was resultsing in a segmentation fault.
                //mat3d_set(&tb, i, j, k, 
                //    interp1d_npts(lut[band[i]+1], lut[0], 
                //        mat3d_get(lrad, i, j, k), nlut_lines));
                lrad_val =  mat3d_get(lrad, i, j, k);
                interp_val =
                    interp1d_npts(lut[band[i]+1], lut[0], lrad_val, nlut_lines);
                mat3d_set(&tb, i, j, k, interp_val);
                debug_lste_lib = 0;
            }
        }
    }

#ifdef DEBUG_HDF
    dbg_stat = 0;
    dbg_id = create_hdf5("TgWVSdebug2.h5");
    if (dbg_id == FAIL)
        ERROR("Could not create TgWVSdebug2.h5", NONE);
    char dbgdset[256];
    int dbgband, dbgcoef;
    dbg_stat = 0;
    for (dbgband = 0; dbgband < MAX_CHANNELS; dbgband++)
    {
        for (dbgcoef = 0; dbgcoef < N_COEFFS_PER_BAND; dbgcoef++)
        {
            // TODO Move this output to above before clearing cband arrays.
            //sprintf(dbgdset, "cband_%d_%d", dbgband+1, dbgcoef+1);
            //dbg_stat += hdf5write_mat3d(dbg_id, dbgdset, &cband[dbgband][dbgcoef]);
            sprintf(dbgdset, "c_wvs_%d_%d", dbgband+1, dbgcoef+1);
            dbg_stat += hdf5write_mat(dbg_id, dbgdset, &c_wvs[dbgband][dbgcoef]);
        }
    }
    close_hdf5(dbg_id);
    if (dbg_stat < 0)
        WARNING("%d error(s) writing datasets to TgWVSdebug2.h5",
            -dbg_stat);
#endif

#ifdef DEBUG
    describe_mat3d("TB", &tb);
#endif

#ifdef DEBUG_PIXEL
    printf("    tb(%d,%d)=[%g, %g, %g]\n",
            DEBUG_LINE, DEBUG_PIXEL,
            mat3d_get(&tb, 0, DEBUG_LINE, DEBUG_PIXEL),
            mat3d_get(&tb, 1, DEBUG_LINE, DEBUG_PIXEL),
            mat3d_get(&tb, 2, DEBUG_LINE, DEBUG_PIXEL));

#endif

    double tg0, tgsum;
    double p, p2;
    int cs;
    int b;
    for (i = 0; i < n_channels; i++)
    {
        b = band[i];
        for (j = 0; j < tg->size2; j++)
        {
            for (k = 0; k < tg->size3; k++)
            {
                p = mat_get(pwv, j, k);
                p2 = p * p;
                tg0 = 
                    mat_get(&c_wvs[b][0], j, k) 
                  + mat_get(&c_wvs[b][1], j, k) * p
                  + mat_get(&c_wvs[b][2], j, k) * p2;
                tgsum = 0.0;
                cs = 0;
                //Accumulate data from each channel
                for (n = 0; n < n_channels; n++)
                {
                    tgsum +=
                      (   mat_get(&c_wvs[b][cs+3], j, k)
                        + mat_get(&c_wvs[b][cs+4], j, k) * p
                        + mat_get(&c_wvs[b][cs+5], j, k) * p2
                       ) * mat3d_get(&tb, n, j, k);
                    cs += 3;
                }
                mat3d_set(tg, i, j, k, tg0 + tgsum);
#ifdef DEBUG_PIXEL
                if (i == 0 && j == DEBUG_LINE && k == DEBUG_PIXEL)
                {
                    printf("    c_wvs array: %d x %d  matrix: %lu x %lu\n",
                            MAX_CHANNELS, N_COEFFS_PER_BAND, c_wvs[0][0].size1,
                            c_wvs[0][0].size2);
                    printf("    Tg[%d,%d,%d]: PWV=%.4f PWV^2=%.4f\n",
                            i, j, k, p, p2);
                    printf("        tg0 = %.4f + %.4f * %.4f + %.4f * %.4f "
                           "= %.4f\n",
                           mat_get(&c_wvs[i][0], j, k),
                           mat_get(&c_wvs[i][1], j, k), p,
                           mat_get(&c_wvs[i][2], j, k), p2,
                           tg0);
                    int idb;
                    double dup_tgsum = 0.0;
                    double dbgamt;
                    int dup_cs = 0;
                    int bdbg;
                    for (idb = 0; idb < n_channels; idb++)
                    {
                        bdbg = band[idb];
                        dbgamt = ( mat_get(&c_wvs[bdbg][dup_cs+3], j, k)
                                 + mat_get(&c_wvs[bdbg][dup_cs+4], j, k) * p
                                 + mat_get(&c_wvs[bdbg][dup_cs+5], j, k) * p2
                                 ) * mat3d_get(&tb, idb, j, k);
                        dup_tgsum += dbgamt;
                        printf("        tgsum += (%.4f + %.4f * %.4f + "
                               "%.4f * %.4f) * %.4f\n"
                               "              += %.4f = %.4f\n",
                               mat_get(&c_wvs[bdbg][dup_cs+3], j, k),
                               mat_get(&c_wvs[bdbg][dup_cs+4], j, k), p,
                               mat_get(&c_wvs[bdbg][dup_cs+5], j, k), p2,
                               mat3d_get(&tb, idb, j, k), dbgamt, dup_tgsum
                               );
                        dup_cs += 3;
                    }
                    printf("    Tg = tg0(%.4f) + tgsum(%.4f) = %.4f\n", 
                            tg0, dup_tgsum, tg0 + dup_tgsum);
                }
#endif
            }
        }
    }

#ifdef DEBUG_HDF
    dbg_id = create_hdf5("TgWVSdebug3.h5");
    if (dbg_id == FAIL)
        ERROR("Could not create TgWVSdebug3.h5", NONE);
    dbg_stat  = hdf5write_mat3d(dbg_id, "TB", &tb);
    dbg_stat += hdf5write_mat3d(dbg_id, "Tg", tg);
    if (dbg_stat < 0)
        WARNING("%d errors out of 2 result datasets written to TgWVSdebug3.h5",
            -dbg_stat);
    close_hdf5(dbg_id);
#endif

    // Release memory used for c_wvs.
    for (i = 0; i < MAX_CHANNELS; i++)
    {
        for (j = 0; j < N_COEFFS_PER_BAND; j++)
        {
            mat_clear(&c_wvs[i][j]);
        }
    }
    mat3d_clear(&tb);
    return 0;
}
