// describe.c
//
// Describe a dataset in a file.

#include "describe.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lste_lib.h"

void describe_data(double *pd, size_t nd)
{
    double maxneg = -1E10;
    double minneg = 0.0;
    double maxpos = 0.0;
    double minpos = 1E10;
    double totpos = 0;
    double totneg = 0;
    int nnan = 0;
    int nneg = 0;
    int npos = 0;
    int nzero = 0;
    int ninfn = 0;
    int ninfp = 0;
    size_t i;
    double d;

    for (i = 0; i < nd; i++)
    {
        d = *pd++;
        if (is_nan(d))
            nnan++;
        else if (isinf(d))
        {
            if (d < 0.0)
                ninfn++;
            else
                ninfp++;
        }
        else if (d < 0.0)
        {
            nneg++;
            totneg += d;
            if (d < minneg) minneg = d;
            if (d > maxneg) maxneg = d;
        }
        else if (d > 0.0)
        {
            npos++;
            totpos += d;
            if (d < minpos) minpos = d;
            if (d > maxpos) maxpos = d;
        }
        else
            nzero++;
    }
    if (ninfn > 0)
        printf("        has %d negative infinity values\n", ninfn);
    if (ninfp > 0)
        printf("        has %d positive infinity values\n", ninfp);
    printf("        has %d nan values\n", nnan);
    printf("        has %d zero values\n", nzero);
    printf("        has %d positive values", npos);
    if (npos > 0) printf(" from %.4f to %.4f, mean=%.4f", minpos, maxpos, totpos/npos);
    printf("\n        has %d negative values", nneg);
    if (nneg > 0) printf(" from %.4f to %.4f, mean=%.4f", minneg, maxneg, totneg/nneg);
    if (npos + nneg > 0)
        printf("\n        total mean=%.4f\n", (totpos + totneg) / (npos + nneg));
    printf("\n");
    fflush(stdout);
}

void describe_uint8(uint8 *pu, size_t nd)
{
    uint8 d;
    int i;
    int count[256];
    for (i = 0; i < 256; i++) 
        count[i] = 0;

    for (i = 0; i < nd; i++)
    {
        d = *pu++;
        ++count[d];
    }

    printf("    contains the following values:\n");
    for (i = 0; i < 256; i++)
    {
        if (count[i] > 0)
            printf("        [%03d] %d times\n", i, count[i]);
    }
    printf("\n");
    fflush(stdout);
}

void describe_mat_u16_bitset(const char *name, MatUint16 *mat)
{
    uint16 d;
    uint16 *pu = mat->vals;
    int i;
    int count[65536];
    for (i = 0; i < 65536; i++) 
        count[i] = 0;

    printf("%s: %lu x %lu\n", name, mat->size1, mat->size2);

    for (i = 0; i < mat->size1 * mat->size2; i++)
    {
        d = *pu++;
        ++count[d];
    }

    printf("    contains the following values:\n");
    for (i = 0; i < 65536; i++)
    {
        if (count[i] > 0)
            printf("        [%03d] %d times\n", i, count[i]);
    }
    printf("\n");
    fflush(stdout);
}

void describe_uint16(uint16 *pd, size_t nd)
{
    unsigned max = 0;
    unsigned min = 0xFFFF;
    unsigned totpos = 0;
    unsigned npos = 0;
    unsigned nzero = 0;
    size_t i;
    uint16 d;

    for (i = 0; i < nd; i++)
    {
        d = *pd++;
        if (d > 0)
        {
            npos++;
            totpos += d;
            if (d < min) min = d;
            if (d > max) max = d;
        }
        else
            nzero++;
    }
    printf("        has %d zero values\n", nzero);
    printf("        has %d positive values", npos);
    if (npos > 0) printf(" from %u to %u, mean=%u", min, max, 
            (unsigned)((double)totpos/(double)npos));
    printf("\n");
    fflush(stdout);
}

void describe_data_2d(const char *name, double *pd, size_t nrows, size_t ncols)
{
    printf("%s: %lu x %lu\n", name, nrows, ncols);
    describe_data(pd, nrows * ncols);
}

void describe_matrix(const char *name, Matrix *mat)
{
    describe_data_2d(name, mat->vals, mat->size1, mat->size2);
}

void describe_mat_u8(const char *name, MatUint8 *mat)
{
    printf("%s: %lu x %lu\n", name, mat->size1, mat->size2);
    describe_uint8(mat->vals, mat->size1 * mat->size2);
}

void describe_mat_u16(const char *name, MatUint16 *mat)
{
    printf("%s: %lu x %lu\n", name, mat->size1, mat->size2);
    describe_uint16(mat->vals, mat->size1 * mat->size2);
}

void describe_mat3d(const char *name, Mat3d *mat)
{
    printf("%s: %lu x %lu x %lu\n", 
            name, mat->size1, mat->size2, mat->size3);
    describe_data(mat->vals, mat->size1 * mat->size2 * mat->size3);
}

void describe_data_3d(const char *name, double *pd[], 
    size_t nlevels, size_t nrows, size_t ncols)
{
    Mat3d flatmat;
    mat3d_init(&flatmat);
    mat3d_alloc(&flatmat, nlevels, nrows, ncols);
    size_t nvals = nrows * ncols;
    size_t lvlsize = nvals * sizeof(double);
    size_t z;
    for (z = 0; z < nlevels; z++)
    {
        memcpy(&flatmat.vals[z * nvals], pd[z], lvlsize);
    }
    describe_mat3d(name, &flatmat);
    mat3d_clear(&flatmat);
}

void describe_vector(const char *name, Vector *vec)
{
    printf("%s: %lu\n", name, vec->size);
    describe_data(vec->vals, vec->size);
}

