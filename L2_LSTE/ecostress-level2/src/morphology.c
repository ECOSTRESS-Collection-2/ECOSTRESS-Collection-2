#include "morphology.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

int enable_check_8 = 1;
const uint64_t all_ones = 0x0101010101010101;

void construct_Kernel(Kernel *kernel)
{
    kernel->noffs = 0;
    kernel->offsets = NULL;
}

void destruct_Kernel(Kernel *kernel)
{
    free(kernel->offsets);
    free(kernel);
}

Kernel *disk_kernel(int radius)
{
    Kernel *kernel = (Kernel*)malloc(sizeof(Kernel));
    assert(kernel != NULL);
    construct_Kernel(kernel);

    int w = radius * 2 + 1;
    kernel->noffs = 0;
    // Allocate the maximum possible number of points.
    // TODO In c++ use vector and push_back.
    kernel->offsets = (OffsetCoordinate*)malloc(w*w*sizeof(OffsetCoordinate));
    assert(kernel->offsets != NULL);

    int midr = radius + 1;
    int midc = radius + 1;
    int r, c;
    int dr, dc;
    double dksq = (double)radius * (double)radius;
    double drsq, dcsq;
    double dsq;
    for (r = 0; r <= midr + radius; r++)
    {
        for (c = 0; c <= midc + radius; c++)
        {
            dr = r - midr;
            drsq = (double)dr * (double)dr;
            dc = c - midc;
            dcsq = (double)dc * (double)dc;
            dsq = drsq + dcsq;
            if (dsq <= dksq)
            {
                kernel->offsets[kernel->noffs].dy = dr;
                kernel->offsets[kernel->noffs].dx = dc;
                kernel->noffs++;
            }
        }
    }
    return kernel;
}

Kernel *square_kernel(int radius)
{
    Kernel *kernel = (Kernel*)malloc(sizeof(Kernel));
    assert(kernel != NULL);
    construct_Kernel(kernel);

    int w = radius * 2 + 1;
    kernel->noffs = 0;
    // Allocate the maximum possible number of points.
    // TODO In c++ use vector and push_back.
    kernel->offsets = (OffsetCoordinate*)malloc(w*w*sizeof(OffsetCoordinate));
    assert(kernel->offsets != NULL);

    int r, c;
    for (r = -radius; r <= radius; r++)
    {
        for (c = -radius; c <= radius; c++)
        {
            kernel->offsets[kernel->noffs].dy = r;
            kernel->offsets[kernel->noffs].dx = c;
            kernel->noffs++;
        }
    }
    return kernel;
}

mpoint_t *dilation(const mpoint_t *points, int nrows, int ncols, const Kernel* kernel)
{
    int r,c,dr,dc,ioff,ipoint,idx;
    int set_flag;
    ipoint = 0;
    // Create the output array
    mpoint_t *outp = (mpoint_t*)malloc(nrows*ncols*sizeof(mpoint_t));
    assert(outp != NULL);
    uint64_t *p8 = (uint64_t*)points; // 8 bytes at once pointer
    

    for (r = 0; r < nrows; r++)
    {
        for (c = 0; c < ncols; c++, ipoint++)
        {
            outp[ipoint] = points[ipoint]; // Init output
            // If the current pixel is a foreground pixel, leave it set.
            if (points[ipoint] != 0)
                continue;
            // Check the points in the kernel ofsets list.
            set_flag = 0;
            for (ioff = 0; ioff < kernel->noffs; ioff++)
            {
                // Convert relative coordinate from kernel to points index.
                dr = r + kernel->offsets[ioff].dy;
                dc = c + kernel->offsets[ioff].dx;
                // Skip points that are off the edge.
                if (dr < 0 || dr >= nrows)
                    continue;
                if (dc < 0 || dc >= ncols)
                    continue;
                // Performance: can check 8 at a time?
                int check_8 = 0;
                if (dc % 8 == 0)
                {
                    if (enable_check_8)
                    {
                        int ioff8 = ioff + 7;
                        if (ioff8 < kernel->noffs && r + kernel->offsets[ioff8].dy == dr)
                        {
                            ioff += 7;
                            check_8 = 1;
                        }
                    }
                }
                if (check_8)
                {
                    idx = dr * (ncols/8) + dc/8;
                    if (p8[idx] != 0)
                    {
                       set_flag = 1;
                       break;
                    }
                }
                else
                {
                    idx = dr * ncols + dc;
                    // Checking for at least one foreground pixel.
                    if (points[idx] != 0)
                    {
                        set_flag = 1;
                        break;
                    }
                }
            }
            // If any foreground pixel was found in kernel, make
            // the current point foreground.
            if (set_flag)
                outp[ipoint] = 1;
        }
    }
    return outp;
}

mpoint_t *erosion(const mpoint_t *points, int nrows, int ncols, const Kernel* kernel)
{
    int r,c,dr,dc,ioff,ipoint,idx;
    int set_flag;
    ipoint = 0;
    // Create the output array
    mpoint_t *outp = (mpoint_t*)malloc(nrows*ncols*sizeof(mpoint_t));
    assert(outp != NULL);
    uint64_t *p8 = (uint64_t*)points; // 8 bytes at once pointer

    for (r = 0; r < nrows; r++)
    {
        for (c = 0; c < ncols; c++, ipoint++)
        {
            outp[ipoint] = points[ipoint]; // Init output
            // If the current pixel is a background pixel, leave it clear.
            if (points[ipoint] == 0)
                continue;
            // Check the points in the kernel ofsets list.
            set_flag = 0;
            for (ioff = 0; ioff < kernel->noffs; ioff++)
            {
                // Convert relative coordinate from kernel to points index.
                dr = r + kernel->offsets[ioff].dy;
                dc = c + kernel->offsets[ioff].dx;
                // Skip points that are off the edge.
                if (dr < 0 || dr >= nrows)
                    continue;
                if (dc < 0 || dc >= ncols)
                    continue;
                // Performance: can check 8 at a time?
                int check_8 = 0;
                if (dc % 8 == 0)
                {
                    if (enable_check_8)
                    {
                        int ioff8 = ioff + 7;
                        if (ioff8 < kernel->noffs && r + kernel->offsets[ioff8].dy == dr)
                        {
                            ioff += 7;
                            check_8 = 1;
                        }
                    }
                }
                if (check_8)
                {
                    idx = dr * (ncols/8) + dc/8;
                    if (p8[idx] != all_ones)
                    {
                       set_flag = 1;
                       break;
                    }
                }
               else
               {
                   idx = dr * ncols + dc;
                   // Checking for at least one background pixel.
                   if (points[idx] == 0)
                   {
                      set_flag = 1;
                       break;
                   }
               }
            }
            // If any background pixel was found in kernel, make
            // the current point background.
            if (set_flag)
                outp[ipoint] = 0;
        }
    }
    return outp;
}

mpoint_t *imclose(const mpoint_t *points, int nrows, int ncols, const Kernel* kernel)
{
    mpoint_t *dilated = dilation(points, nrows, ncols, kernel);
    mpoint_t *closed =  erosion(dilated, nrows, ncols, kernel);
    free(dilated);
    return closed;
}

