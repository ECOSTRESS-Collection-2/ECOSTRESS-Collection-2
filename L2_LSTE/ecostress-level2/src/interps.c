/** interps.c -- see interps.h for documentation
 * Includes C code from the article
 * "Tri-linear Interpolation"
 * by Steve Hill, sah@ukc.ac.uk
 * in "Graphics Gems IV", Academic Press, 1994
 */
#include "interps.h"
#include "lste_lib.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//#define DEBUG 
//#define DEBUG_PIXEL 2032300
#ifdef DEBUG
extern int debug_binsearch;
#endif

#ifndef _NaN_
#define _NaN_ -32768.0
#endif

#ifndef _epsilon_
#define _epsilon_ 1e-9
#endif

/* LERP from "Graphics Gems IV", Academic Press, 1994
 * linear interpolation from l (when a=0) to h (when a=1)
 * (equal to (a*h)+((1-a)*l) */
#define LERP(a,l,h) ((l)+(((h)-(l))*(a)))

double dval2(const double *d, int x, int y, int ysize)
{
    return d[x*ysize + y];
}

double dval3(const double *d, int x, int y, int z, int xsize, int ysize)
{
    return d[(z*xsize + x)*ysize + y];
}

int requal(double v1, double v2)
{
    return (fabs(v1-v2) < _epsilon_);
}

int findIx(double qtarg, double *samp, int nsamps)
{
    int lowix = 0;
    if (samp[nsamps-1] < samp[0])
    {
        // Samples are decreasing
        if (qtarg >= samp[1])
        {
            lowix = 0; // Use the first sample point
        }
        else if (qtarg <= samp[nsamps-2])
        {
            lowix = nsamps - 2; // Use the next-to-the-last sample point
        }
        else
        {
            lowix = dec_binsearch(samp, qtarg, 0, nsamps); // Search
        }
    }
    else
    {
        // Samples are increasing
        if (qtarg <= samp[1])
        {
            lowix = 0;
        }
        else if (qtarg >= samp[nsamps-2])
        {
            lowix = nsamps - 2;
        }
        else
        {
            lowix = inc_binsearch(samp, qtarg, 0, nsamps);
        }
    }
    return lowix;
}

// interp2 sets up a single query set for multi_interp2 so that there is 
// only one implementation of the interpolation algorithm.
void interp2(
    const double *X, const double *Y, double *V,
    int nsampsx, int nsampsy,
    const double *Xq, const double *Yq, double *Vq,
    int nqpoints)
{
    double *insets[1];
    insets[0] = V;
    double *outsets[1];
    outsets[0] = Vq;
    
    multi_interp2(X, Y, nsampsx, nsampsy, Xq, Yq, 
        nqpoints, insets, outsets, 1);
}

void multi_interp2(
    const double *X, const double *Y, int nsampsx, int nsampsy,
    const double *Xq, const double *Yq, int nqpoints,
    double *insets[], double *outsets[], int nsets)
{
    // Index variables
    int i;
    int ix, iy;
    int iq, iset;
#ifdef DEBUG
    int tests = 0;
    int skips = 0;
#endif

    // Surrounding point values and interp results
    double d00, d01, d10, d11, dy0, dy1, dxy;
    // Axis slopes
    double fx, fy;

#ifdef DEBUG
    int debug_on = 0;
    tests = 12;
    skips = nqpoints / (tests + 1);
    printf("multi_interp2(X=%p, Y=%p, nsampsx=%d, nsampsy=%d, "
           "Xq=%p, Yq=%p, nqpoints=%d)\n",
           X, Y, nsampsx, nsampsy, Xq, Yq, nqpoints);
#endif

    // Create linear coordinate vectors
    double *xsamp = malloc(nsampsx * sizeof(double));
    for (i = 0; i < nsampsx; i++)
    {
        xsamp[i] = dval2(X, i, 0, nsampsy);
    }
    double *ysamp = malloc(nsampsy * sizeof(double));
    for (i = 0; i < nsampsy; i++)
    {
        ysamp[i] = dval2(Y, 0, i, nsampsy);
    }

#ifdef DEBUG
    printf("X: [%d]", nsampsx);
    for (i = 0; i < nsampsx; i++)
        printf(" %.4f", xsamp[i]);
    printf("\nY: [%d]", nsampsy);
    for (i = 0; i < nsampsy; i++)
        printf(" %.4f", ysamp[i]);
    printf("\n");
#endif

    // Loop over query points
    for (iq = 0; iq < nqpoints; iq++)
    {
#ifdef DEBUG
        if (iq % skips == 0 && tests-- > 0)
        {
            debug_on = 1;
            debug_binsearch = 1;
        }
#endif
        ix = findIx(Xq[iq], xsamp, nsampsx);
        iy = findIx(Yq[iq], ysamp, nsampsy);
        fx = (Xq[iq] - xsamp[ix])/(xsamp[ix+1]-xsamp[ix]);
        fy = (Yq[iq] - ysamp[iy])/(ysamp[iy+1]-ysamp[iy]);

        for (iset = 0; iset < nsets; iset++)
        {
            // Code based on interp3
            d00 = dval2(insets[iset], ix, iy, nsampsy);
            d10 = dval2(insets[iset], ix, iy+1, nsampsy);
            d01 = dval2(insets[iset], ix+1, iy, nsampsy);
            d11 = dval2(insets[iset], ix+1, iy+1, nsampsy);

            if (requal(d00, _NaN_) ||
                requal(d01, _NaN_) ||
                requal(d10, _NaN_) ||
                requal(d11, _NaN_))
            {
                dxy = _NaN_;
            }
            else
            {
                dy0 = LERP(fy, d00, d10);
                dy1 = LERP(fy, d01, d11);
                dxy = LERP(fx, dy0, dy1);
            }
            outsets[iset][iq] = dxy;
#ifdef DEBUG
            if (debug_on)
            {
                printf("iq=%d, Xq=%.4f, Yq=%.4f\n", iq, Xq[iq], Yq[iq]);
                printf("multi_interp2 ix=%d iy=%d\n", ix, iy);
                printf("Xix=%.4f Yiy=%.4f\n", xsamp[ix], ysamp[iy]);
                printf("Xix+1=%.4f Yiy+1=%.4f\n", xsamp[ix+1], ysamp[iy+1]);
                printf("fx=%.4f, fy=%.4f\n", fx, fy);
                printf("d00=%.4f d01=%.4f d10=%.4f d11=%.4f\n", 
                        d00, d01, d10, d11);
                printf("dy0=%.4f dy1=%.4f dxy=%.4f\n", dy0, dy1, dxy);
                debug_on = 0;
                debug_binsearch = 0;
            }
#endif
        }

    }
    free(xsamp);
    free(ysamp);
}

// interp3 sets up a single query set for multi_interp3 so that there is 
// only one implementation of the interpolation algorithm.
void interp3(
    const double *X, const double *Y, const double *Z, const double *V,
    int nsampsx, int nsampsy, int nsampsz, 
    const double *Xq, const double *Yq, const double *Zq, double *Vq,
    int nqpoints)
{
    const double *insets[1];
    insets[0] = V;
    double *outsets[1];
    outsets[0] = Vq;
    
    multi_interp3(X, Y, Z, nsampsx, nsampsy, nsampsz, Xq, Yq, Zq, 
        nqpoints, insets, outsets, 1);
}

void multi_interp3(
    const double *X, const double *Y, const double *Z, 
    int nsampsx, int nsampsy, int nsampsz, 
    const double *Xq, const double *Yq, const double *Zq, int nqpoints,
    const double *insets[], double *outsets[], int nsets)
{
    // Index variables
    int i;
    int ix, iy, iz;
    int iq, iset;
#ifdef DEBUG
    int tests = 0;
    int skips = 0;
#endif

    // Surrounding point values 
    double d000, d001, d010, d011,
           d100, d101, d110, d111,
           dy00, dy01, dy10, dy11,
           dxy0, dxy1, dxyz;
    double fx, fy, fz;

#ifdef DEBUG
    int debug_on = 0;
    tests = 12;
    skips = nqpoints / (tests + 1);
    printf("multi_interp3(X=%p, Y=%p, Z=%p, nsampsx=%d, nsampsy=%d, nsampsz=%d, "
           "Xq=%p, Yq=%p, Zq=%p, nqpoints=%d, nsets=%d)\n",
           X, Y, Z, nsampsx, nsampsy, nsampsz, Xq, Yq, Zq, nqpoints, nsets);
#endif

    // Create linear coordinate vectors
    double *xsamp = malloc(nsampsx * sizeof(double));
    for (i = 0; i < nsampsx; i++)
    {
        xsamp[i] = dval3(X, i, 0, 0, nsampsx, nsampsy);
    }
    double *ysamp = malloc(nsampsy * sizeof(double));
    for (i = 0; i < nsampsy; i++)
    {
        ysamp[i] = dval3(Y, 0, i, 0, nsampsx, nsampsy);
    }
    double *zsamp = malloc(nsampsz * sizeof(double));
    for (i = 0; i < nsampsz; i++)
    {
        zsamp[i] = dval3(Z, 0, 0, i, nsampsx, nsampsy);
    }

#ifdef DEBUG
    printf("X: [%d]", nsampsx);
    for (i = 0; i < nsampsx; i++)
        printf(" %.4f", xsamp[i]);
    printf("\nY: [%d]", nsampsy);
    for (i = 0; i < nsampsy; i++)
        printf(" %.4f", ysamp[i]);
    printf("\nZ: [%d]", nsampsz);
    for (i = 0; i < nsampsz; i++)
        printf(" %.4f", zsamp[i]);
    printf("\n");
#endif

    // Loop over query points
    for (iq = 0; iq < nqpoints; iq++)
    {
#ifdef DEBUG
        if (iq > 7 && iq < 16) debug_on = 1;
        if (iq % skips == 0 && tests-- > 0)
            debug_on = 1;
#endif
        if (requal(Xq[iq], _NaN_) ||
            requal(Yq[iq], _NaN_) ||
            requal(Zq[iq], _NaN_))
        {
            for (iset = 0; iset < nsets; iset++)
                outsets[iset][iq] = _NaN_;
            continue;
        }

        ix = findIx(Xq[iq], xsamp, nsampsx);
        iy = findIx(Yq[iq], ysamp, nsampsy);
        iz = findIx(Zq[iq], zsamp, nsampsz);

        fx = (Xq[iq] - xsamp[ix])/(xsamp[ix+1]-xsamp[ix]);
        fy = (Yq[iq] - ysamp[iy])/(ysamp[iy+1]-ysamp[iy]);
        fz = (Zq[iq] - zsamp[iz])/(zsamp[iz+1]-zsamp[iz]);

        for (iset = 0; iset < nsets; iset++)
        {
            d000 = dval3(insets[iset], ix,  iy,  iz,  nsampsx, nsampsy);
            d001 = dval3(insets[iset], ix,  iy+1,iz,  nsampsx, nsampsy);
            d010 = dval3(insets[iset], ix+1,iy,  iz,  nsampsx, nsampsy);
            d011 = dval3(insets[iset], ix+1,iy+1,iz,  nsampsx, nsampsy);
            d100 = dval3(insets[iset], ix,  iy,  iz+1,nsampsx, nsampsy);
            d101 = dval3(insets[iset], ix,  iy+1,iz+1,nsampsx, nsampsy);
            d110 = dval3(insets[iset], ix+1,iy,  iz+1,nsampsx, nsampsy);
            d111 = dval3(insets[iset], ix+1,iy+1,iz+1,nsampsx, nsampsy);

            if (requal(d000, _NaN_) ||
                requal(d001, _NaN_) ||
                requal(d010, _NaN_) ||
                requal(d011, _NaN_) ||
                requal(d100, _NaN_) ||
                requal(d101, _NaN_) ||
                requal(d110, _NaN_) ||
                requal(d111, _NaN_))
            {
                dxyz = _NaN_;
            }
            else
            {
                dy00 = LERP(fy, d000, d001);
                dy01 = LERP(fy, d010, d011);
                dy10 = LERP(fy, d100, d101);
                dy11 = LERP(fy, d110, d111);
                dxy0 = LERP(fx, dy00, dy01);
                dxy1 = LERP(fx, dy10, dy11);
                dxyz = LERP(fz, dxy0, dxy1);
            }

            outsets[iset][iq] = dxyz;
#ifdef DEBUG
#ifdef DEBUG_PIXEL
            if (iq == DEBUG_PIXEL) debug_on = 1;
#endif
            if (debug_on > 0)
            {
                printf("iq=%d, Xq=%.4f, Yq=%.4f, Zq=%.4f\n", iq, Xq[iq], Yq[iq], Zq[iq]);
                printf("interp3 ix=%d iy=%d iz=%d\n", ix, iy, iz);
                printf("Xix=%.4f Yiy=%.4f Ziz=%.4f\n", 
                        xsamp[ix], ysamp[iy], zsamp[iz]);
                printf("Xix+1=%.4f Yiy+1=%.4f Ziz+1=%.4f\n", 
                        xsamp[ix+1], ysamp[iy+1], zsamp[iz+1]);
                printf("fx=%.4f, fy=%.4f, fz=%.4f\n", fx, fy, fz);
                printf("d000=%10.4f d010=%10.4f d100=%10.4f d110=%10.4f\n", 
                        d000,       d010,       d100,       d110);
                printf("d001=%10.4f d011=%10.4f d101=%10.4f d111=%10.4f\n", 
                        d001,       d011,       d101,       d111);
                printf("dy00=%10.4f dy10=%10.4f\ndy01=%10.4f dy11=%10.4f\n", 
                        dy00,       dy10,        dy01,       dy11);
                printf("dxy0=%10.4f dxy1=%10.4f\ndxyz=%10.4f\n\n", 
                        dxy0, dxy1, dxyz);
                debug_on = 0;
            }
#endif
        }

    }
    free(xsamp);
    free(ysamp);
    free(zsamp);
}

