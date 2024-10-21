// maptogrid.c
#include "maptogrid.h"
#include <assert.h>
#include <math.h>
#include "bool.h"

//#define DEBUG
#ifdef DEBUG
int debug_flag = 1;
#else
int debug_flag = -1;
#endif
int debug_on = 0;
int debug_mapnearest = 0;
//#define DEBUG_STEPS

/**

    !Description:  
    @brief Given a latitude and longitude, this function determines the 
    line, sample, subline factor and subsample factor of the point 
    in the DAO geographic grid.

    @author P. Fisher, NASA

    !Input Parameters:
    @param[in] lat latitude of point in question
    @param[in] lon longitude of point in question
    @param[in] latres latitude resolution of the input data
    @param[in] lonres longitude resolution of the input data
    @param[in] samps columns in the input data file

    !Output Parameters:
    @param[out] iline line number containing the point in question
    @param[out] isamp sample number containing the point in question
    @param[out] rowfactor vertical distance of the point in question from the 
        center of iline to the center of the next southern line.
    @param[out] colfactor horizontal distance of the point in question from the 
        center of isamp to the center of the next eastern pixel.

    !Input/Output Parameters:

    !Returns:
    0	failure (lat and/or lon is out of range)
    1	success
    @return 0=failure, 1=success
    
    !Revision History:
    Revision 1.1, August 1997, P. Fisher
    1.  Added absolute value call into rowfactor/colfactor calculation
    (bug fix).

    Original version 1997, P. Fisher (fisher@ltpmail.gsfc.nasa.gov).

    !Team-unique Header:

    This software developed by the MODIS Land Science Team for 
    the National Aeronautics and Space Administration, Goddard 
    Space Flight Center, under contract NAS5-96062

    !References and Credits

    Paul Fisher
    Science Systems and Applications Inc.
    NASA's Goddard Space Flight Center Code 614.5
    Greenbelt MD, 20771
    fisher@ltpmail.gsfc.nasa.gov, (301) 286-2980

    !Design Notes
    !END***************************************************************************
    */
int get_cart_linesamp(double lat, double lon, double latres, int samps, 
    double lonres, int *iline,
    int *isamp, double *rowfactor, double *colfactor, 
    double startLat, double startLon)
{
    *iline = *isamp = (int) 0;
    *colfactor = *rowfactor = (double)0.0; /* just in case */

    /* check inputs */
    if(lat>90.0F||lat<startLat||lon>180.0F||lon<startLon)return(0);

    //-- *iline = (int)(((90.0F - lat)/latres)+0.5F);
    *iline = (int)(((lat - startLat)/latres)+0.5F);
    //--  *rowfactor = (double)((90.0F - ((double)(*iline)*latres)) - lat)/latres;
    *rowfactor = ((double)*iline * latres + startLat - lat)/latres;
    if (*rowfactor < 0.0 && *iline > 0){
    	(*iline)--;
    	*rowfactor += 1.0;
    }

    /* special wrap around case */
    if(lon>(180.0F-(lonres*0.5F)))lon=-180.0F-(180.0F-lon);

    *isamp = (int)(((lon - startLon)/lonres) + 0.5F);
    *colfactor = ((double)*isamp * lonres + startLon - lon)/lonres;
    if (*colfactor < 0.0){
    	(*isamp)--;
    	/* special wrap around case */
    	if((*isamp)<0)*isamp=(samps-1);
    	*colfactor += 1.0;
    }

    return(1);
}


/**
!C*****************************************************************************
!Description:  
!@brief bilinearly interpolates a point in question through its 
 isamp, iline, row and column subsample distances 
 (supplied through a call to get_cart_linesamp()) for data 
 in the DAO geographic grid.

@author P. Fisher, NASA
@author Robert Freepaqrtner, JPL

!Input Parameters:
@param iline line number containing the point in question
@param isamp sample number containing the point in question
@param samps number of columns in the input DAO data file
@param rowfactor vertical distance of the point in question from the 
       center of iline to the center of the next southern line.
@param colfactor horizontal distance of the point in question from the 
       center of isamp to the center of the next eastern pixel.
@param indata pointer to input dataset

!Output Parameters: none

!Input/Output Parameters: none

!Returns: interpolated data value
@return interpolated data value

!Revision History:
Original version 1997, P. Fisher (fisher@ltpmail.gsfc.nasa.gov).
Modification 2015 for MOD21 PGE116, Robert Freepartner, JPL

!Team-unique Header:

This software developed by the MODIS Land Science Team for 
the National Aeronautics and Space Administration, Goddard 
Space Flight Center, under contract NAS5-96062

!References and Credits

Paul Fisher
Science Systems and Applications Inc.
NASA's Goddard Space Flight Center Code 614.5
Greenbelt MD, 20771
fisher@ltpmail.gsfc.nasa.gov, (301) 286-2980

!Design Notes
!END****************************************************************************/
double get_cart_interp(int iline, int isamp, int samps, double rowfactor, 
    double colfactor, double *indata)
{
    double v1, v2, v3, v4;
    int i1, i2, i3, i4;
    double interp;

    /* obtain four corners */
    v1= 0.0;
    v2= 0.0;
    v3= 0.0;
    v4= 0.0;
    i1= iline * samps + isamp;
    i2= (iline + 1) * samps + isamp;
    i3= (iline + 1) * samps + isamp+1;
    i4= iline * samps + isamp+1;
    if (i1 >= 0) v1 = indata[i1];
    if (i2 >= 0) v2 = indata[i2];
    if (i3 >= 0) v3 = indata[i3];
    if (i4 >= 0) v4 = indata[i4];

    if (debug_on)
    {
        printf("get_cart_interp(%d, %d, %d, %.4f, %.4f, %p)\n",
                iline, isamp, samps, rowfactor, colfactor, indata);
        printf("    1: [%d, %d] = %.4f factor = %.4f\n",iline, isamp, v1,
               (1.0F-rowfactor)*(1.0F-colfactor)); 
        printf("    2: [%d, %d] = %.4f factor = %.4f\n",iline+1, isamp, v2,
                rowfactor*(1.0F-colfactor)); 
        printf("    3: [%d, %d] = %.4f factor = %.4f\n",iline+1, isamp+1, v3,
                rowfactor*colfactor); 
        printf("    4: [%d, %d] = %.4f factor = %.4f\n",iline, isamp+1, v4,
                colfactor*(1.0F-rowfactor)); 
    }

    /*  v1 = anc->data[iline * samps + isamp];
    v2 = anc->data[(iline + 1) * samps + isamp];
    v3 = anc->data[(iline + 1) * samps + isamp+1];
    v4 = anc->data[iline * samps + isamp+1];
    */
    interp=(1.0F-rowfactor)*(1.0F-colfactor)*v1 +  
    	rowfactor*(1.0F-colfactor)*v2 +
    	rowfactor*colfactor*v3 + 
    	colfactor*(1.0F-rowfactor)*v4;

    if (debug_on)
    {
        printf("    interp: %.4f\n\n", interp);
    }

    return(interp);

}

void meshgrid(Matrix *latgrid, Matrix *longrid, Vector *lats, Vector *lons)
{
    assert(latgrid != NULL);
    assert(longrid != NULL);
    assert(lats != NULL);
    assert(lons != NULL);
    int llrow, llcol;
    mat_init(latgrid);
    mat_init(longrid);
    mat_alloc(latgrid, lats->size, lons->size);
    mat_alloc(longrid, lats->size, lons->size);
    for (llrow = 0; llrow < lats->size; llrow++)
    {
        for (llcol = 0; llcol < lons->size; llcol++)
        {
            mat_set(latgrid, llrow, llcol, vec_get(lats, llrow));
            mat_set(longrid, llrow, llcol, vec_get(lons, llcol));
        }
    }
}

void meshgrat(Matrix *latgrid, Matrix *longrid, double minLat, double minLon, 
        double maxLat, double maxLon, double latRes, double lonRes)
{
    assert(latgrid != NULL);
    assert(longrid != NULL);
    assert(minLat < maxLat);
    assert(minLon < maxLon);
    assert(latRes > 0.0);
    assert(lonRes > 0.0);
    double startLat = floor(minLat);
    double startLon = floor(minLon);
    double endLat = ceil(maxLat);
    double endLon = ceil(maxLon);
    int nLats = (endLat - startLat)/latRes;
    int nLons = (endLon - startLon)/lonRes;
    int ilat, ilon;
    mat_alloc(latgrid, nLats, nLons);
    mat_alloc(longrid, nLats, nLons);
    double *dlat = latgrid->vals;
    double *dlon = longrid->vals;
    double latv = startLat;
    double lonv = startLon;
    for (ilat = 0; ilat < nLats; ilat++)
    {
        lonv = startLon;
        for (ilon = 0; ilon < nLons; ilon++)
        {
            *dlat++ = latv;
            *dlon++ = lonv;
            lonv += lonRes;
        }
        latv += latRes;
    }
}

/**
 * This function is based on map_ancillary_data() 
 *
 *     !References and Credits 
 *  This software was developed by:
 *
 *    MODIS Science Data Support Team for the National Aeronautics and Space 
 *    Administration, Goddard Space Flight Center, under contract NAS5-32373.
 *
 *    Developer:
 *    Jim Ray
 *    Science Systems and Applications, Inc.
 *    MODIS Land Science Team
 *    NASA Goddard Space Flight Center
 *    Greenbelt, Maryland 20771
 *    jim@ltdri.org
 *
 * The modifications from map_ancillary_data to maptogrid for MOD21 PGE116 were
 * developed by:
 *   [BF] Robert Freepartner, JPL, October, 2015
 *   December 2015, Adapted general version to mapworldtogrid [BF]
 *
 */
int maptogrid(GridMapData *gmd, 
    double startLat, double startLon, double latres, double lonres)
{
    int    col, row;
    double colfactor, rowfactor;
    int    id;
    int    iline, isamp;
    int    qindex;
    int    stat;

#ifdef DEBUG
    printf("maptogrid(gmd, %.4f, %.4f, %.4f, %.4f)\n",
        startLat, startLon, latres, lonres);
    printf("    gmd.lat = %.4f, %.4f, %.4f ...\n",
        gmd->lat[0], gmd->lat[1], gmd->lat[2]);
    printf("    gmd.lon = %.4f, %.4f, %.4f ...\n",
        gmd->lon[0], gmd->lon[1], gmd->lon[2]);
    printf("    gmd.nlines = %d\n", gmd->nlines);
    printf("    gmd.ncols  = %d\n", gmd->ncols);
    printf("    gmd.inrows = %d\n", gmd->inrows);
    printf("    gmd.incols = %d\n", gmd->incols);
    printf("    gmd.nsets  = %d\n", gmd->nsets);
#endif

    if (gmd->nsets < 1) return -1;
    
    qindex = 0;		       
    for (row=0; row < gmd->nlines; row++) 
    {
    	for (col=0; col < gmd->ncols; col++)
        {
            if (qindex == debug_flag)
            {
                debug_on = 1;
                printf("Calling get_cart_linesamp(%.4f, %.4f, %.4f, %d, %.4f, ...);\n",
                       gmd->lat[qindex], gmd->lon[qindex], latres, gmd->incols, lonres);
            }
            else
                debug_on = 0;

            // This step does the processing that is common to all the datasets
            // being interpolated. 
            stat = get_cart_linesamp(gmd->lat[qindex], gmd->lon[qindex], 
                latres, gmd->incols, lonres, &iline, &isamp, &rowfactor, &colfactor,
                startLat, startLon);
            if (!stat)
            {
                // If failure, set NaNs.
                for (id = 0; id < gmd->nsets; id++)
                {
                    gmd->dsetout[id][qindex] = nan;
                }
            }
            else
            {
                if (debug_on)
                {
                    printf("Return values from get_cart_linesamp: iline=%d, "
                            "isamp=%d, rowfactor=%.4f, colfactor=%.4f\n",
                            iline, isamp, rowfactor, colfactor);
                }

                // Interpolate all the datasets
                for (id = 0; id < gmd->nsets; id++)
                {
                    // This step is done for each dataset to interpolate the
                    // data for the current grid poisiont: row, col.
                    gmd->dsetout[id][qindex] =
                        get_cart_interp(iline, isamp, gmd->incols, rowfactor, 
                        colfactor, gmd->dsetin[id]);
                }
            }
            qindex++;
        }
    }
    return 0;
}

int mapworldtogrid(GridMapData *gmd)
{
    double latres = 180.0 / (gmd->inrows - 1);
    double lonres = 360.0 / gmd->incols;
    return maptogrid(gmd, -90.0, -180.0, latres, lonres);
}

// Uses a radial search that is modified for bowtie effects. If the 
// row before or after the target point contains a NaN in the sample
// data, the row is selected farther until the data does not
// contain a NaN or is outside the grid.
int mapnearest(GridMapData *gmd, double *latsamp, double *lonsamp)
{
    int    col, row;
    int    id;
    int    samprow, sampcol;
    int    qindex;
    double curlat, curlon;
    double findlat, findlon;
    double distsq;
    double dlat, dlon;
    double nearestdistsq;
    int    nearestrow, nearestcol;
    int    startrow, startcol;
    double closestdlat, closestdlon;
    int    closestrow, closestcol;
    int    ioff, index;
    bool   rowgood;
    const double fill = -9999.0;

    if (gmd->nsets < 1) return -1;

#ifdef DEBUG
    debug_flag = 2 * gmd->ncols;

    printf("mapnearest: inrows=%d incols=%d outrows=%d outcols=%d O(%ld)\n",
            gmd->inrows, gmd->incols, gmd->nlines, gmd->ncols,
            (long)gmd->inrows * (long)gmd->incols * 
            (long)gmd->nlines * (long)gmd->ncols);
#endif

    // Create matrix to flag "already checked"
    MatUint8 checked;
    mat_uint8_init(&checked);
    mat_uint8_alloc(&checked, gmd->inrows, gmd->incols);

    // Loop through the query points
    qindex = -1;		       
    for (row=0; row < gmd->nlines; row++) 
    {
    	for (col=0; col < gmd->ncols; col++)
        {
            qindex++;
            // For each query lat/lon, find closest point in the sample set.
            if (debug_mapnearest != 0 || qindex == debug_flag)
            {
                debug_on = 1;
                printf("  row=%d col=%d qindex=%d lat=%.3f lon=%.3f\n",
                       row, col, qindex, gmd->lat[qindex], gmd->lon[qindex]);
            }
            else
                debug_on = 0;

            // Find a starting point for the radial search.
            findlat = gmd->lat[qindex];
            findlon = gmd->lon[qindex];
            // First find the closest lat in column 0.
            closestdlat = 99.0; // Bigger than any lat diff
            closestrow = gmd->inrows / 2; // Middle row
            for (startrow = 0; startrow < gmd->inrows; startrow++)
            {
                curlat = latsamp[startrow * gmd->incols];
                if (is_nan(curlat) || dequal(curlat, fill))
                    continue;
                dlat = fabs(findlat - curlat);
                if (dlat < closestdlat)
                {
                    closestdlat = dlat;
                    closestrow = startrow;
                }
            }
            startrow = closestrow;

            // Next, find the closest lon in the startrow.
            closestdlon = 999.0;
            closestcol = gmd->incols / 2;
            for (startcol = 0; startcol < gmd->incols; startcol++)
            {
                curlon = lonsamp[startrow * gmd->incols + startcol];
                if (is_nan(curlon) || dequal(curlon, fill))
                    continue;
                dlon = fabs(findlon - curlon);
                if (dlon < closestdlon)
                {
                    closestdlon = dlon;
                    closestcol = startcol;
                }
            }
            startcol = closestcol;

            if (debug_on)
            {
                printf("  startrow=%d startcol=%d lat=%.3f lon=%.3f "
                       "findlat=%.3f findlon=%.3f\n",
                       startrow, startcol, latsamp[startrow * gmd->incols + startcol],
                       lonsamp[startrow * gmd->incols + startcol], findlat, findlon);
            }

            int pass = 0;
            bool foundnearer = false;

            // Find the closest point in the sample data.
            dlat = findlat - latsamp[startrow * gmd->incols + startcol];
            dlon = findlon - lonsamp[startrow * gmd->incols + startcol];
            double startdistsq = dlat * dlat + dlon * dlon;
            // In case of samples full of NaNs, could be Nan.
            if (is_nan(startdistsq))
                startdistsq = 90.0 * 90.0 + 180.0 * 180.0;
            nearestrow = startrow;
            nearestcol = startcol;
            nearestdistsq = startdistsq;
            if (debug_on)
            {
                printf("  dlat=%.6f dlon=%.6f startdistsq=%.6f\n", dlat, dlon, startdistsq);
            }
            //Clear the checked map.
            memset(checked.vals, 0, checked.size1 * checked.size2 * sizeof(uint8));

            do // Repeat until a nearer location is not found.
            {
                startrow = nearestrow;
                startcol = nearestcol;
                startdistsq = nearestdistsq;
                foundnearer = false;
                int prerow = startrow;
                int postrow = startrow;
                // Adjust prerow/portrow to skip over rows with NaNs.
                rowgood = false;
                while (!rowgood)
                {
                    prerow--;
                    if (prerow < 0)
                        break;
                    index = prerow * gmd->incols + startcol;
                    curlat = latsamp[index];
                    if (is_nan(curlat) || dequal(curlat, fill))
                        continue;
                    if (startcol > 0)
                    {
                        curlat = latsamp[index - 1];
                        if (is_nan(curlat) || dequal(curlat, fill))
                            continue;
                    }
                    if (startcol < gmd->incols - 1)
                    {
                        curlat = latsamp[index + 1];
                        if (is_nan(curlat) || dequal(curlat, fill))
                            continue;
                    }
                    rowgood = true;
                }

                rowgood = false;
                while (!rowgood)
                {
                    postrow++;
                    if (postrow >= gmd->inrows)
                        break;
                    index = postrow * gmd->incols + startcol;
                    curlat = latsamp[index];
                    if (is_nan(curlat) || dequal(curlat, fill))
                        continue;
                    if (startcol > 0)
                    {
                        curlat = latsamp[index - 1];
                        if (is_nan(curlat) || dequal(curlat, fill))
                            continue;
                    }
                    if (startcol < gmd->incols - 1)
                    {
                        curlat = latsamp[index + 1];
                        if (is_nan(curlat) || dequal(curlat, fill))
                            continue;
                    }
                    rowgood = true;
                }

                if (debug_on)
                {
                    printf("  prerow=%d postrow=%d midrow=%d midcol=%d\n", 
                            prerow, postrow, startrow, startcol);
                }

                // Loop over up to eight locations around the start point.
                for (samprow = prerow; samprow <= postrow; samprow++)
                {
                    if (samprow < 0) continue;
                    if (samprow >= gmd->inrows) break;
                    ioff = samprow * gmd->incols; // Offset to start of row
                    for (sampcol = startcol - 1; sampcol <= startcol + 1; sampcol++)
                    {
                        if (sampcol < 0) continue;
                        if (sampcol >= gmd->incols) break;
                        // Use checked map to prevent repeat distance calcs.
                        if (mat_uint8_get(&checked, samprow, sampcol) != 0)
                            continue;
                        mat_uint8_set(&checked, samprow, sampcol, 1);

                        // Skip the position that has already been calculated.
                        if (samprow == startrow && sampcol == startcol) continue;

                        curlat = latsamp[ioff + sampcol];
                        if (is_nan(curlat) || dequal(curlat, fill))
                            continue;
                        curlon = lonsamp[ioff + sampcol];
                        if (is_nan(curlon) || dequal(curlon, fill))
                            continue;
                        // Get distance and compare to nearest-so-far
                        dlat = findlat - curlat;
                        dlon = findlon - curlon;
                        distsq = dlat * dlat + dlon * dlon;
                        if (debug_on)
                        {
                            printf("    r=%d c=%d dlat=%.6f dlon=%.6f distsq=%.6f\n",
                                samprow, sampcol, dlat, dlon, distsq);
                        }
                        if (distsq < nearestdistsq)
                        {
                            // Found a closer point
                            nearestrow = samprow;
                            nearestcol = sampcol;
                            nearestdistsq = distsq;
                            foundnearer = true;
                        }
                    }
                }
                if (debug_on)
                {
                    printf("  pass=%d nearestrow=%d nearestcol=%d nearest dsq=%.6f"
                           " foundnearer=%s\n",
                            ++pass, nearestrow, nearestcol, nearestdistsq,
                            (foundnearer ? "yes" : "no"));
                }
            } while (foundnearer);

            if (debug_on)
            {
                printf("    query:  [%d,%d] lat=%.3f lon=%.3f\n"
                       "    nearest:[%d,%d] lat=%.3f lon=%.3f\n",
                       row, col, findlat, findlon, 
                       nearestrow, nearestcol, 
                       latsamp[nearestrow * gmd->incols + nearestcol],
                       lonsamp[nearestrow * gmd->incols + nearestcol]);
            }

            // Found nearestrow, nearestcol to data point in sample.
            // Map each dataset with data from nearest point.
            double vnearest;
            int dir = (nearestrow > gmd->inrows - 12) ? -1 : 1;
            for (id = 0; id < gmd->nsets; id++)
            {
                // This step is done for each dataset to set the nearest
                // sample value for the current query grid position.
                // The loop is used to avoid using NaN values.
                int r = nearestrow;
                do
                {
                    vnearest = gmd->dsetin[id][r * gmd->incols + nearestcol];
                    r += dir;
                } while (is_nan(vnearest) && r >= 0 && r < gmd->inrows);
                gmd->dsetout[id][qindex] = vnearest;

                if (debug_on)
                {
                    printf("    Dataset[%d]=%.4f from [%d,%d]\n", 
                           id, vnearest, r, nearestcol);
                }
#ifdef DEBUG_STEPS
                if (id > 0) continue;
                r -= dir;
                dlat = findlat - latsamp[r * gmd->incols + nearestcol];
                dlon = findlon - lonsamp[r * gmd->incols + nearestcol];
                distsq = dlat * dlat + dlon * dlon;
                printf("r=%4d c=%4d nr=%4d nc=%4d d2=%.4f\n",
                    row, col, r, nearestcol, distsq);
#endif
            }
        }
    }
    mat_uint8_clear(&checked);
    return 0;
}
