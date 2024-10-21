// lstelib.c
// See lstlib.h for credits and documentation
//
// Treat values less than zero as not-a-number
//#define NEG_IS_NAN

#include <assert.h>
#include "lste_lib.h"

int debug_lste_lib = 0;

bool is_nan(double val)
{
    if (isnan(val)) return true;
    return fabs(val - _NaN_) < epsilon;
}

bool dequal(double v1, double v2)
{
    return fabs(v1-v2) < epsilon;
}

int dec_binsearch(double *d, double v, int imin, int imax)
{
    if (imin >= imax)
        return imin;
    else
    {
        int imid = imin + (imax - imin)/2;
        if (imid == imin) imid = imax;
        if (d[imid] < v)
            return dec_binsearch(d, v, imin, imid - 1);
        else if (d[imid] > v)
            //return dec_binsearch(d, v, imid + 1, imax);
            return dec_binsearch(d, v, imid, imax);
        else
            return imid;
    }
}
int debug_binsearch = 0;
int inc_binsearch(double *d, double v, int imin, int imax)
{
    if (debug_binsearch > 0)
        printf("inc_binsearch(d, v=%.4f, imin=%d, imax=%d): ", 
                   v, imin, imax);
    if (imin >= imax)
    {
        if (debug_binsearch > 0)
            printf("returns imin=%d\n", imin);
        return imin;
    }
    else
    {
        int imid = imin + (imax - imin)/2;
        if (imid == imin) imid = imax;
        if (debug_binsearch > 0)
            printf("d[imid=%d] = %.4f\n", 
                   imid, d[imid]);
        if (d[imid] > v)
        {
            return inc_binsearch(d, v, imin, imid - 1);
        }
        else if (d[imid] < v)
        {
            //return inc_binsearch(d, v, imid + 1, imax);
            return inc_binsearch(d, v, imid, imax);
        }
        else
        {
            if (debug_binsearch > 0)
                printf("returns imid=%d\n", imid);
            return imid;
        }
    }
}

int collection_index(const char* MOD_filename)
{
    const char *c6 = (const char*)strstr(MOD_filename, ".006.");
    if (c6 != NULL) return 1;
    return 0;
}

double distance(double lat1, double lon1, double lat2, double lon2)
{
    // haversine formula
    const double emeanrad = 6371229;
    double lat1rad = lat1 * PI / 180.0;
    double lat2rad = lat2 * PI / 180.0;
    double lon1rad = lon1 * PI / 180.0;
    double lon2rad = lon2 * PI / 180.0;
    double latdiff = lat2rad - lat1rad;
    double londiff = lon2rad - lon1rad;
    double sinhalflatdiff = sin(latdiff/2.0);
    double sinhalflondiff = sin(londiff/2.0);
    double coslat1 = cos(lat1rad);
    double coslat2 = cos(lat2rad);
    double a = sinhalflatdiff * sinhalflatdiff +
                coslat1 * coslat2 * sinhalflondiff * sinhalflondiff;
    double c = 2.0 * atan2(sqrt(a), sqrt(1.0-a));
#ifdef DEBUG
    printf("Lat/Lon %.2f %.2f to %.2f %.2f = %.2f meters\n",
        lat1, lon1, lat2, lon2, emeanrad * c);
#endif
    return emeanrad * c;
}

double mean(double *v, int n)
{
  int i;
  double mval;

  mval = 0;
  for(i=0; i<n; i++) mval += v[i];
  return(mval/n);
}

double stddev(double *v, int n) 
{
  int i;
  double sval, mval;
 
  sval = 0;
  mval = mean(v,n);
  for(i=0; i<n; i++) sval += (v[i]-mval)*(v[i]-mval);
  return(sqrt(sval/n));
}

double covvar(double *v1, double *v2, int n)
{
  int i;
  double mval1, mval2, cval;

  mval1 = mean(v1,n);
  mval2 = mean(v2,n);
  cval = 0;
  for(i=0; i<n; i++) cval += (v1[i]-mval1)*(v2[i]-mval2);
  return(cval/n);
} 

double corrcoef(double *v1, double *v2, int n) 
{
  return(covvar(v1,v2,n)/(stddev(v1,n)*stddev(v2,n)));
}

int getnearest(int nx_in, int ny_in, int nx_out, int ny_out, int ix, int iy)
{
  double hx, hy;
  int ii,jj;

  hx = (double)nx_out/(double)nx_in;    
  hy = (double)ny_out/(double)ny_in;    

  ii = (int)(ix*hx);
  if(ii > (nx_out-1)) ii = nx_out-1;
  jj = (int)(iy*hy);
  if(jj > (ny_out-1)) jj = ny_out-1;
  return(jj*nx_out+ii);
}

int getnearestInv(int nx_in, int ny_in, int nx_out, int ny_out, int ix, int iy)
{
  double hx, hy;
  int ii,jj;

  hx = (double)nx_in/(double)nx_out;
  hy = (double)ny_in/(double)ny_out;

  ii = (int)(ix*hx);
  if(ii > (nx_in-1)) ii = nx_in-1;
  jj = (int)(iy*hy);
  if(jj > (ny_in-1)) jj = ny_in-1;
  return(jj*nx_in+ii);
}

int get_nnans(double *in, int nx, int ny)
{
    int n = nx * ny;;
    int nnan = 0;
    while (n-- > 0) 
    {
        if (is_nan(*in++)) nnan++;
    }
    return nnan;
}

bool has_nan(double *in, int nx, int ny)
{
    int n = nx * ny;;
    while (n-- > 0) 
    {
        if (is_nan(*in++)) return true;
    }
    return false;
}

void imresize(double* in, int nx_in, int ny_in, double* out, int nx_out, int ny_out)
{
  int i,j;
  double hx = (double)nx_in/(double)nx_out;
  double hy = (double)ny_in/(double)ny_out;
  int ix, iy;
  for(j=0; j < ny_out; j++) 
  {
    for(i=0; i < nx_out; i++) 
    {
      ix = (int)(i * hx);
      if (ix >= nx_in) ix = nx_in - 1;
      iy = (int)(j * hy);
      if (iy >= ny_in) iy = ny_in - 1;
      out[j * nx_out + i] = in[iy * nx_in + ix];
    } 
  }
  return;
}

void imresize_u8(uint8* in, int nx_in, int ny_in, uint8* out, int nx_out, int ny_out)
{
  int i,j;
  double hx = (double)nx_in/(double)nx_out;
  double hy = (double)ny_in/(double)ny_out;
  int ix, iy;
  for(j=0; j < ny_out; j++) 
  {
    for(i=0; i < nx_out; i++) 
    {
      ix = (int)(i * hx);
      if (ix >= nx_in) ix = nx_in - 1;
      iy = (int)(j * hy);
      if (iy >= ny_in) iy = ny_in - 1;
      out[j * nx_out + i] = in[iy * nx_in + ix];
    } 
  }
  return;
}

double interp1d(double *x, double *y, double xx) {
  double yy, m, b;

  if(xx < x[1]) {
    m = (y[1]-y[0])/(x[1]-x[0]);
  } else {
    m = (y[2]-y[1])/(x[2]-x[1]);
  }
  b = y[1]-m*x[1];
  yy = m*xx+b;
/*
fprintf(stderr,"xx = %lf, slope = %lf, intercept = %lf\n",xx,m,b);
*/
  return(yy);
}

double interp1d_npts(double *x, double *y, double xx, int npts) {
  double m; // slope
  double b; // intercept

if (debug_lste_lib) 
{
    printf("interp1d_npts(%p, %p, %g, %d)\n",
        x, y, xx, npts);
    printf("    x=%g, %g, ... %g, %g\n",
        x[0], x[1], x[npts-2], x[npts-1]);
    printf("    y=%g, %g, ... %g, %g\n",
        y[0], y[1], y[npts-2], y[npts-1]);
}

  if (npts < 2)
  {
	  // Can't interpolate without at least 2 points.
	  return nan;
  }

  int lowix = 0;
  int incdec = 1; // 1 = inclining points, -1 = declining 
  if (x[npts-1] < x[0]) incdec = -1;
  if (incdec < 0)
  {
      if (xx >= x[1])
      {
          lowix = 0;
      }
      else if (xx <= x[npts-2])
      {
          lowix = npts - 2;
      }
      else
      {
          lowix = dec_binsearch(x, xx, 0, npts);
      }
  }
  else
  {
      if (xx <= x[1])
      {
          lowix = 0;
      }
      else if (xx >= x[npts-2])
      {
          lowix = npts - 2;
      }
      else
      {
          lowix = inc_binsearch(x, xx, 0, npts);
      }
  }
  int highix = lowix + 1;
  
  m = (y[highix]-y[lowix]) / (x[highix]-x[lowix]);
  b = y[highix] - m * x[highix];
  return m * xx + b;
}

void interp_id(double *in, int nx, int ny, double r2, double e)
{
  int i,j,iloop;
  int nnan;
  double *vc, *vcc, *D, *D2, *V, *wV, sumV, sumwV;
  int *xg, *yg, *x, *y, nxyg, nxy, index, nD2r2;
  int *xD2r2, *yD2r2;

  x = (int *)malloc(nx*ny*sizeof(int));
  y = (int *)malloc(nx*ny*sizeof(int));
  xg = (int *)malloc(nx*ny*sizeof(int));
  yg = (int *)malloc(nx*ny*sizeof(int));

  nnan = get_nnans(in, nx, ny);

  iloop = 0;
  while(nnan > 0) {
    vc = (double *)malloc((nx*ny-nnan)*sizeof(double));
    j = 0;
    for(i=0; i<nx*ny; i++) {
      if(in[i] < 0) continue;
      // Check for infinite loop situation
      assert(j <= (nx*ny-nnan-1));
      vc[j++] = in[i];
    }
    assert(j == (nx*ny-nnan)); //check for final count problem

    nxyg = 0; nxy = 0;
    for(i=0; i<ny; i++) {
      for(j=0; j<nx; j++) {
        index = i*nx+j;
        if(in[index] > 0) {xg[nxyg] = j; yg[nxyg] = i; nxyg++;}
        if(in[index] < 0) {x[nxy] = j; y[nxy] = i; nxy++;}
    } }
    D = (double *)malloc(nxyg*sizeof(double));
    D2 = (double *)malloc(nxyg*sizeof(double));
    vcc = (double *)malloc(nxyg*sizeof(double));
    xD2r2 = (int *)malloc(nxyg*sizeof(int));
    yD2r2 = (int *)malloc(nxyg*sizeof(int));

    for(j=0; j<nxy; j++) {
      nD2r2 = 0;
      for(i=0; i<nxyg; i++) {
        D[i] = sqrt((x[j]-xg[i])*(x[j]-xg[i])+(y[j]-yg[i])*(y[j]-yg[i]));
        if(D[i] < r2) {
          D2[nD2r2] = D[i];
          xD2r2[nD2r2] = xg[i];
          yD2r2[nD2r2] = yg[i];
          nD2r2++;
        }
      }
      if(nD2r2 == 0) continue;
      V = (double *)malloc(nD2r2*sizeof(double));
      wV = (double *)malloc(nD2r2*sizeof(double));
      for(i=0; i<nD2r2; i++) {
        wV[i] = pow(D2[i],e);
        index = yD2r2[i]*nx+xD2r2[i];
        vcc[i] = in[index];
        V[i] = vcc[i]*wV[i];
      }
      sumV = 0; sumwV = 0;
      for(i=0; i<nD2r2; i++) {
        sumV += V[i]; sumwV += wV[i];
      }
      index = y[j]*nx+x[j];
            in[index] = sumV/sumwV;
      free(V); free(wV);
    }
    free(vc); free(D); free(D2); free(vcc); free(xD2r2); free(yD2r2);
    nnan = get_nnans(in, nx, ny);
    iloop++;
    assert(iloop <= nx*ny); // ERROR("interp_id while endless loop","main");
  }
#ifdef TEST
  fprintf(stderr,"iloop = %d\n",iloop);
#endif
  free(x); free(y); free(xg); free(yg);
}

void interp_id_eneg2(double *in, int nx, int ny, double r2)
{
    int i,j,iloop;
    int dx;
    int dy;
    double dsquared;
    double rsquared = r2 * r2;
    double *vcc, *D2, *V, *wV, sumV, sumwV;
    int *xg, *yg, *x, *y, nxyg, nxy, index, nD2r2;
    int *xD2r2, *yD2r2;

    x = (int *)malloc(nx*ny*sizeof(int));
    y = (int *)malloc(nx*ny*sizeof(int));
    xg = (int *)malloc(nx*ny*sizeof(int));
    yg = (int *)malloc(nx*ny*sizeof(int));

    iloop = 0;
    // Step 1
    // Save the index values in x,y for Nans and xg,yg for good values.
    while(has_nan(in, nx, ny)) {
        nxyg = 0; nxy = 0;
        for(i=0; i<ny; i++) {
            for(j=0; j<nx; j++) {
                index = i*nx+j;
                if(in[index] > 0) {xg[nxyg] = j; yg[nxyg] = i; nxyg++;}
                if(in[index] < 0) {x[nxy] = j; y[nxy] = i; nxy++;}
            }
        }

        D2 = (double *)malloc(nxyg*sizeof(double));  // distance squared
        vcc = (double *)malloc(nxyg*sizeof(double)); // unweighted good values
        xD2r2 = (int *)malloc(nxyg*sizeof(int));     // x of non-NaN in range
        yD2r2 = (int *)malloc(nxyg*sizeof(int));     // y of non-NaN in range

        // For each NaN in the dataset:
        for(j=0; j<nxy; j++) {
            // Save all non-NaN that are in range (r2) of the NaN value.
            nD2r2 = 0;
            for(i=0; i<nxyg; i++) {
                dy = y[j]-yg[i];
                if (fabs(dy) >= r2)
                    continue;
                dx = x[j]-xg[i];
                if (fabs(dx) >= r2)
                    continue;
                dsquared = dx * dx + dy * dy;
                if(dsquared < rsquared) {
                    D2[nD2r2] = dsquared;
                    xD2r2[nD2r2] = xg[i];
                    yD2r2[nD2r2] = yg[i];
                    nD2r2++;
                }
            }
            if(nD2r2 == 0) continue;

            // Perform value weighting by inverse distance
            V = (double *)malloc(nD2r2*sizeof(double));
            wV = (double *)malloc(nD2r2*sizeof(double));
            for(i=0; i<nD2r2; i++) {
                wV[i] = 1.0 / D2[i]; // i.e., pow(distance, -2.0);
                index = yD2r2[i]*nx+xD2r2[i];
                vcc[i] = in[index];
                V[i] = vcc[i]*wV[i];
            }

            // Interpolate the weighted value
            sumV = 0; sumwV = 0;
            for(i=0; i<nD2r2; i++) {
                sumV += V[i]; 
                sumwV += wV[i];
            }

            // Replace the NaN
            index = y[j]*nx+x[j];
            in[index] = sumV/sumwV;
            free(V); 
            free(wV);
        }
        free(D2); free(vcc); free(xD2r2); free(yD2r2);
        iloop++;
        assert(iloop <= nx*ny); //ERROR("interp_id_eneg2 while endless loop","main");
    }
    free(x); free(y); free(xg); free(yg);
}

void multiset_interp_id(double *ds[], int nds, int nx, int ny, double r2)
{
    int i,j,n,iloop;
    int dx;
    int dy;
    double dsquared;
    double rsquared = r2 * r2;
    double *D2, *V[nds], *wV, sumV, sumwV;
    int *xg, *yg, *x, *y, nxyg, nxy, index, nD2r2;
    int *xD2r2, *yD2r2;

    x = (int *)malloc(nx*ny*sizeof(int));
    y = (int *)malloc(nx*ny*sizeof(int));
    xg = (int *)malloc(nx*ny*sizeof(int));
    yg = (int *)malloc(nx*ny*sizeof(int));

    iloop = 0;
    // Step 1
    // Save the index values in x,y for Nans and xg,yg for good values.
    while(has_nan(ds[0], nx, ny)) {
        //fprintf(stderr, "multiset_interp_id iloop=%d\n", iloop);
        nxyg = 0; nxy = 0;
        for(i=0; i<ny; i++) {
            for(j=0; j<nx; j++) {
                index = i*nx+j;
                if(ds[0][index] > 0) {xg[nxyg] = j; yg[nxyg] = i; nxyg++;}
                if(ds[0][index] < 0) {x[nxy] = j; y[nxy] = i; nxy++;}
            }
        }
        //fprintf(stderr, "  good=%d, nan=%d\n", nxyg, nxy);

        D2 = (double *)malloc(nxyg*sizeof(double));  // distance squared
        xD2r2 = (int *)malloc(nxyg*sizeof(int));     // x of non-NaN in range
        yD2r2 = (int *)malloc(nxyg*sizeof(int));     // y of non-NaN in range
        //fprintf(stderr, "  mallocs: D2 = %p, xD2r2=%p, yD2r2=%p\n", D2, xD2r2, yD2r2);

        // For each NaN in the dataset:
        for(j=0; j<nxy; j++) {
            // Save all non-NaN that are in range (r2) of the NaN value.
            nD2r2 = 0;
            for(i=0; i<nxyg; i++) {
                dy = y[j]-yg[i];
                if (dy >= r2)
                    continue;
                dx = x[j]-xg[i];
                if (dx >= r2)
                    continue;
                dsquared = dx * dx + dy * dy;
                if(dsquared < rsquared) {
                    D2[nD2r2] = dsquared;
                    xD2r2[nD2r2] = xg[i];
                    yD2r2[nD2r2] = yg[i];
                    nD2r2++;
                }
            }
            //fprintf(stderr, "  nD2r2=%d\n", nD2r2);
            if(nD2r2 == 0) continue;

            // Perform value weighting by inverse distance
            for (n = 0; n < nds; n++)
            {
                V[n] = (double *)malloc(nD2r2*sizeof(double));
            }
            wV = (double *)malloc(nD2r2*sizeof(double));
            for(i=0; i<nD2r2; i++) {
                wV[i] = 1.0 / D2[i]; // i.e., pow(distance, -2.0);
                index = yD2r2[i]*nx+xD2r2[i];
                for (n = 0; n < nds; n++)
                {
                    V[n][i] = ds[n][index]*wV[i];
                }
            }

            for (n = 0; n < nds; n++)
            {
               // Interpolate the weighted value
                sumV = 0; sumwV = 0;
                for(i=0; i<nD2r2; i++)
                {
                    sumV += V[n][i];
                    sumwV += wV[i];
                }

                // Replace the NaN
                index = y[j]*nx+x[j];
                ds[n][index] = sumV/sumwV;
            }

            for (n = 0; n < nds; n++)
            {
                free(V[n]); 
            }
            free(wV);
        }
        free(D2); free(xD2r2); free(yD2r2);
        iloop++;
        assert(iloop <= nx*ny); // ERROR("interp_id_eneg2 while endless loop","main");
    }
    free(x); free(y); free(xg); free(yg);
}

void minmax2d(double *d, int nx, int ny, double fill, double *min, double *max)
{
    int i = 0;
    *min = d[0];
    *max = d[0];
    for (i = 1; i < nx * ny; i++)
    {
	    if (is_nan(*min))
		{
		    *min = d[i];
			*max = d[i];
			continue;
		}
        if (is_nan(d[i]) || dequal(d[i], fill))
            continue;
        if (dequal(*min, fill))
        {
            *min = d[i];
	        *max = d[i];
            continue;
        }
        if (d[i] <  *min)
            *min = d[i];
        if (d[i] > *max)
            *max = d[i];
    }
    
}

double round2(double x, double y) {return round (x/y)*y;}

/* This version came from MOD_LST_main1.c */
void smooth2a(double* in, int nx, int ny, int dx, int dy)
{
  double *local[ny];
  double sum, bufferx[nx], buffery[ny];
  int i,j,k,ix,iy,m;
  /// @TBD Wy was ij set but not used? Check older versions.
  //int ij; 

  for(i=0; i<ny; i++) local[i] = in+i*nx;

// do averaging along each row:
  for(j=0; j<ny; j++) {
    //ij = j*nx;
    for(i=0; i<nx; i++) {
//
//    if(fabs(local[j][i]-nan) < epsilon) {
//      bufferx[i] = nan;
//      continue;
//    }
      m = 0; sum = 0;
      for(k=-dx; k<=dx; k++) {
        ix = i+k;
        if(ix < 0) continue;
        if(ix > (nx-1)) continue;
        if(is_nan(local[j][ix])) continue;
        sum += local[j][ix]; m++;
      }
      if(m > 0) bufferx[i] = sum/(double)m;
      else bufferx[i] = nan;
    }
    for(i=0; i<nx; i++) local[j][i] = bufferx[i];
  }

// do averaging along each column: 
  for(i=0; i<nx; i++) {
    for(j=0; j<ny; j++) {
      //ij = j*nx+i;
//
//    if(fabs(local[j][i]-nan) < epsilon) {
//      buffery[j] = nan;
//      continue;
//    }

      m = 0; sum = 0;
      for(k=-dy; k<=dy; k++) {
        iy = j+k;
        if(iy < 0) continue;
        if(iy > (ny-1)) continue;
        if(is_nan(local[iy][i])) continue;
        sum += local[iy][i]; m++;
      }
      if(m > 0) buffery[j] = sum/(double)m;
      else buffery[j] = nan;
    }
    for(j=0; j<ny; j++) local[j][i] = buffery[j];
  }
  return;
}

void smooth2a_i16(int16 *in, int nx, int ny, int dx, int dy)
{
  int16 *local[ny];
  double sum, bufferx[nx], buffery[ny];
  int i,j,k,ix,iy,m;
  //int ij;

  for(i=0; i<ny; i++) local[i] = in+i*nx;

// do averaging along each row:
  for(j=0; j<ny; j++) {
    //ij = j*nx;
    for(i=0; i<nx; i++) {
//
//    if(fabs(local[j][i]-nan) < epsilon) {
//      bufferx[i] = nan;
//      continue;
//    }
      m = 0; sum = 0;
      for(k=-dx; k<=dx; k++) {
        ix = i+k;
        if(ix < 0) continue;
        if(ix > (nx-1)) continue;
        if(is_nan(local[j][ix])) continue;
        sum += local[j][ix]; m++;
      }
      if(m > 0) bufferx[i] = sum/(double)m;
      else bufferx[i] = nan;
    }
    for(i=0; i<nx; i++) local[j][i] = bufferx[i];
  }

// do averaging along each column: 
  for(i=0; i<nx; i++) {
    for(j=0; j<ny; j++) {
      //ij = j*nx+i;
//
//    if(fabs(local[j][i]-nan) < epsilon) {
//      buffery[j] = nan;
//      continue;
//    }

      m = 0; sum = 0;
      for(k=-dy; k<=dy; k++) {
        iy = j+k;
        if(iy < 0) continue;
        if(iy > (ny-1)) continue;
        if(is_nan(local[iy][i])) continue;
        sum += local[iy][i]; m++;
      }
      if(m > 0) buffery[j] = sum/(double)m;
      else buffery[j] = nan;
    }
    for(j=0; j<ny; j++) local[j][i] = buffery[j];
  }
  return;
}

