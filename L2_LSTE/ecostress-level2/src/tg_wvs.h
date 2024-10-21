#pragma once
/** 
 * @file This header file declares the datatypes and functions 
 * for the tg_wvs calculation. The calculation uses a coefficient
 * data file selected by sensor type (Aqua, Terra) and other inputs
 * to calculate Tg = Surface Brightness temperature for 3 bands.
 *
 * @author Dr. Glynn Hulley, October 2015, JPL 
 * @author Robert Freepartner, JPL, January 2016
 *
 * Algoritm design by Glynn Hulley, JPL
 *
 * @copyright Copyright (c) JPL, All rights reserved
 */
 
 #include "matrix.h"
 
/**
 * @brief Computes the Tg_WVS Surface Brightness temperature for 3 bands.
 *
 * @param[out]  tg      3D matrix of result Tg values (3 bands x 1 km)
 * @param[in]   pwv     2D matrix of 1 km total water vapor (cm)
 * @param[in]   senszen 2D matrix of 1 km sensor zenith
 * @param[in]   emis_aster 2D matrix of adjusted ASTER GED data
 * @param[in]   lrad    3D matrix of 1 km Sensor Radiance (3 bands x 1 km)
 * @param[in]   solzen  2D matrix of solar zenith (1 km or 5 km)
 * @param[in]   coef_file 	Pathname of file containing coefficients.
 * @param[in]   lut  	    Brightness Temperature Lookup Table 
 * @param[in]   nlut_lines  Number of lines in the LUT file
 * @return      0=success -1=error
 */
 
int tg_wvs(Mat3d *tg, Matrix *pwv, Matrix *senszen,
           Matrix *emis_aster, Mat3d *lrad, Matrix *solzen, 
           const char *coef_file, double *lut[], int nlut_lines);
