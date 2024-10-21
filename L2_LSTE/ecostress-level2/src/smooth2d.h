#pragma once
/**
 * @file This header file contains definitions and prototypes for 
 * the smooth2d function.
 *
 * @author Robert Freepartner, JPL/Raytheon/JaDa
 * @date March 2016
 * @copyright (c) Copyright 2016, Jet Propulsion Laboratories, Pasadena, CA
 */
 
/**
 * @brief This is a more accurate replication of the MatLab smooth2a function
 * that implements mean filter smoothing over a rectangular area.
 *
 * It contains a performance improvement over the simpler method of
 * O(nrows * ncols * 2*(Nr+1) * 2*(Nc+1)) to O(nrows * ncols * 2).
 *
 * @param[in,out]   matrixIn    data matrix sized nrows x ncols
 * @param           nrows       number of rows of data
 * @param           ncols       number of cols per each row of data
 * @param           Nr          number of rows up and down for smoothing
 */
void smooth2d(double *matrixIn, int nrows, int ncols, int Nr, int Nc);
