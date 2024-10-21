// smooth2d
//
/// @file This is a more accurate replication of the MatLab smooth2a function
/// that implements a performance improvement over 
/// O(nrows * ncols * 2*(Nr+1) * 2*(Nc+1)) to O(nrows * ncols * 2).
///
/// From MatLab: smooth2a
/// This function smooths the data in matrixIn using a mean filter over a
/// rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
/// element "i" by the mean of the rectange centered on "i".  Any NaN
/// elements are ignored in the averaging.  If element "i" is a NaN, then it
/// will be preserved as NaN in the output.  At the edges of the matrix,
/// where you cannot build a full rectangle, as much of the rectangle that
/// fits on your matrix is used (similar to the default on Matlab's builtin
/// function "smooth").
///
/// @author Robert Freepartner, JPL/Raytheon/JaDa
/// @date March 2016
/// @copyright (c) Copyright 2016, Jet Propulsion Laboratories, Pasadena, CA

#include <assert.h>
#include <stdlib.h>
#include "lste_lib.h"

typedef struct
{
    double  total;  ///< total of all non-NaN values in partial column
    int     count;  ///< number of non-NaN values included in total
} SmoothingCol;

void constructSmoothingCol(SmoothingCol *sc)
{
    sc->total = 0;
    sc->count = 0;
}

typedef struct
{
    int      ndatacols; ///< number of columns per row of data
    SmoothingCol *cols; ///< Array of ndatacols SmoothingCol structs
} SmoothingSquare;

void constructSmoothingSquare(SmoothingSquare *ss, int cols)
{
    ss->ndatacols = cols;
    ss->cols = (SmoothingCol*)malloc(cols * sizeof(SmoothingCol));
    assert(ss->cols != NULL);
    int c;
    for (c = 0; c < cols; c++)
    {
        constructSmoothingCol(&ss->cols[c]);
    }
}

void destructSmoothingSquare(SmoothingSquare *ss)
{
    free(ss->cols);
}

// "matrixIn": original matrix sized nrows x ncols
// matrixIn is replaced by smoothed matrix for output.
// "Nr": number of points used to smooth rows
// "Nc": number of points to smooth columns.
void smooth2d(double *matrixIn, int nrows, int ncols, int Nr, int Nc)
{
    int r, c;
    double d;
    double sqtotal;
    int sqcount;
    size_t matsize = nrows * ncols * sizeof(double);
    double *matrixOut = (double*)malloc(matsize);
    assert(matrixOut != NULL);

    // Initialize the SmoothingSquare.
    SmoothingSquare sq;
    constructSmoothingSquare(&sq, ncols);

    // Initialize each SmoothingCol to represent row = -1.
    for (c = 0; c < ncols; c++)
    {
        for (r = 0; r < Nr && r < nrows; r++)
        {
            d = matrixIn[r*ncols + c];
            if (!is_nan(d))
            {
                sq.cols[c].total += d;
                sq.cols[c].count++;
            }
        }
    }

    // Loop over rows of data.
    for (r = 0; r < nrows; r++)
    {
        // Add the next row of data to the SmoothingCols.
        int addrow = r + Nr;
        // If the next row is beyond the last, no update needed.
        if (addrow < nrows)
        {
            // Loop over all matrixIn columns.
            for (c = 0; c < ncols; c++)
            {
                d = matrixIn[addrow*ncols + c];
                // If NaN, do not add.
                if (!is_nan(d))
                {
                    sq.cols[c].total += d;
                    sq.cols[c].count++;
                }
            }
        }
        
        // Initialize the square total and count to represent col = -1;
        sqtotal = 0.0;
        sqcount = 0;
        for (c = 0; c < Nc; c++)
        {
            sqtotal += sq.cols[c].total;
            sqcount += sq.cols[c].count;
        }
        
        // Loop over the columns.
        for (c = 0; c < ncols; c++)
        {
            // Add the next column to the square total and count.
            int addcol = c + Nc;
            // No need to add if beyond the last column.
            if (addcol < ncols)
            {
                sqtotal += sq.cols[addcol].total;
                sqcount += sq.cols[addcol].count;
            }

            // Get the data value of the current row, col.
            d = matrixIn[r*ncols + c];

            // If current data value is NaN, leave it NaN.
            if (is_nan(d))
            {
                matrixOut[r*ncols + c] = d;
            }
            
            // Else, if square count is zero, leave the value.
            // Note: this is to prevent divide by zero, which 
            // would result in inf. The case should only happen 
            // when the whole square is NaN, in which case d
            // would be NaN anyway. Just being safe here.
            else if (sqcount == 0)
            {
                matrixOut[r*ncols + c] = d;
            }
            
            // Else, replace value with mean (sqtotal/sqcount).
            else
            {
                matrixOut[r*ncols + c] = sqtotal / (double)sqcount;
            }
            
            // Subtract the leftmost column value from sqtotal and the count 
            // from sqcount, but sqcount must not be < zero. If sqcount is
            // zero, sq.cols[leftc] must contain zeroes.
            if (sqcount > 0)
            {
                int leftc = c - Nc;
                // No need to adjust if leftc is < zero.
                if (leftc >= 0)
                {
                    sqtotal -= sq.cols[leftc].total;
                    sqcount -= sq.cols[leftc].count;
                }
            }
        } // end loop over columns

        // Prepare for the next row by subtracting the top row 
        // from the SmoothingCols. No need if last row.
        if (r < nrows - 1)
        {
            int toprow = r - Nr;
            // No need if toprow is not in dataset.
            if (toprow >= 0)
            {
                for (c = 0; c < ncols; c++)
                {
                    d = matrixIn[toprow * ncols + c];
                    if (!is_nan(d))
                    {
                        sq.cols[c].total -= d;
                        if (sq.cols[c].count > 0) // just being safe here
                        {
                            sq.cols[c].count--;
                        }
                    }
                }
            }
        }
    } // end loop over rows

    // Copy matrixOut to matrixIn.
    memcpy(matrixIn, matrixOut, matsize);
    // Delete memory for matrixOut.
    free(matrixOut);
    // destruct square.
    destructSmoothingSquare(&sq);
}
