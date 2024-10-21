#pragma once
/**
 * @file Contains definitions and functions used in order to process
 * L2_Cloud for ECOSTRESS.
 *
 * @author Robert Freepartner, JPL/Raytheon/JaDa Systems
 * @copyright (c) 2021 JPL, All rights reserved
 */
#include "matrix.h"
#include "tes_util.h"

/// @brief ECOSTRESS L2 Cloud Processing
/// @param l1b_rad Radiance data from the L1B_RAD file
/// @param bt11_lut_file full pathname of the LUT_Cloud_BT11.h5 file
/// @param btdiff_file full pathname of the Cloud BTdiff h5 file
/// @param rad_lut_data RAD_BT lookup file data in 6 columns
/// @param rad_lut_nlines number of lines in rad_lut_data
/// @param cloud_filename pathname to use to output the L2_CLOUD data
void process_cloud(RAD *l1b_rad,
                   const char* bt11_lut_file,
                   const char* btdiff_file,
                   double* rad_lut_data[],
                   int rad_lut_nlines,
                   const char* cloud_filename);
