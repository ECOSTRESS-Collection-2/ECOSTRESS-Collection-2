#ifndef BOOL_H
#define BOOL_H

/*
!C*****************************************************************************

File: bool.h

!Description: 
    Bool (logical) type declaration.

!Input Parameters: (not applicable)

!Input/output Parameters: (not applicable)

!Output Parameters: (not applicable)

!Revision History:
  Revision 01.00  1998/11/06
  Robert E. Wolfe
  1. Original version.
  
!Team-unique Header:
  This software was developed by:

    MODIS Land Science Team Support Group for the National Aeronautics and 
    Space Administration, Goddard Space Flight Center, under 
    contract NAS5-32373.

!References and Credits:
  
    MODIS Science Team Members:
      Eric Vermote
      MODIS Land Science Team           University of Maryland
      eric@kratmos.gsfc.nasa.gov        NASA's GSFC Code 923
      phone: 301-614-5512               Greenbelt, MD 20771

    Developers:
      Robert E. Wolfe
      MODIS Land Team Support Group     Raytheon TSC
      robert.e.wolfe.1@gsfc.nasa.gov    1616 McCormick Dr.
      phone: 301-614-5508               Upper Marlboro, MD 20772
      NASA/GSFC Code 922, Greenbelt Rd., Greenbelt, MD 20771

!Design Notes: (none)

!END
****************************************************************************/

/* Define the bool (logical) type */

#ifdef true
#undef true
#endif
#ifdef false
#undef false
#endif
typedef enum {false = 0, true = 1} bool;

#endif
