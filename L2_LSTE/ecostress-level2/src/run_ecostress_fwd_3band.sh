#!/bin/sh
# Script to run the EXAMPLE_FWD code (RTTOV Interface)
# Dr. Tanvir Islam, NASA/JPL
#
########################################################################
# Set BIN directory if supplied
#BIN=`perl -e 'for(@ARGV){m/^BIN=(\S+)$/o && print "$1";}' $*`
#if [ "x$BIN" = "x" ]
#then
#  BIN=bin
#fi

######## Edit this section for your pathnames and test case input ######

# Path relative to the rttov_test directory:
#TEST=/data/local/users/tislam/tes_sandbox/data

# Paths relative to the rttov_test/${TEST} directory:
#BINDIR=/data/local/users/tislam/rttov_11.2/bin                    # BIN directory (may be set with BIN= argument)
#DATART=/data/local/users/tislam/rttov_11.2/rtcoef_rttov11/rttov9pred54L   # Coefficients directory

# Args: <rttov-exe-file> <rttov-coef-file>
EXE=$1
COEF_FILENAME=$2

# Test case input data
PROF_FILENAME="prof_in.bin"
NPROF=0        # Number of profiles defined in profile file (not used in binary case)
NLEVELS=0      # Number of profile levels (not used in binary case)
DO_SOLAR=0     # 0 = solar off / 1 = solar on
FTYPE=0        # 0 = binary / 1 = ascii (1 is not tested)

# It is possible to specify input emissivity and BRDF values below.
# Alternatively set them to zero to use RTTOV internal defaults.

########################################################################
#ARG_ARCH=`perl -e 'for(@ARGV){m/^ARCH=(\S+)$/o && print "$1";}' $*`
#if [ ! "x$ARG_ARCH" = "x" ]; then
#  ARCH=$ARG_ARCH
#fi
#if [ "x$ARCH" = "x" ];
#then
#  echo 'Please supply ARCH'
#  exit 1
#fi

#CWD=`pwd`
#cd $TEST

echo " "
echo $EXE
echo " "
echo  "Coef filename:      ${COEF_FILENAME}"
echo  "Input profile file: ${PROF_FILENAME}"
echo  "Number of profiles: ${NPROF}"
echo  "Number of levels:   ${NLEVELS}"
echo  "Do solar:           ${DO_SOLAR}"
echo  "Binary/Ascii:       ${FTYPE}"

########################################################################
# The RTTOV program requires the coefficient file to be in the default
# directory when the program is run. To make it work, a symbolic link
# name rttov_coef_file.lnk is created in the local direwctory.
rm -f rttov_coef_file.lnk
ln -s $COEF_FILENAME rttov_coef_file.lnk

########################################################################

$EXE << EOF
rttov_coef_file.lnk, Coefficient filename
${PROF_FILENAME}, Input profile filename
${NPROF}        , Number of profiles
${NLEVELS}      , Number of levels
${DO_SOLAR}     , Turn solar radiation on/off
3               , Number of channels
2 4 5           , Channel Number(s)
1               , Number of threads
${FTYPE}        , Binary/Ascii
EOF

if [ $? -ne 0 ]; then
  echo " "
  echo "RTTOV FAILED"
  echo " "
  exit 1
fi
rm -f rttov_coef_file.lnk
exit
# Previous channel list:
# 1, 1, 0.00001, 0.0 , channel 2 - 8.3μm TIR,  valid, emissivity=0.0, BRDF=0.0
# 2, 1, 0.00001, 0.0 , channel 3 - 8.6μm TIR,  valid, emissivity=0.0, BRDF=0.0
# 3, 1, 0.00001, 0.0 , channel 4 - 9.1μm TIR,  valid, emissivity=0.0, BRDF=0.0
# 4, 1, 0.00001, 0.0 , channel 5 - 11.2μm TIR, valid, emissivity=0.0, BRDF=0.0
# 5, 1, 0.00001, 0.0 , channel 6 - 12.0μm TIR, valid, emissivity=0.0, BRDF=0.0
