#!/bin/bash

# Splits the input NetCDF files into two files each, one containing bathymetry and one containing displacements.
#
# Usage:
# ./split_bath_disp.sh <file1> [file2 [...]]
#
# A file "file.nc" will be split into "file-disp.nc" and "file-bathy.nc" in the same directory.
# The file itself remains untouched.

module load nco

FILES=$@
for fn in $FILES
do
    # Cut ".nc" extension off filena,e
    FILE=$(echo $fn | rev | cut -f 2- -d '.' | rev)

    echo "Working on $FILE.nc, outputting to $FILE-disp.nc..."
    ncks -v d $FILE.nc $FILE-disp.nc
    ncrename -v d,z $FILE-disp.nc
    echo "DISP done"

    echo "Working on $FILE.nc, outputting to $FILE-bathy.nc..."
    ncks -v b $FILE.nc $FILE-bathy.nc
    ncrename -v b,z $FILE-bathy.nc
    echo "BATHY done"
done
