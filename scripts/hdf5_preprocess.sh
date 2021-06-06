#!/bin/bash

# Converts chunked HDF5 files and such with unlimited dimensions to unchunked, limited-dimension HDF5 files.
# Needed for some SeisSol outputs because the HDF5 library in Julia can't handle chunked inputs.
#
# Usage:
# ./hdf5_preprocess.sh <in_xdmf_stem> [out_dir]
#
# The in_xdmf_stem parameter is the common part of the hdf5 files' filenames. (e.g. "out" for "out_cell.h5" and "out_vertex.h5")
# When specifying an out_dir, the corresponding XDMF file will also be copied to there.

module load nco/4.6.2

if [ -z "$2" ]; then
    echo "Preprocessing and copying to same directory"
    ncks -O --fix_rec_dmn=all --cnk_plc=uck $1_cell.h5 $1_cell.uck.h5
    ncks -O --fix_rec_dmn=all --cnk_plc=uck $1_vertex.h5 $1_vertex.uck.h5
#    mv $1_cell.h5 $1_cell.h5.bak
#    mv $1_vertex.h5 $1_vertex.h5.bak
#    mv $1_cell.uck.h5 $1_cell.h5
#    mv $1_vertex.uck.h5 $1_vertex.h5
else
    echo "Preprocessing and copying to target directory"
    ncks -O --fix_rec_dmn=all --cnk_plc=uck $1_cell.h5 $2_cell.h5
    ncks -O --fix_rec_dmn=all --cnk_plc=uck $1_vertex.h5 $2_vertex.h5
    cp -f $1.xdmf $2.xdmf
fi
