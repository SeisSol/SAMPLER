#!/bin/bash

module load nco/4.6.2

if [ -z "$2" ]; then
    echo "Preprocessing and copying to same directory"
    ncks -O --fix_rec_dmn=all --cnk_plc=uck $1_cell.h5 $1_cell.uck.h5
    mv $1_cell.h5 $1_cell.h5.bak
    mv $1_cell.uck.h5 $1_cell.h5

    ncks -O --fix_rec_dmn=all --cnk_plc=uck $1_vertex.h5 $1_vertex.uck.h5
    mv $1_vertex.h5 $1_vertex.h5.bak
    mv $1_vertex.uck.h5 $1_vertex.h5
    echo "  $1_....h5 has been backed up and the original files are now preprocessed for use with SAMPLER"
else
    echo "Preprocessing and copying to target directory"
    ncks -O --fix_rec_dmn=all --cnk_plc=uck $1_cell.h5 $2_cell.h5
    ncks -O --fix_rec_dmn=all --cnk_plc=uck $1_vertex.h5 $2_vertex.h5
    cp -f $1.xdmf $2.xdmf
    echo "  $1_... has been copied to $2_... and preprocessed for use with SAMPLER"
fi
