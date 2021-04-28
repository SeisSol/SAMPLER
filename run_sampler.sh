#!/bin/bash
#SBATCH -J SAMPLER
#SBATCH -o ./%A_%a.out
#SBATCH -e ./%A_%a.err
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --partition=micro
#SBATCH --account=pr63qo
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.schmeller@tum.de
#SBATCH --time=00:30:00
module load slurm_setup
. ~/.bashrc

# exit when any command fails
set -e
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

export JULIA_NUM_THREADS=48
date

SRCDIR="/hppfs/work/pn68fi/ga24dib3/shared/ascete_waterlayer_wo_gravity/run_0/output"
DSTDIR="$WORK_pn68fi/lukas"
OUTNAME='out-surface'
mkdir -p $DSTDIR
echo "=== Working on $SRCDIR"

#cp $SRCDIR/${OUTNAME}.xdmf $DSTDIR
#echo '=== Copied XDMF to own directory'

#echo '=== Pre-processing HDF5 files'
#$HOME/ascete-tools/hdf5_preprocess.sh $SRCDIR/$OUTNAME $DSTDIR/$OUTNAME
#echo '=== Pre-processed HDF5 files'

echo "=== Starting SAMPLER with $JULIA_NUM_THREADS threads:"
fvars="U=>fU,V=>fV,W=>fW,u=>fu,v=>fv,w=>fw"
svars="U=>sU,V=>sV,W=>sW,u=>su,v=>sv,w=>sw"

julia ./sampler.jl --water-height=2000 --seafloor-vars $fvars --surface-vars $svars -m 25G -o "$DSTDIR/rasterized.nc" "$SRCDIR/${OUTNAME}.xdmf"
date
echo "=== Complete."
