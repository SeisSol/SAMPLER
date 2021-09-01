#!/bin/bash
#SBATCH -J SAMPLER
#SBATCH -o ./%A_%a.out
#SBATCH -e ./%A_%a.err
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --partition=micro
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --time=00:30:00
# TODO: Specify cluster / account / mail settings / etc.
module load slurm_setup
. ~/.bashrc

export JULIA_NUM_THREADS=48
date

SRCDIR="/path/to/seissol/output/folder"
DSTDIR="/path/to/sampler/output/folder"
OUTNAME='name-of-xdmf-file' # without ".xdmf" attached
mkdir -p $DSTDIR
echo "=== Working on $SRCDIR"

echo "=== Starting SAMPLER with $JULIA_NUM_THREADS threads:"
# Just exmaples, replace or leave out.
fvars="U=>fU,V=>fV,W=>fW,u=>fu,v=>fv,w=>fw"
svars="U=>sU,V=>sV,W=>sW,u=>su,v=>sv,w=>sw"

julia --project=. ./src/main.jl --water-height=2000 --seafloor-vars $fvars --surface-vars $svars -m 25G -o "$DSTDIR/rasterized.nc" "$SRCDIR/${OUTNAME}.xdmf"
date
echo "=== Complete."
