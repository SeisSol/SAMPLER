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

SRCDIR="/path/to/seissol/output/folder"
DSTDIR="/path/to/sampler/output/folder"
OUTNAME='name-of-xdmf-file' # without ".xdmf" attached
mkdir -p $DSTDIR
echo "=== Working on $SRCDIR"

echo "=== Starting SAMPLER with $JULIA_NUM_THREADS threads:"
fvars="U=>fU,V=>fV,W=>fW,u=>fu,v=>fv,w=>fw"
svars="U=>sU,V=>sV,W=>sW,u=>su,v=>sv,w=>sw"

julia --project=. ./src/SAMPLER.jl --water-height=2000 --seafloor-vars $fvars --surface-vars $svars -m 25G -o "$DSTDIR/rasterized.nc" "$SRCDIR/${OUTNAME}.xdmf"
date
echo "=== Complete."
