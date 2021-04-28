#!/bin/bash
#SBATCH -J SAMPLER
#SBATCH -o ./%A_%a.out
#SBATCH -e ./%A_%a.err
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --partition=micro
#SBATCH --account=pr63qo
#SBATCH --nodes=1-12
#SBATCH --export=NONE
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.schmeller@tum.de
#SBATCH --time=00:30:00
#SBATCH --array=1-12
module load slurm_setup
. ~/.bashrc
export JULIA_NUM_THREADS=48
date

case $SLURM_ARRAY_TASK_ID in
    1)  MAJ=25; MIN=25;;
    2)  MAJ=25; MIN=50;;
    3)  MAJ=25; MIN=75;;
    4)  MAJ=30; MIN=25;;
    5)  MAJ=30; MIN=50;;
    6)  MAJ=30; MIN=75;;
    7)  MAJ=40; MIN=25;;
    8)  MAJ=40; MIN=50;;
    9)  MAJ=40; MIN=75;;
    10) MAJ=45; MIN=25;;
    11) MAJ=45; MIN=50;;
    12) MAJ=45; MIN=75;;
esac

RELDIR="${MAJ}km/$MIN"
mkdir -p $WORK_pr63qo/$RELDIR
echo "=== Working on $RELDIR"

cp /hppfs/work/pr63qo/ru94haf3/share/Max/$RELDIR/output-surface.xdmf $WORK_pr63qo/$RELDIR/
echo '=== Copied XDMF to own directory'

echo '=== Pre-processing HDF5 files'
$HOME/ascete-tools/hdf5_preprocess.sh /hppfs/work/pr63qo/ru94haf3/share/Max/$RELDIR/output-surface $WORK_pr63qo/$RELDIR/output-surface
echo '=== Pre-processed HDF5 files'

echo "=== Starting SAMPLER with $JULIA_NUM_THREADS threads:"
julia ./sampler.jl --water-height=2000 --kajiura -m 40G -r 500,500,500 -o $WORK_pr63qo/$RELDIR/displacement.nc $WORK_pr63qo/$RELDIR/output-surface.xdmf
date
echo "=== Complete."
