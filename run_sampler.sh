#!/bin/bash
#SBATCH -J SAMPLER
#SBATCH -o ./full-run.%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --nodes=1-1
#SBATCH --cpus-per-task=28
#SBATCH --export=NONE
#SBATCH --time=02:00:00
module load slurm_setup
export JULIA_NUM_THREADS=28
echo "Starting SAMPLER with $JULIA_NUM_THREADS threads:"
date
julia ./sampler.jl --water-height=2000 -m 25G -o $SCRATCH/dummy.nc $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf
date
echo "Done."
