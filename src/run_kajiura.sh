#!/bin/bash
#SBATCH -J KAJIURA
#SBATCH -o ./kaji-run.%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --nodes=1-1
#SBATCH --cpus-per-task=28
#SBATCH --export=NONE
#SBATCH --time=02:00:00
module load slurm_setup
export JULIA_NUM_THREADS=56
echo "Starting Kajiura with $JULIA_NUM_THREADS threads:"
~/julia-1.4.1/bin/julia -O3 ./kajiura.jl $SCRATCH/full_out.nc $SCRATCH/kaj_out_300.nc 300
