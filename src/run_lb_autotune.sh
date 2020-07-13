#!/bin/bash
#SBATCH -J SAMPLER
#SBATCH -o ./lb-autotune.%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=mpp3
#SBATCH --partition=mpp3_batch
#SBATCH --nodes=1-1
#SBATCH --cpus-per-task=64
#SBATCH --constraint=cache
# 56 is the maximum reasonable value for CooLMUC-2
#SBATCH --mail-type=all
#SBATCH --mail-user=m.schmeller@tum.de
#SBATCH --export=NONE
#SBATCH --time=00:45:00
module load slurm_setup
export JULIA_NUM_THREADS=64
echo "Starting SAMPLER with $JULIA_NUM_THREADS threads:"
~/julia-1.4.1/bin/julia -O3 ./sampler.jl --lb-autotune -m 56G -o $SCRATCH/dummy.nc $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf
echo "Done."
