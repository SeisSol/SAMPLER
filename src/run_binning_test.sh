#!/bin/bash
#SBATCH -J SAMPLER
#SBATCH -o ./binning.%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=mpp3
#SBATCH --partition=mpp3_batch
#SBATCH --nodes=1-1
#SBATCH --cpus-per-task=256
#SBATCH --hint=multithread
# 56 is the maximum reasonable value for CooLMUC-2
#SBATCH --mail-type=all
#SBATCH --mail-user=m.schmeller@tum.de
#SBATCH --export=NONE
#SBATCH --time=02:00:00
module load slurm_setup
export JULIA_NUM_THREADS=8
echo "Starting SAMPLER with $JULIA_NUM_THREADS threads:"
#echo "LB=naive:"
#~/julia-1.4.1/bin/julia -O3 ./sampler.jl --load-balancer=naive -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf
#echo "LB=count:"
#~/julia-1.4.1/bin/julia -O3 ./sampler.jl --load-balancer=count -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf
echo "LB=workload:"
~/julia-1.4.1/bin/julia -O3 ./sampler.jl -m 64G --load-balancer=workload -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf
export JULIA_NUM_THREADS=28
echo "Starting SAMPLER with $JULIA_NUM_THREADS threads:"
#echo "LB=naive:"
#~/julia-1.4.1/bin/julia -O3 ./sampler.jl --load-balancer=naive -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf 
#echo "LB=count:"
#~/julia-1.4.1/bin/julia -O3 ./sampler.jl --load-balancer=count -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf 
echo "LB=workload:"
~/julia-1.4.1/bin/julia -O3 ./sampler.jl -m 64G --load-balancer=workload -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf
export JULIA_NUM_THREADS=64
echo "Starting SAMPLER with $JULIA_NUM_THREADS threads:"
#echo "LB=naive:"
#~/julia-1.4.1/bin/julia -O3 ./sampler.jl --load-balancer=naive -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf 
#echo "LB=count:"
#~/julia-1.4.1/bin/julia -O3 ./sampler.jl --load-balancer=count -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf 
echo "LB=workload:"
~/julia-1.4.1/bin/julia -O3 ./sampler.jl -m 64G --load-balancer=workload -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf 
export JULIA_NUM_THREADS=128
echo "Starting SAMPLER with $JULIA_NUM_THREADS threads:"
#echo "LB=naive:"
#~/julia-1.4.1/bin/julia -O3 ./sampler.jl --load-balancer=naive -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf 
#echo "LB=count:"
#~/julia-1.4.1/bin/julia -O3 ./sampler.jl --load-balancer=count -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf 
echo "LB=workload:"
~/julia-1.4.1/bin/julia -O3 ./sampler.jl -m 64G --load-balancer=workload -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf 
export JULIA_NUM_THREADS=256
echo "Starting SAMPLER with $JULIA_NUM_THREADS threads:"
#echo "LB=naive:"
#~/julia-1.4.1/bin/julia -O3 ./sampler.jl --load-balancer=naive -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf 
#echo "LB=count:"
#~/julia-1.4.1/bin/julia -O3 ./sampler.jl --load-balancer=count -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf 
echo "LB=workload:"
~/julia-1.4.1/bin/julia -O3 ./sampler.jl -m 64G --load-balancer=workload -o $SCRATCH/dummy.nc -t 5 $SCRATCH/output2/out-surface.xdmf $SCRATCH/output2/out.xdmf
echo "Done."
