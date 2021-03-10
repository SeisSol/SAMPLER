# SAMPLER
**S**eismic **A**ctivity **M**apping & **P**rojection for **L**ethal maritime **E**vent **R**econstruction (working title)

SAMPLER is a collection of scripts for rasterizing unstructured meshes from [SeisSol][1] into a format compatible with [sam(oa)²][2].

Additionally, it offers postprocessing options, namely depth filtering via [Kajiura's filter][3].

## Installation
First, set up [Julia][4] on your machine.
Then, in the command line, type:

    julia
    ]
    add NetCDF Interpolations FFTW ArgParse XMLDict
    ←
    exit()

Note: "←" is the backspace key.
Note: The `julia` command might not be available right after installation and you may need to manually make it available (you could add this line to your `~/.bashrc` file):

    export PATH=$PATH:$HOME/julia-1.4.1/bin

Afterwards, clone the SAMPLER repository to your machine (if you read this on the machine, you have already completed this step!):

    git clone https://gitlab.lrz.de/ge73tes/sampler

Note: In the following it is assumed that `path/to/sampler` is the installation path (the sampler folder contains this README.md file). 

## Usage
Once you have installed SAMPLER and you have the SeisSol outputs in XDMF format you can run the script:

    cd path/to/sampler/src
    export JULIA_NUM_THREADS = <number_of_hardware_threads_available>
    julia sampler.jl <options>

Note: For an explanation of the available command line options, type:

    julia sampler.jl --help

You can also have a look at the `run_sampler.sh` and `run_kajiura.sh` scripts to see how SAMPLER can be executed on the LRZ compute clusters.
More help for running on the cluster can be fount [here][5].

Once you have rasterized the SeisSol output files, you can _optionally_ use Kajiura's filter on the outputs:

    julia kajiura.jl <in_file.nc> <out_file.nc> <timestep_end>

Note: Kajiura's Filter needs to process all timesteps prior to "timestep_end".

## Examples
### Basic example: Rasterize 2D grids (seafloor and sea surface data is available)

    julia sampler.jl -m 8G --water-height=2000 -o ~/sampler-output.nc ~/seissol-outputs/out-surface.xdmf

Only rasterizes the 2D surface outputs (seafloor and sea surface elevation).

### Rasterize only seafloor (and optionally apply Kajiura's Filter)

    julia sampler.jl -m 8G --water-height=2000 -o ~/sampler-output.nc --seafloor-only ~/seissol-outputs/out-surface.xdmf
    
Only rasterizes the 2D _seafloor_ elevation over time. Optionally thereafter:
    
    julia kajiura.jl ~/sampler-output.nc ~/kajiura-output.nc 300

Applies Kajiura's Filter for 300 timesteps.
Note: You can also rasterize the SeisSol outputs fully (without `--seafloor-only`) and then apply Kajiura's Filter to them. But that usually does not make much sense.

### Rasterize seafloor using Tanioka's method

    julia sampler.jl -m 8G --water-height=2000 -o ~/sampler-output.nc --seafloor-only --tanioka ~/seissol-outputs/out-surface.xdmf

Applies Tanioka's method while rasterizing the seafloor. Needed if bathymetry is not flat.
Note: You can also use Tanioka's method when rasterizing more than just the seafloor outputs. 
It will only affect the seafloor uplift, though.

### Rasterize fully (with 3D velocity grid)

    julia sampler.jl -m 8G --water-height=2000 -o ~/sampler-output.nc -t 300 ~/seissol-outputs/out-surface.xdmf

Rasterizes all grids and variables needed for tsunami simulation, including the 3D velocity grid.
Only rasterizes timestep 300 (for example).

## Using the outputs with sam(oa)²
You need the [max-bachelor][6] branch of sam(oa)² in order to seamlessly use the output files of SAMPLER.
The other versions of sam(oa)² cannot handle multiple variables in one input file yet.
Refer to its README.md for mor information.

## Known Issues
* The `--memory-limit` or `-m` argument imposes a soft limit on the memory used. Thus, choose about half of the memory available on your machine / cluster node.


[1]: http://www.seissol.org/
[2]: https://gitlab.lrz.de/samoa/samoa
[3]: https://ci.nii.ac.jp/naid/120000866529/
[4]: https://julialang.org/downloads/
[5]: https://doku.lrz.de/display/PUBLIC/Running+serial+jobs+on+the+Linux-Cluster#RunningserialjobsontheLinuxCluster-Script-drivenSLURMjobs
[6]: https://gitlab.lrz.de/samoa/samoa/-/tree/max-bachelor