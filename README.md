# SAMPLER
**S**eismic **A**ctivity **M**apping & **P**rojection for **L**ethal maritime **E**vent **R**econstruction (working title)

SAMPLER is a program for rasterizing unstructured meshes output by [SeisSol][1] into a format compatible with [sam(oa)²][2].

Additionally, it offers postprocessing options, namely depth filtering via [Kajiura's filter][3] and converting horizontal displacements to vertical ones using [Tanioka's method][7].

## Installation
First, set up [Julia][4] on your machine. Please stick to version **1.4.x** for now as some libraries have changed behavior in newer versions.
Then, in the command line, type:

    julia
    ]
    add NetCDF Interpolations FFTW ArgParse XMLDict HDF5 Infiltrator
    ←
    exit()

_Note: `←` is the backspace key._

_Note: The `julia` command might not be available right after installation and you may need to manually make it available (you could add this line to your `~/.bashrc` file):_

    export PATH=$PATH:$HOME/julia-1.4.1/bin

_Note: In Windows, your path should already be set correctly after the installation of Julia._

Afterwards, clone the SAMPLER repository to your machine (if you read this on the machine, you have already completed this step!):

    git clone https://github.com/SeisSol/SAMPLER

_Note: In the following it is assumed that `path/to/sampler` is the installation path (the sampler folder contains this README.md file)._

## Usage
Once you have installed SAMPLER and you have the SeisSol outputs in XDMF format you can run the script:

    cd path/to/sampler/
    export JULIA_NUM_THREADS=<num_threads>
    julia sampler.jl <options>

_Note: In Windows (PowerShell), replace `export ...` with:_

    $env:JULIA_NUM_THREADS = <num_threads>

_Note: For an explanation of the available command line options, type:_

    julia sampler.jl --help

You can also have a look at the `run_sampler.sh` and `run_kajiura.sh` scripts to see how SAMPLER can be executed on the LRZ compute clusters.
More help for running on the cluster can be found [here][5].

Once you have rasterized the SeisSol output files, you can optionally use Kajiura's filter on the outputs:

    julia kajiura.jl <in_file.nc> <out_file.nc> <timestep_end>

_Note: Kajiura's Filter needs to process all timesteps prior to `timestep_end`._

## Examples
### Basic example: Rasterize 2D grids (seafloor and sea surface data is available)

    julia sampler.jl -m 8G --water-height=2000 -o ~/sampler-output.nc ~/seissol-outputs/out-surface.xdmf

Only rasterizes the 2D surface outputs (seafloor and sea surface elevation).

### Rasterize only seafloor (and optionally apply Kajiura's Filter)

    julia sampler.jl -m 8G --water-height=2000 -o ~/sampler-output.nc --seafloor-only ~/seissol-outputs/out-surface.xdmf
    
Only rasterizes the 2D _seafloor_ elevation over time. Optionally thereafter:
    
    julia kajiura.jl ~/sampler-output.nc ~/kajiura-output.nc 300

Applies Kajiura's Filter for 300 timesteps.

_Note: You can also rasterize the SeisSol outputs fully (without `--seafloor-only`) and then apply Kajiura's Filter to them. But that usually does not make much sense._

### Rasterize seafloor using Tanioka's method

    julia sampler.jl -m 8G --water-height=2000 -o ~/sampler-output.nc --seafloor-only --tanioka ~/seissol-outputs/out-surface.xdmf

Applies Tanioka's method while rasterizing the seafloor. Needed if bathymetry is not flat.

_Note: You can also use Tanioka's method when rasterizing more than just the seafloor outputs. It will only affect the seafloor uplift._

### Rasterize fully (with 3D velocity grid)

    julia sampler.jl -m 8G --water-height=2000 -o ~/sampler-output.nc -s 300 ~/seissol-outputs/out-surface.xdmf

Rasterizes all grids and variables needed for tsunami simulation, including the 3D velocity grid.
Only rasterizes timestep 300 (for example).

## Using the outputs with sam(oa)²
You need the [max-bachelor][6] branch of sam(oa)² in order to seamlessly use the output files of SAMPLER.
The other versions of sam(oa)² cannot handle multiple variables in one input file yet.
Refer to its README.md for mor information.

## Output format
SAMPLER will always produce a regular grid in NetCDF format, with point data according to the command line parameters given.
You can control which input variables get mapped to which outputs via the `--seafloor-vars`, `--surface-vars`, and `--volumetric-vars` arguments.
The standard mapping is:

| Mesh region | Mappings    |
|-------------|-------------|
| Seafloor    | `W => d`    |
|             | `b => b`*   |
| Surface     | `W => eta`  |
| Volumetric  | `u => u`    |
|             | `v => v`    |

\* `b` is not a variable in SeisSol outputs but is understood by SAMPLER as being the geometry of the 2D mesh. Thus, you can use the `b => ...` mapping to output mesh geometry.

Additionally, `x`, `y` and `time` will always be output as dimensions and variables and cannot be renamed.

Example mapping:

    --seafloor-vars "W,b=>bathy"

Output `W` with the same name, remap mesh geometry height to `bathy`.

## Known Issues
* The `--memory-limit` or `-m` argument imposes a soft limit on the memory used. Thus, choose about half of the memory available on your machine / cluster node.
* Currently, chunking in HDF5 inputs is not supported due to library constraints. If you get an error when using such inputs, use `scripts/hdf5_preprocess.sh` on them first.


[1]: http://www.seissol.org/
[2]: https://gitlab.lrz.de/samoa/samoa
[3]: https://ci.nii.ac.jp/naid/120000866529/
[4]: https://julialang.org/downloads/
[5]: https://doku.lrz.de/display/PUBLIC/Running+serial+jobs+on+the+Linux-Cluster#RunningserialjobsontheLinuxCluster-Script-drivenSLURMjobs
[6]: https://gitlab.lrz.de/samoa/samoa/-/tree/max-bachelor
[7]: https://dx.doi.org/10.1029/96GL00736