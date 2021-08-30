# SAMPLER
**S**eismic **A**ctivity **M**apping & **P**rojection for **L**ethal maritime **E**vent **R**econstruction (working title)

SAMPLER is a program for rasterizing unstructured meshes output by [SeisSol][1] into a format compatible with [sam(oa)²][2].

Additionally, it offers postprocessing options, namely depth filtering via [Kajiura's filter][3] and converting horizontal displacements to vertical ones using [Tanioka's method][7].

## Installation
First, set up Julia on your machine.
The recommended ways are:

* On SuperMUC-NG: [Julia on SuperMUC-NG and HPC Systems, Getting Started][8]
* On Linux and Windows, locally: [Julia Downloads][4]

Afterwards, clone the SAMPLER repository to your machine (if you read this on the machine, you have already completed this step!):

```bash
git clone https://github.com/SeisSol/SAMPLER
```

_Note: In the following it is assumed that `path/to/sampler` is the installation path (the sampler folder contains this README.md file)._

Then, in `path/to/sampler`, run:

```bash
path/to/sampler> julia
julia> ]
(@v1.x) pkg> activate .
(SAMPLER) pkg> instantiate
(SAMPLER) pkg> test
```

When all tests pass, SAMPLER is installed correctly.

## Usage
Once you have installed SAMPLER and you have the SeisSol outputs in XDMF format you can run the script:

```bash
cd path/to/sampler/
export JULIA_NUM_THREADS=<num_threads>
julia --project=. src/SAMPLER.jl <options>
```

_Note: In Windows (PowerShell), replace `export ...` with:_

```powershell
$env:JULIA_NUM_THREADS = <num_threads>
```
_Note: For an explanation of the available command line options, type:_

```bash
julia --project=. src/SAMPLER.jl --help
```

You can also have a look at the `run_sampler.sh` and `run_kajiura.sh` scripts to see how SAMPLER can be executed on the LRZ compute clusters.
More help for running on the cluster can be found [here][5].

Once you have rasterized the SeisSol output files, you can optionally use Kajiura's filter on the outputs:

```bash
julia scripts/kajiura.jl <in_file.nc> <out_file.nc> <timestep_end>
```

_Note: Kajiura's Filter needs to process all timesteps prior to `timestep_end`._

## Examples
### Basic example: Rasterize 2D grids (seafloor and sea surface data is available)

```bash
julia --project=. src/SAMPLER.jl -m 8G --water-height=2000 -o ~/sampler-output.nc ~/seissol-outputs/out-surface.xdmf
```

Only rasterizes the 2D surface outputs (seafloor and sea surface elevation).

### Rasterize only seafloor (and optionally apply Kajiura's Filter)

```bash
julia --project=. src/SAMPLER.jl -m 8G --water-height=2000 -o ~/sampler-output.nc --seafloor-only ~/seissol-outputs/out-surface.xdmf
```

Only rasterizes the 2D _seafloor_ elevation over time. Optionally thereafter:

```bash    
julia --project=. src/SAMPLER.jl ~/sampler-output.nc ~/kajiura-output.nc 300
```

Applies Kajiura's Filter for 300 timesteps.

_Note: You can also rasterize the SeisSol outputs fully (without `--seafloor-only`) and then apply Kajiura's Filter to them. But that usually does not make much sense._

### Rasterize seafloor using Tanioka's method

```bash
julia --project=. src/SAMPLER.jl -m 8G --water-height=2000 -o ~/sampler-output.nc --seafloor-only --tanioka ~/seissol-outputs/out-surface.xdmf
```

Applies Tanioka's method while rasterizing the seafloor. Needed if bathymetry is not flat.

_Note: You can also use Tanioka's method when rasterizing more than just the seafloor outputs. It will only affect the seafloor uplift._

### Rasterize fully (with 3D velocity grid)

```bash
julia --project=. src/SAMPLER.jl -m 8G --water-height=2000 -o ~/sampler-output.nc -s 300 ~/seissol-outputs/out-surface.xdmf
```

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

```bash
--seafloor-vars "W,b=>bathy"
```

Output `W` with the same name, remap mesh geometry height to `bathy`.

## Known Issues
* The `--memory-limit` or `-m` argument imposes a soft limit on the memory used. Thus, choose about half of the memory available on your machine / cluster node.


[1]: http://www.seissol.org/
[2]: https://gitlab.lrz.de/samoa/samoa
[3]: https://ci.nii.ac.jp/naid/120000866529/
[4]: https://julialang.org/downloads/
[5]: https://doku.lrz.de/display/PUBLIC/Running+serial+jobs+on+the+Linux-Cluster#RunningserialjobsontheLinuxCluster-Script-drivenSLURMjobs
[6]: https://gitlab.lrz.de/samoa/samoa/-/tree/max-bachelor
[7]: https://dx.doi.org/10.1029/96GL00736
[8]: https://doku.lrz.de/display/PUBLIC/FAQ%3A+Julia+on+SuperMUC-NG+and+HPC+Systems#FAQ:JuliaonSuperMUCNGandHPCSystems-Gettingstarted