# SAMPLER Supplementary Scripts

## `split_bath_disp.sh`

Splits one or more NetCDF files into two files each, one containing only bathymetry and one only containing displacements in z-direction.

This is needed since SAMPLER outputs _one_ NetCDF file containing all rasterized variables while sam(oa)² requires each variable to be in a separate file.

Usage:
```bash
./split_bath_disp.sh <file1.nc> [file2.nc [...]]
```

The input NetCDF files have to have a variable called `b` (bathymetry) and one called `d` (displacement).
The default settings of SAMPLER also produce those output variables.

The output files will be named `file1-disp.nc`, `file1-bathy.nc` and so on, respectively.
The NetCDF variable in each output file will be renamed to `z` as required by sam(oa)², the input files remain unchanged.

The `split_bath_disp.sh` script requires [NetCDF Operators (NCO)][1] to be installed.


[1]: http://nco.sourceforge.net/