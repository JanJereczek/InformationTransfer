# InformationTransfer

## Reproducibility

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> InformationTransfer

It is authored by Jan Swierczek-Jereczek.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

## Information transfer

The relative information transfer allows to determine causal links between variables in multi-dimensional time series. It allows to make quantitative statements, as it returns a value between -100 and 100%. The methodology used here was introduced in [Liang (2021)](https://doi.org/10.3390/e23060679) and applied in [Docquier et al. (2022)](https://esd.copernicus.org/articles/14/577/2023/). The implementation of the present code was largely eased by a [repository](https://github.com/Climdyn/Liang_Index_climdyn) associated with the latter.

Following observational products are used in this repository:
 - [IMBIE data used in IPCC AR6](https://ramadda.data.bas.ac.uk/repository/entry/show?entryid=77b64c55-7166-4a06-9def-2e400398e452)
 - [ENSO-3.4 index from NOAA](https://psl.noaa.gov/data/timeseries/monthly/NINO34/)
 - [SAM index from NOAA](https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/aao.shtml)