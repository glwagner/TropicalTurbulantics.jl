# TropicalTurbulantics.jl

Data and code to setup a large eddy simulation of the equatorial ocean using Oceananigans,
following Whitt et al. 2022.

## Downloading and processing setup data

The initial conditions, boundary conditions, and forcing for large eddy simulations
are extracted from a simulation of the equatorial ocean using the
Regional Ocean Modeling System (ROMS).

1. Download the [ROMS data](https://figshare.com/ndownloader/files/28415004).
   The forcing, boundary conditions, and initial conditions data are archived
   in `data_les.tar`.

2. Untar `data_les.tar` to the directory `TropicalTurbulantics.jl/data_les`.

3. Run the julia script `process_data.jl`, which loads the NetCDF data,
   removes spurious values, plots the cleaned data, and saves the cleaned data
   in a `JLD2` file:

```julia
julia --project process_data.jl
```

## Running the large eddy simulation

```julia
julia --project deep_cycle_turbulence.jl
```

