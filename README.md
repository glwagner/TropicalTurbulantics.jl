# TropicalTurbulantics.jl

Data and code to setup a large eddy simulation of the equatorial ocean using [Oceananigans](https://github.com/CliMA/Oceananigans.jl),
following [Whitt et al. 2022](https://journals.ametsoc.org/view/journals/phoc/52/5/JPO-D-21-0153.1.xml).

## Downloading and processing setup data

The initial conditions, boundary conditions, and forcing for large eddy simulations
are extracted from a simulation of the equatorial ocean using the
Regional Ocean Modeling System (ROMS).

1. Download the [ROMS data](https://figshare.com/ndownloader/files/28415004) provided by Whitt et al. 2022.
   The forcing, boundary conditions, and initial conditions data are archived in `data_les.tar`.

2. Untar `data_les.tar` to the directory `TropicalTurbulantics.jl/data_les`.

3. Run the julia script `process_data.jl`, which loads the NetCDF data,
   removes spurious values, plots the cleaned data, and saves the cleaned data
   in a `JLD2` file:

```julia
julia --project process_data.jl
```

The resulting plot should look like this:

![forcing_and_bcs_and_ics_0N140W](https://user-images.githubusercontent.com/15271942/205716011-31131754-71a2-4cb5-bb49-fa4fa8047e67.png)

## Running the large eddy simulation

```julia
julia --project deep_cycle_turbulence.jl
```

