# TropicalTurbulantics.jl

Data and code to setup single column models and large eddy simulations of the equatorial ocean using [Oceananigans](https://github.com/CliMA/Oceananigans.jl),
following [Whitt et al. 2022](https://journals.ametsoc.org/view/journals/phoc/52/5/JPO-D-21-0153.1.xml).

## Downloading and processing setup data

The initial conditions, boundary conditions, and forcing data for large eddy simulations
are extracted from a simulation of the equatorial ocean using the
Regional Ocean Modeling System (ROMS). To download and process the setup data,

1. Download the [ROMS data](https://figshare.com/ndownloader/files/28415004) provided by Whitt et al. 2022.
   The forcing, boundary conditions, and initial conditions data are archived in `data_les.tar`.
   Place `data_les.tar` in the root directory `TropicalTurbulantics.jl/`.

2. Untar `data_les.tar` to the directory `TropicalTurbulantics.jl/data_les`:

```
tar xvf data_les.tar
```

3. Run the julia script `process_data.jl`, which loads the NetCDF data,
   removes spurious values, plots the cleaned data, and saves the cleaned data
   in a `JLD2` file:

```julia
julia --project process_data.jl
```

The resulting plot should look like this:

![forcing_and_bcs_and_ics_0N140W](https://github.com/glwagner/TropicalTurbulantics.jl/blob/main/forcing_and_bcs_and_ics_0N140W.png)

## Running and analyzing the single column model

```julia
julia --project tropical_turbulence_single_column_model.jl
```

The simulation takes 4-5 minutes on an M1 MacBook pro with 1 thread.
After the simulation is run, the results may be plotted with


```julia
julia --project plot_single_column_simulation.jl
```

The result should look something like


![tropical_turbulence_single_column_simulation](https://github.com/glwagner/TropicalTurbulantics.jl/blob/main/tropical_turbulence_single_column_simulation.png)

