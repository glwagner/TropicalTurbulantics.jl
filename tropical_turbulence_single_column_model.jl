using Oceananigans
using Oceananigans.Units
using Printf

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
     CATKEVerticalDiffusivity,
     TurbulentKineticEnergyEquation,
     MixingLength

include("tropical_turbulence_setup.jl")

arch = CPU()
hi_res_setup = tropical_turbulence_setup(arch; Nz=216)
z = hi_res_setup.z[1:8:end]
Nz = length(z) - 1
setup = tropical_turbulence_setup(arch; Nz, z)

grid = RectilinearGrid(arch,
                       size = Nz,
                       z = setup.z,
                       topology = (Flat, Flat, Bounded))

default_catke = CATKEVerticalDiffusivity()

model = HydrostaticFreeSurfaceModel(; grid,
                                    tracers = (:T, :S, :e),
                                    buoyancy = setup.buoyancy,
                                    forcing = setup.forcing,
                                    boundary_conditions = setup.boundary_conditions,
                                    closure = default_catke)

set!(model; e=1e-9, setup.initial_conditions...)

#####
##### Simulation with adaptive time-stepping
##### + callback to update the forcing time index every iteration
#####

stop_time = 34days
simulation = Simulation(model; Δt=10minute, stop_time)

simulation.callbacks[:update_time_index] = setup.update_time_index

function progress(sim)
    u, v, w = sim.model.velocities
    e = sim.model.tracers.e

    msg = @sprintf("Iter: %d, time: %s, forcing time index: %d, max|u|: (%.2e, %.2e), max e: %.2e",
                   iteration(sim), prettytime(sim), setup.forcing_time_index[],
                   maximum(abs, u), maximum(abs, v), maximum(abs, e))
                   
    @info msg

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

u, v, w = model.velocities
b = Oceananigans.BuoyancyModels.buoyancy(model)

Ri = ∂z(b) / (∂z(u)^2 + ∂z(v)^2)
outputs = merge(model.velocities, model.tracers, (; Ri, b)) 
filename = string("tropical_turbulence_single_column_model_Nz", Nz, ".jld2")

simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs,
                                                    schedule = TimeInterval(20minutes);
                                                    filename,
                                                    overwrite_existing = true)

run!(simulation)

