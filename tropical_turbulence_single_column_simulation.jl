using Oceananigans
using Oceananigans.Units
using Printf
using JLD2

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
     CATKEVerticalDiffusivity,
     TurbulentKineticEnergyEquation,
     MixingLength

include("tropical_turbulence_setup.jl")

arch = CPU()
hi_res_setup = tropical_turbulence_setup(arch; Nz=216)
z = hi_res_setup.z[1:4:end]
Nz = length(z) - 1
setup = tropical_turbulence_setup(arch; Nz, z)

grid = RectilinearGrid(arch,
                       size = Nz,
                       z = setup.z,
                       topology = (Flat, Flat, Bounded))

@load "optimal_catke.jld2" optimal_catke
@show optimal_catke

optimal_mixing_length = optimal_catke.mixing_length
mixing_length_parameters = Dict(name => getproperty(optimal_mixing_length, name)
                                for name in propertynames(optimal_mixing_length))

optimal_turbulent_kinetic_energy_equation = optimal_catke.turbulent_kinetic_energy_equation
turbulent_kinetic_energy_equation_parameters =
    Dict(name => getproperty(optimal_turbulent_kinetic_energy_equation, name)
         for name in propertynames(optimal_turbulent_kinetic_energy_equation))

mixing_length_parameters[:Cᵇ] = 1.0
turbulent_kinetic_energy_equation_parameters[:Cᵂϵ] = 0.0

mixing_length = MixingLength(; mixing_length_parameters...)
turbulent_kinetic_energy_equation =
    TurbulentKineticEnergyEquation(; turbulent_kinetic_energy_equation_parameters...)

closure = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)

# Default:
# closure = CATKEVerticalDiffusivity()

model = HydrostaticFreeSurfaceModel(; grid, closure,
                                    tracers = (:T, :S, :e),
                                    buoyancy = setup.buoyancy,
                                    forcing = setup.forcing,
                                    boundary_conditions = setup.boundary_conditions)

set!(model; e=1e-9, setup.initial_conditions...)

#####
##### Simulation with adaptive time-stepping
##### + callback to update the forcing time index every iteration
#####

Δt = 30.0 #5minute
stop_time = 34days
simulation = Simulation(model; Δt, stop_time)

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

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

u, v, w = model.velocities
T = model.tracers.T
b = Oceananigans.BuoyancyModels.buoyancy(model)
N² = @at (Center, Center, Center) ∂z(b)
S² = @at (Center, Center, Center) ∂z(u)^2 + ∂z(v)^2
κc = model.diffusivity_fields.κc
κu = model.diffusivity_fields.κu

# Note the location for consistency w/ LES
wT = @at (Center, Center, Center) κc * ∂z(T) * -1
wu = @at (Center, Center, Center) κu * ∂z(u) * -1
Ri = @at (Center, Center, Center) ∂z(b) / (∂z(u)^2 + ∂z(v)^2)

# Let's save these too...
wT_ccf = @at (Center, Center, Face) κc * ∂z(T) * -1 
wu_ccf = @at (Center, Center, Face) κu * ∂z(u) * -1
Ri_ccf = @at (Center, Center, Face) ∂z(b) / (∂z(u)^2 + ∂z(v)^2)

outputs = merge(model.velocities, model.tracers, (; Ri, b, N², S², κc, κu, wT, wu, wT_ccf, wu_ccf, Ri_ccf)) 
filename = string("single_column_tropical_turbulence_tiny_time_step_Nz", Nz, ".jld2")

simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs,
                                                    schedule = TimeInterval(20minutes);
                                                    filename,
                                                    overwrite_existing = true)

run!(simulation)

