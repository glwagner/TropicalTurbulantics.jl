using Oceananigans
using Printf

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
     CATKEVerticalDiffusivity,
     TurbulentKineticEnergyEquation,
     MixingLength

include("tropical_turbulence_setup.jl")

arch = CPU()
Nz = 27
setup = tropical_turbulence_setup(arch; Nz)

#Nz = length(setup.z) - 1
grid = RectilinearGrid(arch,
                       size = Nz,
                       z = setup.z,
                       topology = (Flat, Flat, Bounded))


turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(CᵂwΔ = 2.858323560177654,
                                                                   Cᵂu★ = 6.429661459350047,
                                                                   C⁺D  = 0.6332855351727852,
                                                                   C⁻D  = 0.776966635028129,
                                                                   CᶜD  = 3.031920275020819,
                                                                   CᵉD  = 0.04035710595894642)

mixing_length = MixingLength(C⁺c  = 0.13681143163095005,
                             C⁺u  = 0.2635185903882231,
                             C⁺e  = 9.82562434662029,
                             Cᵇ   = 0.8322522991885453,
                             C⁻c  = 0.4175147980679791,
                             C⁻u  = 0.49250780910121605,
                             C⁻e  = 7.410282104867351,
                             CRiᶜ = 0.4933951366483328,
                             CRiʷ = 0.400084096820722,
                             Cᶜc  = 8.23976241523068,
                             Cᶜe  = 3.695385457245675,
                             Cᵉc  = 0.4685456558627356,
                             Cᵉe  = 1.8793640062659647,
                             Cˢᶜ  = 0.19659617642661853)

optimal_catke = CATKEVerticalDiffusivity(; turbulent_kinetic_energy_equation, mixing_length)

model = HydrostaticFreeSurfaceModel(; grid,
                                    tracers = (:T, :S, :e),
                                    buoyancy = setup.buoyancy,
                                    forcing = setup.forcing,
                                    boundary_conditions = setup.boundary_conditions,
                                    closure = optimal_catke)

set!(model; e=1e-9, setup.initial_conditions...)

#####
##### Simulation with adaptive time-stepping
##### + callback to update the forcing time index every iteration
#####

simulation = Simulation(model, Δt=1minute, stop_time=33days)

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

simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs,
                                                    schedule = TimeInterval(20minutes),
                                                    filename = "tropical_turbulence_single_column_model.jld2",
                                                    overwrite_existing = true)

run!(simulation)

