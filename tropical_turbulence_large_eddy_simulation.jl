using Oceananigans
using Printf

include("tropical_turbulence_setup.jl")

arch = GPU()
setup = tropical_turbulence_setup(arch)

Lh = 216
Nh = 216 # horizontal resolution is uniform Lh / Nh

Nz = length(setup.z) - 1

grid = RectilinearGrid(arch,
                       size = (Nh, Nh, Nz),
                       x = (0, Nh),
                       y = (0, Nh),
                       z = setup.z,
                       topology = (Periodic, Periodic, Bounded))

model = NonhydrostaticModel(; grid,
                            advection = WENO(),
                            timestepper = :RungeKutta3,
                            tracers = (:T, :S),
                            buoyancy = setup.buoyancy,
                            forcing = setup.forcing,
                            boundary_conditions = setup.boundary_conditions,
                            closure = AnisotropicMinimumDissipation())

set!(model; setup.initial_conditions...)

#####
##### Simulation with adaptive time-stepping
##### + callback to update the forcing time index every iteration
#####

simulation = Simulation(model, Δt=1e-3, stop_time=1day)

simulation.callbacks[:update_time_index] = setup.update_time_index

wizard = TimeStepWizard(cfl=0.5, max_change=1.1, max_Δt=1minute)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

function progress(sim)
    u, v, w = sim.model.velocities

    msg = @sprintf("Iter: %d, time: %s, forcing time index: %d, max|u|: (%.2e, %.2e, %.2e)",
                   iteration(sim), prettytime(sim), setup.forcing_time_index[],
                   maximum(abs, u), maximum(abs, v), maximum(abs, w))
                   
    @info msg

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

run!(simulation)



