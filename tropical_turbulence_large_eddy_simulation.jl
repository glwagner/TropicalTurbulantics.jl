using Oceananigans
using Printf

using Oceanostics.TKEBudgetTerms: TurbulentKineticEnergy

include("tropical_turbulence_setup.jl")

arch = GPU()
Nz = 216
Lh = 306
Nh = 360

setup = tropical_turbulence_setup(arch; Nz)

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

# Seed with a random velocity field
wᵢ(x, y, z) = 1e-4 * (1 - 2rand())
set!(model, w=wᵢ)

#####
##### Simulation with adaptive time-stepping
##### + callback to update the forcing time index every iteration
#####

simulation = Simulation(model, Δt=1.0, stop_time=1day)

simulation.callbacks[:update_time_index] = setup.update_time_index

wizard = TimeStepWizard(cfl=0.7, max_change=1.1, max_Δt=1minute)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

wall_time = Ref(time_ns())

function progress(sim)
    u, v, w = sim.model.velocities
    elapsed = 1e-9 * (time_ns() - wall_time[])

    msg = @sprintf("Iter: %d, time: %s, Δt: %s, wall time: %s, forcing time index: %d, max|u|: (%.2e, %.2e, %.2e)",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt),
                   prettytime(elapsed), setup.forcing_time_index[],
                   maximum(abs, u), maximum(abs, v), maximum(abs, w))
                   
    @info msg

    wall_time[] = time_ns()

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

prefix = string("tropical_turbulence_Nz", Nz)

u, v, w = model.velocities
T, S = model.tracers
b = Oceananigans.BuoyancyModels.buoyancy(model)

S²_op = @at (Center, Center, Face) ∂z(u)^2 + ∂z(v)^2
S² = Field(S²_op)
Ri = ∂z(b) / S²

U = Field(Average(u, dims=(1, 2)))
V = Field(Average(v, dims=(1, 2)))
T_avg = Average(T, dims=(1, 2))
S_avg = Average(S, dims=(1, 2))
Ri_avg = Average(Ri, dims=(1, 2))

e = TurbulentKineticEnergy(model, U=U, V=V)
E = Average(e, dims=(1, 2))

averaged_output = (; u=U, v=V, T=T_avg, S=S_avg, e=E, Ri=Ri_avg)
fields_output = (; u, v, w, T, S, e)

simulation.output_writers[:avg] = JLD2OutputWriter(model, averaged_output;
                                                   schedule = TimeInterval(20minutes),
                                                   filename = string(prefix, "_averages.jld2"),
                                                   overwrite_existing = false)

simulation.output_writers[:xy] = JLD2OutputWriter(model, fields_output;
                                                  indices = (:, :, Nz),
                                                  schedule = TimeInterval(20minutes),
                                                  filename = string(prefix, "_xy.jld2"),
                                                  overwrite_existing = false)

simulation.output_writers[:yz] = JLD2OutputWriter(model, fields_output;
                                                  indices = (1, :, :),
                                                  schedule = TimeInterval(20minutes),
                                                  filename = string(prefix, "_yz.jld2"),
                                                  overwrite_existing = false)

simulation.output_writers[:xz] = JLD2OutputWriter(model, fields_output;
                                                  indices = (:, 1, :),
                                                  schedule = TimeInterval(20minutes),
                                                  filename = string(prefix, "_xz.jld2"),
                                                  overwrite_existing = false)

run!(simulation)

