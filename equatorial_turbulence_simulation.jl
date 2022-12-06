# Equatorial turublence simulation forced by 1/20ᵒ ROMS output following
#
# Whitt, D. B. et al, "Simulation and Scaling of the Turbulence Vertical Heat Transport
# and Deep-Cycle Turbulence across the Equatorial Pacific Cold Tongue", JPO, 2022
#
# See README.md for instructions for downloading and pre-processing the data.

using Oceananigans
using Oceananigans.Units
using Oceananigans.Architectures: architecture, arch_array
using SeawaterPolynomials
using JLD2
using Printf

#####
##### Load the data that was processed by `process_data.jl`
#####

setup_filepath = "forcing_and_ics_0N140W.jld2"

file = jldopen(setup_filepath)

z = file["z"]
forcing_times = file["time_seconds"]
U = file["U"]
V = file["V"]
T = file["T"]
S = file["S"]
Fᵁ = file["Fᵁ"]
Fⱽ = file["Fⱽ"]
Fᵀ = file["Fᵀ"]
Fᴵ = file["Fᴵ"]
thermal_expansion = file["thermal_expansion"]
haline_contraction = file["haline_contraction"]

close(file)

# Append the last vertical level interface (?)
z = vcat(-108, z)
Nz = length(z) - 1 # Whitt et al. 2022 have Nz = 216

# 1m horizontal spacing
Lh = 216
Nh = 1 #216

arch = CPU() # GPU()

grid = RectilinearGrid(arch,
                       size = (Nh, Nh, Nz),
                       x = (0, Lh),
                       y = (0, Lh),
                       z = z,
                       topology = (Periodic, Periodic, Bounded))

equation_of_state = LinearEquationOfState(; thermal_expansion, haline_contraction)

#####
##### Set up forcing
#####

n = Ref(1)

# Callback to add to Simulation
function update_forcing_time_index(sim)
    nn = findfirst(t -> t > time(simulation), forcing_times)
    n[] = nn
    return nothing
end

@inline function tendency_forcing(i, j, k, grid, clock, model_fields, params)
    n_reference, F, tᶠ = params
    n = n_reference[]

    @inbounds begin
        F₁ = F[k, n-1]
        F₂ = F[k, n]
        t₁ = tᶠ[n-1]
        t₂ = tᶠ[n]
    end

    # Linear interpolation
    t = clock.time
    dFdt = (F₂ - F₁) / (t₂ - t₁)

    return 0 #F₁ + dFdt * (t - t₁)
end

Fᵁ = arch_array(architecture(grid), Fᵁ)
Fⱽ = arch_array(architecture(grid), Fⱽ)
Fᵀ = arch_array(architecture(grid), Fᵀ .+ Fᴵ)
tᶠ = arch_array(architecture(grid), forcing_times)

u_forcing = Forcing(tendency_forcing, discrete_form=true, parameters=(n, Fᵁ, tᶠ))
v_forcing = Forcing(tendency_forcing, discrete_form=true, parameters=(n, Fⱽ, tᶠ))
T_forcing = Forcing(tendency_forcing, discrete_form=true, parameters=(n, Fᵀ, tᶠ))

#####
##### LES setup with AMD closure and RK3 time-stepping
#####

model = NonhydrostaticModel(; grid,
                            advection = WENO(),
                            timestepper = :RungeKutta3,
                            tracers = (:T, :S),
                            buoyancy = SeawaterBuoyancy(; equation_of_state),
                            forcing = (u=u_forcing, v=v_forcing, T=T_forcing),
                            closure = AnisotropicMinimumDissipation())

Uᵢ = reshape(U[:, 1], 1, 1, Nz)
Vᵢ = reshape(V[:, 1], 1, 1, Nz)
Tᵢ = reshape(T[:, 1], 1, 1, Nz)
Sᵢ = reshape(S[:, 1], 1, 1, Nz)

set!(model; u=Uᵢ, v=Vᵢ, T=Tᵢ, S=Sᵢ)

#####
##### Simulation with adaptive time-stepping
##### + callback to update the forcing time index every iteration
#####

simulation = Simulation(model, Δt=1e-3, stop_iteration=10)

simulation.callbacks[:update_forcing_time] = Callback(update_forcing_time_index)

wizard = TimeStepWizard(cfl=0.5, max_change=1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

function progress(sim)
    u, v, w = sim.model.velocities

    msg = @sprintf("Iter: %d, time: %s, forcing time index: %d, max|u|: (%.2e, %.2e, %.2e)",
                   iteration(sim), prettytime(sim), n[],
                   maximum(abs, u), maximum(abs, v), maximum(abs, w))
                   
    @info msg

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

run!(simulation)

