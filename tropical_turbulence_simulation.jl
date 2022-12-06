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

Qᵁ_surface = file["Qᵁ_surface"]
Qⱽ_surface = file["Qⱽ_surface"]
Qᵀ_surface = file["Qᵀ_surface"]
Qˢ_surface = file["Qˢ_surface"]

Qᵁ_bottom = file["Qᵁ_bottom"]
Qⱽ_bottom = file["Qⱽ_bottom"]
Qᵀ_bottom = file["Qᵀ_bottom"]
Qˢ_bottom = file["Qˢ_bottom"]

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
function update_time_index(sim)
    nn = findfirst(t -> t > time(simulation), forcing_times)
    n[] = nn
    return nothing
end

@inline function interp_forcing(i, j, k, grid, clock, model_fields, params)
    n_reference, tᶠ, F = params
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

    F = F₁ + dFdt * (t - t₁)

    #return F
    return ifelse(F > 1.0, 0.0, F)
end

@inline function interp_bc(i, j, grid, clock, model_fields, params)
    n_reference, tᶠ, Q = params
    n = n_reference[]

    @inbounds begin
        Q₁ = Q[n-1]
        Q₂ = Q[n]
        t₁ = tᶠ[n-1]
        t₂ = tᶠ[n]
    end

    # Linear interpolation
    t = clock.time
    dQdt = (Q₂ - Q₁) / (t₂ - t₁)

    Q = Q₁ + dQdt * (t - t₁)

    return Q
end

# Convert CPU arrays to GPU arrays if necessary:
tᶠ = arch_array(architecture(grid), forcing_times)

Fᵁ = arch_array(architecture(grid), Fᵁ)
Fⱽ = arch_array(architecture(grid), Fⱽ)

# Note: combine ROMS tendency + insolution
Fᵀ = arch_array(architecture(grid), Fᵀ .+ Fᴵ)

Qᵁ_surface = arch_array(architecture(grid), Qᵁ_surface)
Qⱽ_surface = arch_array(architecture(grid), Qⱽ_surface)
Qᵀ_surface = arch_array(architecture(grid), Qᵀ_surface)
Qˢ_surface = arch_array(architecture(grid), Qˢ_surface)
Qᵁ_bottom  = arch_array(architecture(grid), Qᵁ_bottom)
Qⱽ_bottom  = arch_array(architecture(grid), Qⱽ_bottom)
Qᵀ_bottom  = arch_array(architecture(grid), Qᵀ_bottom)
Qˢ_bottom  = arch_array(architecture(grid), Qˢ_bottom)

u_forcing = Forcing(interp_forcing, discrete_form=true, parameters=(n, tᶠ, Fᵁ))
v_forcing = Forcing(interp_forcing, discrete_form=true, parameters=(n, tᶠ, Fⱽ))
T_forcing = Forcing(interp_forcing, discrete_form=true, parameters=(n, tᶠ, Fᵀ))

u_surface_bc = FluxBoundaryCondition(interp_bc, discrete_form=true, parameters=(n, tᶠ, Qᵁ_surface))
v_surface_bc = FluxBoundaryCondition(interp_bc, discrete_form=true, parameters=(n, tᶠ, Qⱽ_surface))
T_surface_bc = FluxBoundaryCondition(interp_bc, discrete_form=true, parameters=(n, tᶠ, Qᵀ_surface))
S_surface_bc = FluxBoundaryCondition(interp_bc, discrete_form=true, parameters=(n, tᶠ, Qˢ_surface))
u_bottom_bc  = FluxBoundaryCondition(interp_bc, discrete_form=true, parameters=(n, tᶠ, Qᵁ_bottom))
v_bottom_bc  = FluxBoundaryCondition(interp_bc, discrete_form=true, parameters=(n, tᶠ, Qⱽ_bottom))
T_bottom_bc  = FluxBoundaryCondition(interp_bc, discrete_form=true, parameters=(n, tᶠ, Qᵀ_bottom))
S_bottom_bc  = FluxBoundaryCondition(interp_bc, discrete_form=true, parameters=(n, tᶠ, Qˢ_bottom))

u_bcs = FieldBoundaryConditions(top=u_surface_bc, bottom=u_bottom_bc)
v_bcs = FieldBoundaryConditions(top=v_surface_bc, bottom=v_bottom_bc)
T_bcs = FieldBoundaryConditions(top=T_surface_bc, bottom=T_bottom_bc)
S_bcs = FieldBoundaryConditions(top=S_surface_bc, bottom=S_bottom_bc)

#####
##### LES setup with AMD closure and RK3 time-stepping
#####

model = NonhydrostaticModel(; grid,
                            advection = WENO(),
                            timestepper = :RungeKutta3,
                            tracers = (:T, :S),
                            buoyancy = SeawaterBuoyancy(; equation_of_state),
                            forcing = (u=u_forcing, v=v_forcing, T=T_forcing),
                            boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs),
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

simulation.callbacks[:update_time_index] = Callback(update_time_index)

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

