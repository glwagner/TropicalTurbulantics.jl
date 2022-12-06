using Oceananigans
using Oceananigans.Units
using Oceananigans.Architectures: architecture, arch_array
using SeawaterPolynomials
using JLD2
using Printf

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

"""
    tropical_turbulence_setup(arch; Nz=216, datapath = "forcing_and_ics_0N140W.jld2")

Return a setup for simulating tropical deep cycle turbulence.
"""
function tropical_turbulence_setup(arch = CPU();
                                   z = nothing,
                                   Nz = 216,
                                   datapath = "forcing_and_bcs_and_ics_0N140W.jld2")
                                     

    # Load data
    file = jldopen(datapath)

    # Load the "original" z-grid for Whitt et al. 2022.
    user_z = z # save the user spec

    z = file["z"]

    # Add one more cell interface at z = -108 meters.
    # Not sure if this is correct. However, we must have
    # that length(z) = Nz + 1, where Nz is the number of cells.
    # Otherwise the grid is underdetermined.
    @info string("First cell interface is z[1] = ", z[1])
    @info string("Adding one more interface at -108 meters")

    z = vcat(-108, z)

    @info string("First cell interface is now z[1] = ", z[1])
    @assert length(z) == 217

    # Forcing time
    forcing_times = file["time_seconds"]

    # Large-scale state
    # Can use this to implement restoring --- not implemented yet.
    U = file["U"]
    V = file["V"]
    T = file["T"]
    S = file["S"]

    # Thermodynamic properties
    # Linear equation of state:
    #
    #   (ρ - ρᵣ) / ρᵣ ∼ α T - β S
    #
    # where α=thermal_expansion and β=haline_contraction
    thermal_expansion = file["thermal_expansion"]
    haline_contraction = file["haline_contraction"]

    # Large-scale tendency forcing + solar insolation
    Fᵁ = file["Fᵁ"]
    Fⱽ = file["Fⱽ"]

    # Note: combine ROMS tendency + insolution
    Fᵀ = file["Fᵀ"] .+ file["Fᴵ"]

    # Surface fluxes
    Qᵁ_surface = file["Qᵁ_surface"]
    Qⱽ_surface = file["Qⱽ_surface"]
    Qᵀ_surface = file["Qᵀ_surface"]
    Qˢ_surface = file["Qˢ_surface"]

    # Bottom fluxes
    Qᵁ_bottom = file["Qᵁ_bottom"]
    Qⱽ_bottom = file["Qⱽ_bottom"]
    Qᵀ_bottom = file["Qᵀ_bottom"]
    Qˢ_bottom = file["Qˢ_bottom"]

    close(file)

    equation_of_state = LinearEquationOfState(; thermal_expansion, haline_contraction)
    buoyancy = SeawaterBuoyancy(; equation_of_state)

    #####
    ##### Set up forcing
    #####

    # Callback that updates the forcing time index.
    # This callback runs on the CPU.
    forcing_time_index = n = Ref(1)

    function update_time_index_func(sim)
        n = findfirst(t -> t > time(simulation), forcing_times)

        if isnothing(n)
            msg = @sprintf("Simulation time %.2e was not found in forcing times.", time(simulation))
            error(msg)
        end

        forcing_time_index[] = n
        return nothing
    end

    # regrid if necessary
    if Nz != 216 || !isnothing(user_z)
        if isnothing(user_z)
            user_z = (-108, 0)
        else
            Nz = length(user_z) - 1
        end

        original_vertical_grid = RectilinearGrid(size=216; z, topology=(Flat, Flat, Bounded))
        new_vertical_grid      = RectilinearGrid(size=Nz; z=user_z, topology=(Flat, Flat, Bounded))

        orig_field = CenterField(original_vertical_grid)
        new_field = CenterField(new_vertical_grid)

        Nt = length(forcing_times)
        new_Fᵁ = zeros(Nz, Nt)
        new_Fⱽ = zeros(Nz, Nt)
        new_Fᵀ = zeros(Nz, Nt)
        new_U = zeros(Nz, Nt)
        new_V = zeros(Nz, Nt)
        new_T = zeros(Nz, Nt)
        new_S = zeros(Nz, Nt)

        for (new, orig) in [(new_Fᵁ, Fᵁ),
                            (new_Fⱽ, Fⱽ),
                            (new_Fᵀ, Fᵀ),
                            (new_U,  U),
                            (new_V,  V),
                            (new_T,  T),
                            (new_S,  S)]

            for n = 1:Nt
                orig_field .= reshape(orig[:, n], 1, 1, 216)
                regrid!(new_field, orig_field)
                new[:, n] .= interior(new_field, :)
            end
        end

        # Just pretend nothing happened
        Fᵁ = new_Fᵁ
        Fⱽ = new_Fⱽ 
        Fᵀ = new_Fᵀ
        U  = new_U  
        V  = new_V  
        T  = new_T  
        S  = new_S  
        z  = znodes(Face, new_vertical_grid)
    end

    # Convert CPU arrays to GPU arrays if necessary:
    tᶠ = arch_array(arch, forcing_times)

    Qᵁ_surface = arch_array(arch, Qᵁ_surface)
    Qⱽ_surface = arch_array(arch, Qⱽ_surface)
    Qᵀ_surface = arch_array(arch, Qᵀ_surface)
    Qˢ_surface = arch_array(arch, Qˢ_surface)
    Qᵁ_bottom  = arch_array(arch, Qᵁ_bottom)
    Qⱽ_bottom  = arch_array(arch, Qⱽ_bottom)
    Qᵀ_bottom  = arch_array(arch, Qᵀ_bottom)
    Qˢ_bottom  = arch_array(arch, Qˢ_bottom)

    Fᵁ = arch_array(arch, Fᵁ)
    Fⱽ = arch_array(arch, Fⱽ)

    # Note: combine ROMS tendency + insolution
    Fᵀ = arch_array(arch, Fᵀ)
    
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

    u_bcs = FieldBoundaryConditions(top=u_surface_bc) #, bottom=u_bottom_bc)
    v_bcs = FieldBoundaryConditions(top=v_surface_bc) #, bottom=v_bottom_bc)
    T_bcs = FieldBoundaryConditions(top=T_surface_bc) #, bottom=T_bottom_bc)
    S_bcs = FieldBoundaryConditions(top=S_surface_bc) #, bottom=S_bottom_bc)

    #####
    ##### LES setup with AMD closure and RK3 time-stepping
    #####

    forcing = (u=u_forcing, v=v_forcing, T=T_forcing)
    boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs)
    update_time_index = Callback(update_time_index_func)

    Uᵢ = reshape(U[:, 1], 1, 1, Nz)
    Vᵢ = reshape(V[:, 1], 1, 1, Nz)
    Tᵢ = reshape(T[:, 1], 1, 1, Nz)
    Sᵢ = reshape(S[:, 1], 1, 1, Nz)

    initial_conditions = (u=Uᵢ, v=Vᵢ, T=Tᵢ, S=Sᵢ)

    return (; update_time_index,
              forcing_time_index,
              z,
              buoyancy,
              forcing,
              boundary_conditions,
              initial_conditions)
end
