# Data processing for an equatorial turbulence simulation forced by 1/20ᵒ ROMS output, following
#
# Whitt, D. B. et al, "Simulation and Scaling of the Turbulence Vertical Heat Transport
# and Deep-Cycle Turbulence across the Equatorial Pacific Cold Tongue", JPO, 2022
#
# See README.md for more.
#
# Notes:
#
# * The variable "z" extends from 0 to 107.5 m, so we add the "bottom" interface at 108m.
# * The bottommost cells in the ROMS data and forcing data are problematic, so we overwrite.

using NCDatasets
using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using JLD2
using Dates
using Printf

function average_out!(Q, badvalue=0)
    for n = 1:length(Q)
        if Q[n] == badvalue
            if n == length(Q)
                Q[n] = Q[n-1]
            elseif n == 1
                Q[n] = Q[n+1]
            else
                Q[n] = 1/2 * (Q[n-1] + Q[n+1])
            end
        end
    end
    return nothing
end

make_plot = false
dir = "data_les"
filename = "ROMS_PSH_6HRLIN_0N140W_360x360x216_22OCT2020.nc"
jld2filename = "forcing_and_bcs_and_ics_0N140W.jld2"
figname = "forcing_and_bcs_and_ics_0N140W.png"
filepath = joinpath(dir, filename)
dataset = Dataset(filepath)

z = dataset["z"][:]
z[end] = 0
time_dates = dataset["time"][:]
Nt = length(time_dates)

time_seconds = datetime2unix.(time_dates)
time_seconds .-= time_seconds[1]
time_hours = time_seconds ./ hours

U  = dataset["ume"][:, :]
V  = dataset["vme"][:, :]
T  = dataset["tempme"][:, :]
S  = dataset["saltme"][:, :]
Ri = dataset["RIG"][:, :]
N² = dataset["N2"][:, :]
S² = dataset["S2"][:, :]
wu = dataset["uw"][:, :]
wv = dataset["vw"][:, :]
wT = dataset["tempw"][:, :]
wS = dataset["saltw"][:, :]

# Fix bottommost row
for ϕ in (U, V, T, S, Ri, N², S², wu, wv, wT, wS)
    ϕ[216, :] .= ϕ[215, :]
end

@printf("extrema(U) = (%.2e, %.2e) \n", extrema(U)...) 
@printf("extrema(V) = (%.2e, %.2e) \n", extrema(V)...)
@printf("extrema(T) = (%.2e, %.2e) \n", extrema(T)...)
@printf("extrema(S) = (%.2e, %.2e) \n", extrema(S)...)

@show dataset

Fᵁ = dataset["dUdtFORCE"][:, :]
Fⱽ = dataset["dVdtFORCE"][:, :]
Fᵀ = dataset["dTdtFORCE"][:, :]
Fᴵ = dataset["dTdtSOLAR"][:, :]

Fᵁ[216, :] .= Fᵁ[215, :]  
Fⱽ[216, :] .= Fⱽ[215, :] 
Fᵀ[216, :] .= Fᵀ[215, :] 
Fᴵ[216, :] .= Fᴵ[215, :] 

@printf("extrema(Fᵁ) = (%.2e, %.2e) \n", extrema(Fᵁ)...) 
@printf("extrema(Fⱽ) = (%.2e, %.2e) \n", extrema(Fⱽ)...)
@printf("extrema(Fᵀ) = (%.2e, %.2e) \n", extrema(Fᵀ)...)
@printf("extrema(Fᴵ) = (%.2e, %.2e) \n", extrema(Fᴵ)...)

# Note: Dan's LES data doesn't include dSdtFORCE
# Fˢ = dataset["dSdtFORCE"][:, :]
# Fˢ[216, :] .= Fˢ[215, :] 
# @printf("extrema(Fˢ) = (%.2e, %.2e) \n", extrema(Fˢ)...)

ρᵣ = dataset["rho0"]

Jᵁ_bottom = - dataset["nududzbot"][:]
Jⱽ_bottom = - dataset["nudvdzbot"][:]
Jᵀ_bottom = - dataset["kappadtdzbot"][:]
Jˢ_bottom = - dataset["kappadsdzbot"][:]

Jᵁ_surface = - dataset["nududztop"][:]
Jⱽ_surface = - dataset["nudvdztop"][:]
Jᵀ_surface = - dataset["kappadtdztop"][:]
Jˢ_surface = - dataset["kappadsdztop"][:]

# Fix first value
Jᵁ_surface[1] = Jᵁ_surface[2]
Jⱽ_surface[1] = Jⱽ_surface[2]
Jᵀ_surface[1] = Jᵀ_surface[2]
Jˢ_surface[1] = Jˢ_surface[2]

for Q in [Jᵁ_bottom, 
          Jⱽ_bottom, 
          Jᵀ_bottom, 
          Jˢ_bottom, 
          Jᵁ_surface, 
          Jⱽ_surface,
          Jᵀ_surface,
          Jˢ_surface]

    average_out!(Q, 0)
end

@printf("extrema(Jᵁ_surface) = (%.2e, %.2e) \n", extrema(Jᵁ_surface)...) 
@printf("extrema(Jⱽ_surface) = (%.2e, %.2e) \n", extrema(Jⱽ_surface)...)
@printf("extrema(Jᵀ_surface) = (%.2e, %.2e) \n", extrema(Jᵀ_surface)...)
@printf("extrema(Jˢ_surface) = (%.2e, %.2e) \n", extrema(Jˢ_surface)...)

@printf("extrema(Jᵁ_bottom) = (%.2e, %.2e) \n", extrema(Jᵁ_bottom)...) 
@printf("extrema(Jⱽ_bottom) = (%.2e, %.2e) \n", extrema(Jⱽ_bottom)...)
@printf("extrema(Jᵀ_bottom) = (%.2e, %.2e) \n", extrema(Jᵀ_bottom)...)
@printf("extrema(Jˢ_bottom) = (%.2e, %.2e) \n", extrema(Jˢ_bottom)...)

thermal_expansion = first(dataset["alpha"][:])

# Note sign convention
haline_contraction = - first(dataset["beta"][:])

close(dataset)

@printf("thermal expansion = %.2e \n", thermal_expansion)
@printf("haline contraction = %.2e \n", haline_contraction)

if make_plot
    using GLMakie
    
    fig = Figure(size=(1600, 1600))
    
    axU = Axis(fig[1, 1], xlabel="Time (hr)", ylabel="z (m)")
    axV = Axis(fig[1, 2], xlabel="Time (hr)", ylabel="z (m)")
    axT = Axis(fig[3, 1], xlabel="Time (hr)", ylabel="z (m)")
    axS = Axis(fig[3, 2], xlabel="Time (hr)", ylabel="z (m)")
    
    hmU = heatmap!(axU, time_hours, z, permutedims(U), colormap=:balance, colorrange=(-1, 1))
    hmV = heatmap!(axV, time_hours, z, permutedims(V), colormap=:balance, colorrange=(-1, 1))
    hmT = heatmap!(axT, time_hours, z, permutedims(T), colormap=:thermal, colorrange=(20, 26.1))
    hmS = heatmap!(axS, time_hours, z, permutedims(S), colormap=:haline, colorrange=(34.9, 35.3))
    
    Colorbar(fig[1, 0], hmU, vertical=true, flipaxis=false, label="Zonal velocity (m s⁻¹)")
    Colorbar(fig[1, 3], hmV, vertical=true, label="Meridional velocity (m s⁻¹)")
    Colorbar(fig[3, 0], hmT, vertical=true, flipaxis=false, label="Temperature (ᵒC)")
    Colorbar(fig[3, 3], hmS, vertical=true, label="Salinity (g kg⁻¹)")
    
    axJU = Axis(fig[2, 1], xlabel="Time (hrs)", ylabel="Zonal momentum flux (m² s⁻²)")
    axJV = Axis(fig[2, 2], xlabel="Time (hrs)", ylabel="Meridional momentum flux (m² s⁻²)")
    axJT = Axis(fig[4, 1], xlabel="Time (hrs)", ylabel="Temperature flux (m² s⁻²)")
    axJS = Axis(fig[4, 2], xlabel="Time (hrs)", ylabel="Salt flux (m² s⁻²)")
    
    lines!(axJU, time_hours, Jᵁ_surface, label="surface")
    lines!(axJV, time_hours, Jⱽ_surface, label="surface")
    lines!(axJT, time_hours, Jᵀ_surface, label="surface")
    lines!(axJS, time_hours, Jˢ_surface, label="surface")
    
    lines!(axJU, time_hours, Jᵁ_bottom, label="bottom")
    lines!(axJV, time_hours, Jⱽ_bottom, label="bottom")
    lines!(axJT, time_hours, Jᵀ_bottom, label="bottom")
    lines!(axJS, time_hours, Jˢ_bottom, label="bottom")
    
    axislegend(axJU)
    
    axI = Axis(fig[5, 1], xlabel="Time (hrs)", ylabel="z (m)")
    hmI = heatmap!(axI, time_hours, z, permutedims(Fᴵ), colorrange=(0, 1e-5))
    Colorbar(fig[5, 0], hmI, vertical=true, flipaxis=false,
             label="Solar insolation temperature flux divergence (ᵒC m s⁻¹)")
    
    display(fig)
    
    save(figname, fig)
end

#####
##### Save all the data
#####

rm(jld2filename, force=true)
file = jldopen(jld2filename, "a+")

file["z"] = z
file["time_dates"] = time_dates
file["time_seconds"] = time_seconds
file["ρᵣ"] = ρᵣ
file["thermal_expansion"]  = thermal_expansion
file["haline_contraction"] = haline_contraction

file["U"] = U
file["V"] = V
file["T"] = T
file["S"] = S

file["Jᵁ_surface"] = Jᵁ_surface
file["Jⱽ_surface"] = Jⱽ_surface
file["Jᵀ_surface"] = Jᵀ_surface
file["Jˢ_surface"] = Jˢ_surface

file["Jᵁ_bottom"] = Jᵁ_bottom
file["Jⱽ_bottom"] = Jⱽ_bottom
file["Jᵀ_bottom"] = Jᵀ_bottom
file["Jˢ_bottom"] = Jˢ_bottom

file["Fᵁ"] = Fᵁ 
file["Fⱽ"] = Fⱽ 
file["Fᵀ"] = Fᵀ 
file["Fᴵ"] = Fᴵ

close(file)

Nz = length(z)

grid = RectilinearGrid(size = (1, 1, Nz),
                       halo = (3, 3, 3),
                       x = (0, 306),
                       y = (0, 306),
                       z = vcat(-108, z),
                       topology = (Periodic, Periodic, Bounded))


ψ = Field{Nothing, Nothing, Center}(grid)

function set_from_data!(ψt, data)
    @info "Setting $(summary(ψt))..."
    start_time = time_ns()
    Nz, Nt = size(data)
    for n = 1:Nt
        set!(ψ, reshape(data[:, n], 1, 1, Nz))
        fill_halo_regions!(ψ)
        set!(ψt, ψ, n)
    end
    elapsed = 1e-9 * (time_ns() - start_time)
    @info string("    ... done (", prettytime(elapsed), ")")
    return nothing
end

function set_from_data!(many_ψt::Tuple, many_data::Tuple)
    @info "Setting fields..."
    for ψt in many_ψt
        @info summary(ψt)
    end

    start_time = time_ns()
    Nz, Nt = size(many_data[1])
    for n = 1:Nt
        for (ψ, data) in zip(many_ψt, many_data)
            set!(ψ, reshape(data[:, n], 1, 1, Nz))
            fill_halo_regions!(ψ)
            set!(ψt, ψ, n)
        end
    end

    elapsed = 1e-9 * (time_ns() - start_time)
    @info string("    ... done (", prettytime(elapsed), ")")

    return nothing
end



backend = OnDisk()
path = "tropical_turbulence_whitt2022_0N140W.jld2"
rm(path; force=true)

# Fields
Ut  = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="u")
Vt  = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="v")
Tt  = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="T")
St  = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="S")
Rit = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="Ri")
N²t = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="N²")
S²t = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="S²")
wut = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="wu")
wvt = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="wv")
wTt = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="wT")
wSt = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="wS")

#=
# Surface bcs
τxtt = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, time_seconds; backend, path, name="τxt")
τytt = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, time_seconds; backend, path, name="τyt")
JTtt = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, time_seconds; backend, path, name="JTt")
JStt = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, time_seconds; backend, path, name="JSt")

# Bottom bcs
τxbt = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, time_seconds; backend, path, name="τxb")
τybt = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, time_seconds; backend, path, name="τyb")
JTbt = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, time_seconds; backend, path, name="JTb")
JSbt = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, time_seconds; backend, path, name="JSb")

Fut = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="Fu")
Fvt = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="Fv")
FTt = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="FT")
FSt = FieldTimeSeries{Nothing, Nothing, Center}(grid, time_seconds; backend, path, name="FS")
=#

# set_from_data!(Ut, U)
# set_from_data!(Vt, V)
# set_from_data!(Tt, T)
# set_from_data!(St, S)
# set_from_data!(Rit, Ri)
# set_from_data!(N²t, N²)
# set_from_data!(S²t, S²)
# set_from_data!(wut, wu)
# set_from_data!(wvt, wv)
# set_from_data!(wTt, wT)
# set_from_data!(wSt, wS)

data = (U, V, T, S, Ri, N², S², wu, wv, wT, wS)
fts  = (Ut, Vt, Tt, St, Rit, N²t, S²t, wut, wvt, wTt, wSt)
set_from_data!(fts, data)

