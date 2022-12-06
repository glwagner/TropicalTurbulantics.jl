# Data processing for an equatorial turbulence simulation forced by 1/20ᵒ ROMS output, following
#
# Whitt, D. B. et al, "Simulation and Scaling of the Turbulence Vertical Heat Transport
# and Deep-Cycle Turbulence across the Equatorial Pacific Cold Tongue", JPO, 2022
#
# See README.md for more.
#
# Greg's questions:
#
# 1. The variable "z" extends from 0 to 107.5 m.
#    Are these cell centers, or interfaces (since the point "0" is included).
#    If interfaces, is the bottom interface missing?
#
# 2. Is "beta" defined as the "haline expansion coefficient" (rather than haline _contraction_),
#    and is therefore negative rather than positive?
#
#    More precisely, is the equation of state:
#
#    b ∼ α T + β S
#
#    or 
#
#    b ∼ α T - β S
#
# 3. It seems the only "problem" cells are the bottommost -- true? (Why?)


using NCDatasets
using GLMakie
using Oceananigans.Units
using JLD2
using Dates

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

U = dataset["ume"][:, :]
V = dataset["vme"][:, :]
T = dataset["tempme"][:, :]
S = dataset["saltme"][:, :]

# Fix bottommost row
U[216, :] .= U[215, :]
V[216, :] .= V[215, :]
T[216, :] .= T[215, :]
S[216, :] .= S[215, :]

@printf("extrema(U) = (%.2e, %.2e) \n", extrema(U)...) 
@printf("extrema(V) = (%.2e, %.2e) \n", extrema(V)...)
@printf("extrema(T) = (%.2e, %.2e) \n", extrema(T)...)
@printf("extrema(S) = (%.2e, %.2e) \n", extrema(S)...)

Fᵁ = dataset["dUdtFORCE"][:, :]
Fⱽ = dataset["dVdtFORCE"][:, :]
Fᵀ = dataset["dTdtFORCE"][:, :]
#Fˢ = dataset["dSdtFORCE"]
Fᴵ = dataset["dTdtSOLAR"][:, :]

Fᵁ[216, :] .= Fᵁ[215, :]  
Fⱽ[216, :] .= Fⱽ[215, :] 
Fᵀ[216, :] .= Fᵀ[215, :] 
Fᴵ[216, :] .= Fᴵ[215, :] 

@printf("extrema(Fᵁ) = (%.2e, %.2e) \n", extrema(Fᵁ)...) 
@printf("extrema(Fⱽ) = (%.2e, %.2e) \n", extrema(Fⱽ)...)
@printf("extrema(Fᵀ) = (%.2e, %.2e) \n", extrema(Fᵀ)...)
@printf("extrema(Fᴵ) = (%.2e, %.2e) \n", extrema(Fᴵ)...)

ρᵣ = dataset["rho0"]

Qᵁ_bottom = dataset["nududzbot"][:]
Qⱽ_bottom = dataset["nududzbot"][:]
Qᵀ_bottom = dataset["kappadsdzbot"][:]
Qˢ_bottom = dataset["kappadtdzbot"][:]

Qᵁ_surface = dataset["nududztop"][:]
Qⱽ_surface = dataset["nududztop"][:]
Qᵀ_surface = dataset["kappadsdztop"][:]
Qˢ_surface = dataset["kappadtdztop"][:]

# Fix first value
Qᵁ_surface[1] = Qᵁ_surface[2]
Qⱽ_surface[1] = Qⱽ_surface[2]
Qᵀ_surface[1] = Qᵀ_surface[2]
Qˢ_surface[1] = Qˢ_surface[2]

for Q in [Qᵁ_bottom, 
          Qⱽ_bottom, 
          Qᵀ_bottom, 
          Qˢ_bottom, 
          Qᵁ_surface, 
          Qⱽ_surface,
          Qᵀ_surface,
          Qˢ_surface]

    average_out!(Q, 0)
end

@printf("extrema(Qᵁ_surface) = (%.2e, %.2e) \n", extrema(Qᵁ_surface)...) 
@printf("extrema(Qⱽ_surface) = (%.2e, %.2e) \n", extrema(Qⱽ_surface)...)
@printf("extrema(Qᵀ_surface) = (%.2e, %.2e) \n", extrema(Qᵀ_surface)...)
@printf("extrema(Qˢ_surface) = (%.2e, %.2e) \n", extrema(Qˢ_surface)...)

@printf("extrema(Qᵁ_bottom) = (%.2e, %.2e) \n", extrema(Qᵁ_bottom)...) 
@printf("extrema(Qⱽ_bottom) = (%.2e, %.2e) \n", extrema(Qⱽ_bottom)...)
@printf("extrema(Qᵀ_bottom) = (%.2e, %.2e) \n", extrema(Qᵀ_bottom)...)
@printf("extrema(Qˢ_bottom) = (%.2e, %.2e) \n", extrema(Qˢ_bottom)...)

thermal_expansion = dataset["alpha"][:]
haline_contraction = - dataset["beta"][:]

@printf("thermal expansion = %.2e \n", thermal_expansion)
@printf("haline contraction = %.2e \n", haline_contraction)

fig = Figure(resolution=(1600, 1600))

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

axQU = Axis(fig[2, 1], xlabel="Time (hrs)", ylabel="Zonal momentum flux (m² s⁻²)")
axQV = Axis(fig[2, 2], xlabel="Time (hrs)", ylabel="Meridional momentum flux (m² s⁻²)")
axQT = Axis(fig[4, 1], xlabel="Time (hrs)", ylabel="Temperature flux (m² s⁻²)")
axQS = Axis(fig[4, 2], xlabel="Time (hrs)", ylabel="Salt flux (m² s⁻²)")

lines!(axQU, time_hours, Qᵁ_surface, label="surface")
lines!(axQV, time_hours, Qⱽ_surface, label="surface")
lines!(axQT, time_hours, Qᵀ_surface, label="surface")
lines!(axQS, time_hours, Qˢ_surface, label="surface")

lines!(axQU, time_hours, Qᵁ_bottom, label="bottom")
lines!(axQV, time_hours, Qⱽ_bottom, label="bottom")
lines!(axQT, time_hours, Qᵀ_bottom, label="bottom")
lines!(axQS, time_hours, Qˢ_bottom, label="bottom")

axislegend(axQU)

axI = Axis(fig[5, 1], xlabel="Time (hrs)", ylabel="z (m)")
hmI = heatmap!(axI, time_hours, z, permutedims(Fᴵ), colorrange=(0, 1e-5))
Colorbar(fig[5, 0], hmI, vertical=true, flipaxis=false, label="Solar insolation temperature flux divergence (ᵒC m s⁻¹)")

display(fig)

save(figname, fig)

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

file["Qᵁ_surface"] = Qᵁ_surface
file["Qⱽ_surface"] = Qⱽ_surface
file["Qᵀ_surface"] = Qᵀ_surface
file["Qˢ_surface"] = Qˢ_surface

file["Qᵁ_bottom"] = Qᵁ_bottom
file["Qⱽ_bottom"] = Qⱽ_bottom
file["Qᵀ_bottom"] = Qᵀ_bottom
file["Qˢ_bottom"] = Qˢ_bottom

file["Fᵁ"] = Fᵁ 
file["Fⱽ"] = Fⱽ 
file["Fᵀ"] = Fᵀ 
file["Fᴵ"] = Fᴵ

close(file)

