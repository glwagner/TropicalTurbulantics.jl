using NCDatasets
using GLMakie
using Oceananigans.Units
using JLD2

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
time = dataset["time"][:]
Nt = length(time)

time_seconds = datetime2unix.(time)
time_seconds .-= time_seconds[1]
time_hours = time_seconds ./ hours

U = dataset["ume"][:, :]
V = dataset["vme"][:, :]
T = dataset["tempme"][:, :]
S = dataset["saltme"][:, :]

U[216, :] .= U[215, :]
V[216, :] .= V[215, :]
T[216, :] .= T[215, :]
S[216, :] .= S[215, :]

# Fix bottommost row
@show extrema(U)
@show extrema(V)
@show extrema(T)
@show extrema(S)

Fᵁ = dataset["dUdtFORCE"][:, :]
Fⱽ = dataset["dVdtFORCE"][:, :]
Fᵀ = dataset["dTdtFORCE"][:, :]
#Fˢ = dataset["dSdtFORCE"]

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

Fᴵ = dataset["dTdtSOLAR"][:, :]
thermal_expansion = dataset["alpha"][:]
haline_contraction = dataset["beta"][:]

fig = Figure(resolution=(1600, 1600))

axU = Axis(fig[1, 1], xlabel="Time (hr)", ylabel="z (m)")
axV = Axis(fig[1, 2], xlabel="Time (hr)", ylabel="z (m)")
axT = Axis(fig[3, 1], xlabel="Time (hr)", ylabel="z (m)")
axS = Axis(fig[3, 2], xlabel="Time (hr)", ylabel="z (m)")

hmU = heatmap!(axU, time_hours, z, permutedims(U), colormap=:balance, colorrange=(-1, 1))
hmV = heatmap!(axV, time_hours, z, permutedims(V), colormap=:balance, colorrange=(-1, 1))
hmT = heatmap!(axT, time_hours, z, permutedims(T), colormap=:thermal, colorrange=(20, 26.1))
hmS = heatmap!(axS, time_hours, z, permutedims(S), colormap=:haline, colorrange=(34.9, 35.3))

Colorbar(fig[0, 1], hmU, vertical=false, label="Zonal velocity (m s⁻¹)")
Colorbar(fig[0, 2], hmV, vertical=false, label="Meridional velocity (m s⁻¹)")
Colorbar(fig[4, 1], hmT, vertical=false, label="Temperature (ᵒC)")
Colorbar(fig[4, 2], hmS, vertical=false, label="Salinity (g kg⁻¹)")

axQU = Axis(fig[2, 1], xlabel="Time (hrs)", ylabel="Zonal momentum flux (m² s⁻²)")
axQV = Axis(fig[2, 2], xlabel="Time (hrs)", ylabel="Meridional momentum flux (m² s⁻²)")
axQT = Axis(fig[5, 1], xlabel="Time (hrs)", ylabel="Temperature flux (m² s⁻²)")
axQS = Axis(fig[5, 2], xlabel="Time (hrs)", ylabel="Salt flux (m² s⁻²)")

lines!(axQU, time_hours, Qᵁ_surface)
lines!(axQV, time_hours, Qⱽ_surface)
lines!(axQT, time_hours, Qᵀ_surface)
lines!(axQS, time_hours, Qˢ_surface)

lines!(axQU, time_hours, Qᵁ_bottom)
lines!(axQV, time_hours, Qⱽ_bottom)
lines!(axQT, time_hours, Qᵀ_bottom)
lines!(axQS, time_hours, Qˢ_bottom)

axI = Axis(fig[6, 1], xlabel="Time (hrs)", ylabel="z (m)")
heatmap!(axI, time_hours, z, permutedims(insolation), colorrange=(0, 1e-5))

display(fig)

save(figname, fig)

#####
##### Save all the data
#####

rm(jld2filename, force=true)
file = jldopen(jld2filename, "a+")

file["z"] = z
file["time"] = time
file["time_seconds"] = time_seconds
file["ρᵣ"] = ρᵣ
file["thermal_expansion"]  = thermal_expansion
file["haline_contraction"] = haline_contraction

file["U"] = U
file["V"] = V
file["T"] = T
file["S"] = S

file["Qᵁ_surface"] =  Qᵁ_surface
file["Qⱽ_surface"] =  Qⱽ_surface
file["Qᵀ_surface"] =  Qᵀ_surface
file["Qˢ_surface"] =  Qˢ_surface

file["Qᵁ_bottom"]  =  Qᵁ_bottom
file["Qⱽ_bottom"]  =  Qⱽ_bottom
file["Qᵀ_bottom"]  =  Qᵀ_bottom
file["Qˢ_bottom"]  =  Qˢ_bottom

file["Fᵁ"] = Fᵁ 
file["Fⱽ"] = Fⱽ 
file["Fᵀ"] = Fᵀ 
file["Fᴵ"] = Fᴵ

close(file)

