using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
# using CairoMakie
using GLMakie
using Statistics
using JLD2
using NCDatasets
using MathTeXEngine

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

# Load insolation
dir = "data_les"
filename = "ROMS_PSH_6HRLIN_0N140W_360x360x216_22OCT2020.nc"
filepath = joinpath(dir, filename)
dataset = Dataset(filepath)
Jᵀ = - dataset["kappadtdztop"][:]
Fᴵ = dataset["dTdtSOLAR"][:, :]
Fᵀ = dataset["dTdtFORCE"][:, :]

wrms = dataset["wrms"][:, :]
w²les = wrms.^2
wules = dataset["uw"][:, :]

close(dataset)

lesfilename = "tropical_turbulence_whitt2022_0N140W.jld2"
scmfilename = "single_column_tropical_turbulence_tiny_time_step_Nz108.jld2"

etscm = FieldTimeSeries(scmfilename, "e")
Rtles = FieldTimeSeries(lesfilename, "Ri")
Rtscm = FieldTimeSeries(scmfilename, "Ri")
wutscm = FieldTimeSeries(scmfilename, "wu")

average_out!(Jᵀ, 0)

tles = Rtles.times
tlesdays = tles ./ days

tscm = Rtscm.times
tscmdays = tscm ./ days

Ntles = length(Rtles)
Ntscm = length(Rtscm)

lesgrid = Rtles.grid
scmgrid = Rtscm.grid

Nzles = size(lesgrid, 3)
Nzscm = size(scmgrid, 3)

zles = znodes(Rtles)
zscm = znodes(Rtscm)

# Uniform spacing.
Δzles = zles[2] - zles[1]
Δzscm = zscm[2] - zscm[1]

time_series_to_array(ts) = interior(ts, 1, 1, :, :) |> permutedims |> Array

eles = permutedims(w²les)
escm = time_series_to_array(etscm)

Rscm = time_series_to_array(Rtscm)
Rles = time_series_to_array(Rtles)
wules = permutedims(wules)
wuscm = time_series_to_array(wutscm)

#set_theme!(Theme(fontsize=16, font="Helvectica"))
#set_theme!(Theme(fontsize=16, font=texfont()))
set_theme!(Theme(fontsize=16))
fig = Figure(size=(1000, 800))

Qlims = (-1600, 300)
zlims = (-100, 0)

#=
slider = Slider(fig[7, 1], range=1:0.1:34, startvalue=0)
n = slider.value

limits = @lift begin
    t1 = max(0,  $n - 0.5)
    t2 = min(34, $n + 5.5)
    ((t1, t2), zlims)
end

Qtlimits = @lift begin
    t1 = max(0,  $n - 0.5)
    t2 = min(34, $n + 5.5)
    ((t1, t2), Qlims)
end
=#

tlims = (7.6, 13.8)
limits = (tlims, zlims)
Qtlimits = (tlims, Qlims)

#ax_Q = Axis(fig[2, 1]; limits=Qtlimits, xlabel="Time (days)", ylabel=L"Heat flux ($\mathrm{W \, m^{-2}}$)", xaxisposition=:top)
ax_Q = Axis(fig[2, 1]; limits=Qtlimits, xlabel="Time (days)", ylabel="Heat flux (W m⁻²)", xaxisposition=:top)

ax_e_les = Axis(fig[3, 1]; limits, xlabel="Time (days)", ylabel="z (m)")
ax_e_scm = Axis(fig[4, 1]; limits, xlabel="Time (days)", ylabel="z (m)")

ax_R_les = Axis(fig[5, 1]; limits, xlabel="Time (days)", ylabel="z (m)")
ax_R_scm = Axis(fig[6, 1]; limits, xlabel="Time (days)", ylabel="z (m)")

# ax_τ_les = Axis(fig[7, 1]; limits, xlabel="Time (days)", ylabel="z (m)")
# ax_τ_scm = Axis(fig[8, 1]; limits, xlabel="Time (days)", ylabel="z (m)")

ρₒ = 1024
cₚ = 3991

∫Fᴵ = sum(Fᴵ[1:end-1, :], dims=1)[:]
∫Fᵀ = sum(Fᵀ[1:end-1, :], dims=1)[:]
average_out!(∫Fᴵ, 0)
average_out!(∫Fᵀ, 0)

Jᴵ = Fᴵ[end-1, :]
average_out!(Jᴵ, 0)

ΣQ = ρₒ .* cₚ .* Jᵀ
Qs = - ρₒ .* cₚ .* ∫Fᴵ
#Qs = - ρₒ .* cₚ .* Jᴵ
#Qs = - ρₒ .* cₚ .* (∫Fᴵ .+ ∫Fᵀ)

d = 4.0
densedash = Linestyle([0.0, d, 1.7d, 2.7d])
lines!(ax_Q, tlesdays, ΣQ .+ Qs,                color=(:black, 0.6),     linewidth=4, label="Surface heat flux + column-integrated solar insolation") 
lines!(ax_Q, tlesdays, Qs, linestyle=densedash, color=(:tomato1, 0.9),   linewidth=2, label="Column-integrated solar insolation")
lines!(ax_Q, tlesdays, ΣQ, linestyle=densedash, color=(:royalblue, 0.9), linewidth=2, label="Surface heat flux (sensible, latent, and longwave radiation)")

hlines!(ax_Q, 0, color=:black, linewidth=1)
Legend(fig[1, 1], ax_Q, tellwidth=false, nbanks=1, framevisible=false)

colorrange = (0, 2e-4)
colormap = :solar
hm = heatmap!(ax_e_les, tlesdays, zles, eles; colormap, colorrange)
#Colorbar(fig[3, 2], hm, label="LES w² \n (m² s⁻²)")

colorrange = (0, 2e-4)
hm = heatmap!(ax_e_scm, tscmdays, zscm, escm; colormap, colorrange)
ticks = ([0, 1e-4, 2e-4], ["0", "10⁻⁴", "2 × 10⁻⁴"])
#label = L"CATKE $e$ and LES $\overline{w′²}$ \n (both $\mathrm{m² \, s^{-2}}$)"
label = "CATKE e and LES w′² \n (both m² s⁻²)"
Colorbar(fig[3:4, 2], hm; ticks, label)

Rles[Rles .< 0] .= NaN
Rscm[Rscm .< 0] .= NaN
colorrange = (0.15, 0.3)
colormap = :tempo
hm = heatmap!(ax_R_les, tlesdays, zles, Rles; colormap, colorrange, nan_color=:red)
hm = heatmap!(ax_R_scm, tscmdays, zscm, Rscm; colormap, colorrange, nan_color=:red)
#Colorbar(fig[5:6, 2], hm, label=L"Ri = \partial_z b / \left (\partial_z u \right )^2")
Colorbar(fig[5:6, 2], hm, label="Richardson number")

#=
colorrange = (-1e-4, 1e-4)
colormap = :balance
hm = heatmap!(ax_τ_les, tlesdays, zles, wules; colormap, colorrange, nan_color=:red)
hm = heatmap!(ax_τ_scm, tscmdays, zscm, wuscm; colormap, colorrange, nan_color=:red)
Colorbar(fig[7:8, 2], hm, label=L"\overline{uw}")
=#

#=
xlims!(ax_e_les, tlims...)
ylims!(ax_e_les, zlims...)
xlims!(ax_e_scm, tlims...)
ylims!(ax_e_scm, zlims...)

xlims!(ax_R_les, tlims...)
ylims!(ax_R_les, zlims...)
xlims!(ax_R_scm, tlims...)
ylims!(ax_R_scm, zlims...)
=#

tt = tlims[1] + 0.01 * (tlims[end] - tlims[1])
zt = -95
align = (:left, :bottom)
text!(ax_e_les, tt, zt; align, color=:white, text="(a) LES (Whitt et al. 2022)")
text!(ax_e_scm, tt, zt; align, color=:white, text="(b) CATKE")
text!(ax_R_les, tt, zt; align, color=:white, text="(c) LES (Whitt et al. 2022)")
text!(ax_R_scm, tt, zt; align, color=:white, text="(d) CATKE")
#text!(ax_τ_les, tt, zt; align, color=:black, text="LES (Whitt et al. 2022)")
#text!(ax_τ_scm, tt, zt; align, color=:black, text="CATKE")

rowsize!(fig.layout, 1, Relative(0.15))
rowsize!(fig.layout, 2, Relative(0.1))

hidexdecorations!(ax_e_les)
hidexdecorations!(ax_e_scm)
hidexdecorations!(ax_R_les)
#hidexdecorations!(ax_R_scm)
#hidexdecorations!(ax_τ_les)

display(fig)

save("tropical_turbulence_overview.png", fig)

