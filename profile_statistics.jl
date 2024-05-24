using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using CairoMakie
#using GLMakie
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

N²tles = FieldTimeSeries(lesfilename, "N²")
N²tscm = FieldTimeSeries(scmfilename, "N²")
S²tles = FieldTimeSeries(lesfilename, "S²")
S²tscm = FieldTimeSeries(scmfilename, "S²")
Ritles = FieldTimeSeries(lesfilename, "Ri")
Ritscm = FieldTimeSeries(scmfilename, "Ri")

average_out!(Jᵀ, 0)

tles = Ritles.times
tlesdays = tles ./ days

tscm = Ritscm.times
tscmdays = tscm ./ days

Ntles = length(Ritles)
Ntscm = length(Ritscm)

lesgrid = Ritles.grid
scmgrid = Ritscm.grid

Nzles = size(lesgrid, 3)
Nzscm = size(scmgrid, 3)

zles = znodes(Ritles)
zscm = znodes(Ritscm)

# Uniform spacing.
Δzles = zles[2] - zles[1]
Δzscm = zscm[2] - zscm[1]

time_series_to_array(ts) = interior(ts, 1, 1, :, :) |> permutedims |> Array

Riscm = time_series_to_array(Ritscm)
Riles = time_series_to_array(Ritles)
N²scm = time_series_to_array(N²tscm)
N²les = time_series_to_array(N²tles)
S²scm = time_series_to_array(S²tscm)
S²les = time_series_to_array(S²tles)

Riscm_med = median(Riscm, dims=1)[:]   
Riles_med = median(Riles, dims=1)[:]
N²scm_med = median(N²scm, dims=1)[:]
N²les_med = median(N²les, dims=1)[:]
S²scm_med = median(S²scm, dims=1)[:]
S²les_med = median(S²les, dims=1)[:]

Riscm_avg = mean(Riscm, dims=1)[:] 
Riles_avg = mean(Riles, dims=1)[:]
N²scm_avg = mean(N²scm, dims=1)[:]
N²les_avg = mean(N²les, dims=1)[:]
S²scm_avg = mean(S²scm, dims=1)[:]
S²les_avg = mean(S²les, dims=1)[:]

Ntscm, Nzscm = size(Riscm)
Ntles, Nzles = size(Riles)

Riscm_q1 = zeros(Nzscm)
N²scm_q1 = zeros(Nzscm)
S²scm_q1 = zeros(Nzscm)

Riscm_q3 = zeros(Nzscm)
N²scm_q3 = zeros(Nzscm)
S²scm_q3 = zeros(Nzscm)

Riles_q1 = zeros(Nzles)
N²les_q1 = zeros(Nzles)
S²les_q1 = zeros(Nzles)

Riles_q3 = zeros(Nzles)
N²les_q3 = zeros(Nzles)
S²les_q3 = zeros(Nzles)

q1 = 0.35
q3 = 0.65

for k = 1:Nzscm
    Riscm_q1[k] = quantile(view(Riscm, :, k), q1)
    N²scm_q1[k] = quantile(view(N²scm, :, k), q1)
    S²scm_q1[k] = quantile(view(S²scm, :, k), q1)

    Riscm_q3[k] = quantile(view(Riscm, :, k), q3)
    N²scm_q3[k] = quantile(view(N²scm, :, k), q3)
    S²scm_q3[k] = quantile(view(S²scm, :, k), q3)
end

for k = 2:Nzles-2
    Riles_q1[k] = quantile(view(Riles, 2:Ntles, k), q1)
    N²les_q1[k] = quantile(view(N²les, 2:Ntles, k), q1)
    S²les_q1[k] = quantile(view(S²les, 2:Ntles, k), q1)

    Riles_q3[k] = quantile(view(Riles, 2:Ntles, k), q3)
    N²les_q3[k] = quantile(view(N²les, 2:Ntles, k), q3)
    S²les_q3[k] = quantile(view(S²les, 2:Ntles, k), q3)
end

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(size=(1200, 600))

axRi = Axis(fig[1, 1], ylabel="z (m)", xlabel="Ri")
axN² = Axis(fig[1, 2], ylabel="z (m)", xlabel="N²")
axS² = Axis(fig[1, 3], ylabel="z (m)", xlabel="S²", yaxisposition=:right)

lescolor = :black
scmcolor = :royalblue
αles = 0.4
αscm = 0.3
d = 4.0
densedash = Linestyle([0.0, d, 1.7d, 2.7d])

lines!(axRi, Riles_med, zles, color=lescolor, label="LES")
lines!(axRi, Riscm_med, zscm, color=scmcolor, label="CATKE")

Legend(fig[0, 1:3], axRi, nbanks=2)

band!(axRi, Point2f.(Riles_q1, zles), Point2f.(Riles_q3, zles), color=(lescolor, αles))
band!(axRi, Point2f.(Riscm_q1, zscm), Point2f.(Riscm_q3, zscm), color=(scmcolor, αscm))

#lines!(axRi, Riles_avg, zles, color=lescolor, linestyle=densedash)
#lines!(axRi, Riscm_avg, zscm, color=scmcolor, linestyle=densedash)

lines!(axN², N²scm_med, zscm, color=scmcolor)
lines!(axN², N²les_med, zles, color=lescolor)

band!(axN², Point2f.(N²les_q1, zles), Point2f.(N²les_q3, zles), color=(lescolor, αles))
band!(axN², Point2f.(N²scm_q1, zscm), Point2f.(N²scm_q3, zscm), color=(scmcolor, αscm))

#lines!(axN², N²les_avg, zles, color=lescolor, linestyle=densedash)
#lines!(axN², N²scm_avg, zscm, color=scmcolor, linestyle=densedash)

lines!(axS², S²scm_med, zscm, color=scmcolor)
lines!(axS², S²les_med, zles, color=lescolor)

band!(axS², Point2f.(S²les_q1, zles), Point2f.(S²les_q3, zles), color=(lescolor, αles))
band!(axS², Point2f.(S²scm_q1, zscm), Point2f.(S²scm_q3, zscm), color=(scmcolor, αscm))

#lines!(axS², S²les_avg, zles, color=lescolor, linestyle=densedash)
#lines!(axS², S²scm_avg, zscm, color=scmcolor, linestyle=densedash)

xlims!(axRi, 0.17, 0.35)
xlims!(axN², -5e-5, 2.9e-4)
xlims!(axS², 0, 1.1e-3)

ylims!(axRi, -90, 0)
ylims!(axN², -90, 0)
ylims!(axS², -90, 0)

hidespines!(axRi, :t, :r)
hidespines!(axN², :t, :r, :l)
hidespines!(axS², :t, :l)

hideydecorations!(axN², grid=false)

save("tropical_turbulence_profile_statistics.pdf", fig)

display(fig)

