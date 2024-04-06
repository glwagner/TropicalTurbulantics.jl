using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using GLMakie
using Statistics
using JLD2

lesfilename = "tropical_turbulence_Nz54_averages.jld2"
#scmfilename = "tropical_turbulence_single_column_model_Nz108.jld2"
scmfilename = "tropical_turbulence_single_column_model_tiny_time_step_Nz54.jld2"

utles = FieldTimeSeries(lesfilename, "u")
vtles = FieldTimeSeries(lesfilename, "v")
Ttles = FieldTimeSeries(lesfilename, "T")
etles = FieldTimeSeries(lesfilename, "e")
Rtles = FieldTimeSeries(lesfilename, "Ri")

utscm = FieldTimeSeries(scmfilename, "u")
vtscm = FieldTimeSeries(scmfilename, "v")
Ttscm = FieldTimeSeries(scmfilename, "T")
etscm = FieldTimeSeries(scmfilename, "e")
Rtscm = FieldTimeSeries(scmfilename, "Ri")

uscm = Array(interior(utscm, 1, 1, :, 1:Nt))
ules = Array(interior(utles, 1, 1, :, 1:Nt))
Tscm = Array(interior(Ttscm, 1, 1, :, 1:Nt))
Tles = Array(interior(Ttles, 1, 1, :, 1:Nt))
escm = Array(interior(etscm, 1, 1, :, 1:Nt))
eles = Array(interior(etles, 1, 1, :, 1:Nt))

Ri⁻¹scm = 1 ./ interior(Rtscm, 1, 1, :, 1:Nt)

lesgrid = utles.grid
scmgrid = utscm.grid
Nt = length(utscm)
Nzles = size(lesgrid, 3)
Nzscm = size(scmgrid, 3)

Tzles = zeros(Nzles+1, Nt)
Tzscm = zeros(Nzscm+1, Nt)
uzles = zeros(Nzles+1, Nt)
uzscm = zeros(Nzscm+1, Nt)
Ri⁻¹les = zeros(Nzles+1, Nt)
ψles = Field{Nothing, Nothing, Face}(lesgrid)
ψscm = Field{Nothing, Nothing, Face}(scmgrid)

jld2filename = "forcing_and_bcs_and_ics_0N140W.jld2"
file = jldopen(jld2filename)
α = file["thermal_expansion"]
β = file["haline_contraction"]
close(file)
g = 9.81

for n = 1:Nt
    un = utles[n]
    vn = vtles[n]
    Tn = Ttles[n]
    Sn = Stles[n]

    ψles .= ∂z(Tn)
    fill_halo_regions!(ψles)
    Tzles[:, n] .= interior(ψles, 1, 1, :)

    ψles .= ∂z(un)
    fill_halo_regions!(ψles)
    uzles[:, n] .= interior(ψles, 1, 1, :)

    ψscm .= ∂z(Ttscm[n])
    fill_halo_regions!(ψscm)
    Tzscm[:, n] .= interior(ψscm, 1, 1, :)

    ψscm .= ∂z(utscm[n])
    fill_halo_regions!(ψscm)
    uzscm[:, n] .= interior(ψscm, 1, 1, :)

    N² = g * (α * ∂z(Tn) - β * ∂z(Sn))
    S² = ∂z(un)^2 + ∂z(vn)^2
    ψles .= S² / N²
    fill_halo_regions!(ψles)
    Ri⁻¹les[:, n] .= interior(ψles, 1, 1, :)
end

Ri⁻¹les[Ri⁻¹les .< 0] .= NaN
Ri⁻¹scm[Ri⁻¹scm .< 0] .= NaN

fig = Figure(size=(1000, 1000))

ax_eles = Axis(fig[1, 1], xlabel="Time (days)", ylabel="z (m)", xaxisposition=:top)
ax_escm = Axis(fig[2, 1], xlabel="Time (days)", ylabel="z (m)")

ax_Tzles = Axis(fig[3, 1], xlabel="Time (days)", ylabel="z (m)")
ax_Tzscm = Axis(fig[4, 1], xlabel="Time (days)", ylabel="z (m)")

ax_Rles = Axis(fig[5, 1], xlabel="Time (days)", ylabel="z (m)")
ax_Rscm = Axis(fig[6, 1], xlabel="Time (days)", ylabel="z (m)")

hidexdecorations!(ax_escm)
hidexdecorations!(ax_Tzles)
hidexdecorations!(ax_Tzscm)
hidexdecorations!(ax_Rles)

t = utles.times
tdays = t ./ days

zles_u  = znodes(utles)
zscm_u  = znodes(utscm)
zles_Ri = znodes(Rtles)
zscm_Ri = znodes(Rtscm)

colorrange = (0, 0.7)
colormap = :solar
hm = heatmap!(ax_eles, tdays, zscm_u, permutedims(escm) ./ maximum(escm); colorrange, colormap)
hm = heatmap!(ax_escm, tdays, zles_u, permutedims(eles) ./ maximum(eles); colorrange, colormap)
Colorbar(fig[1:2, 2], hm, label="e / max(e)")

colorrange = (-0.01, 0.01)
colormap = :balance
hm = heatmap!(ax_Tzles, tdays, zles_Ri, permutedims(Tzles); colorrange, colormap)
hm = heatmap!(ax_Tzscm, tdays, zscm_Ri, permutedims(Tzscm); colorrange, colormap)
Colorbar(fig[3:4, 2], hm, label="∂z T (ᵒC m⁻¹)")

colorrange = (0, 5)
colormap = :viridis
levels = 0:0.5:5
hm = heatmap!(ax_Rles, tdays, zscm_Ri, permutedims(Ri⁻¹scm); colorrange, colormap, nan_color=:gray) #, levels)
hm = heatmap!(ax_Rscm, tdays, zles_Ri, permutedims(Ri⁻¹les); colorrange, colormap, nan_color=:gray) #, levels)
Colorbar(fig[5:6, 2], hm, label="1 / Ri = (∂z u)² / N²")

tlims = (6.5, 11.0)
xlims!(ax_Tzles, tlims...)
xlims!(ax_Tzscm, tlims...)
xlims!(ax_uzles, tlims...)
xlims!(ax_uzscm, tlims...)
xlims!(ax_eles, tlims...)
xlims!(ax_escm, tlims...)
xlims!(ax_Rles, tlims...)
xlims!(ax_Rscm, tlims...)

zlims = (-100, 0)
ylims!(ax_Tzles, zlims...)
ylims!(ax_Tzscm, zlims...)
ylims!(ax_uzles, zlims...)
ylims!(ax_uzscm, zlims...)
ylims!(ax_eles,  zlims...)
ylims!(ax_escm,  zlims...)
ylims!(ax_Rles,  zlims...)
ylims!(ax_Rscm,  zlims...)

display(fig)

save("compare_les_scm_turbulence.png", fig)

