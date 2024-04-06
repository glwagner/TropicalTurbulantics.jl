using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using GLMakie
using JLD2

lesfilename = "tropical_turbulence_Nz54_averages.jld2"
whifilename = "tropical_turbulence_original_les_data.jld2"
#scmfilename = "tropical_turbulence_single_column_model_tiny_time_step_Nz54.jld2"
scmfilename = "tropical_turbulence_single_column_model_Nz54.jld2"

#scmfilename = "tropical_turbulence_single_column_model_Nz108.jld2"
#scmfilename = "tropical_turbulence_single_column_model_non_conservative_Nz108.jld2"
#scmfilename = "tropical_turbulence_single_column_model_Nz216.jld2"
#scmfilename = "tropical_turbulence_single_column_model_Nz54.jld2"
#scmfilename = "tropical_turbulence_single_column_model_Nz27.jld2"

utles = FieldTimeSeries(lesfilename, "u")
vtles = FieldTimeSeries(lesfilename, "v")
Ttles = FieldTimeSeries(lesfilename, "T")
Stles = FieldTimeSeries(lesfilename, "S")
etles = FieldTimeSeries(lesfilename, "e")
Rtles = FieldTimeSeries(lesfilename, "Ri")

utwhi = FieldTimeSeries(whifilename, "u")
vtwhi = FieldTimeSeries(whifilename, "v")
Ttwhi = FieldTimeSeries(whifilename, "T")
Stwhi = FieldTimeSeries(whifilename, "S")

utscm = FieldTimeSeries(scmfilename, "u")
vtscm = FieldTimeSeries(scmfilename, "v")
Ttscm = FieldTimeSeries(scmfilename, "T")
Stscm = FieldTimeSeries(scmfilename, "S")
etscm = FieldTimeSeries(scmfilename, "e")
Rtscm = FieldTimeSeries(scmfilename, "Ri")

Nt = length(utscm)

n = Observable(1)
unles = @lift interior(utles[$n], 1, 1, :)
Tnles = @lift interior(Ttles[$n], 1, 1, :)
Snles = @lift interior(Stles[$n], 1, 1, :)
enles = @lift interior(etles[$n], 1, 1, :)

unwhi = @lift interior(utwhi[$n], 1, 1, :)
Tnwhi = @lift interior(Ttwhi[$n], 1, 1, :)
Snwhi = @lift interior(Stwhi[$n], 1, 1, :)

unscm = @lift interior(utscm[$n], 1, 1, :)
Tnscm = @lift interior(Ttscm[$n], 1, 1, :)
Snscm = @lift interior(Stscm[$n], 1, 1, :)
enscm = @lift interior(etscm[$n], 1, 1, :)

Ri⁻¹scm = 1 ./ interior(Rtscm, 1, 1, :, 1:Nt)
Ri⁻¹les = 1 ./ interior(Rtles, 1, 1, :, 1:Nt)

#Ri⁻¹les = 0 * Ri⁻¹scm

lesgrid = utles.grid
uzles = Field{Nothing, Nothing, Face}(lesgrid)
vzles = Field{Nothing, Nothing, Face}(lesgrid)
Tzles = Field{Nothing, Nothing, Face}(lesgrid)
Nzles = size(uzles, 3)

Uznles = @lift begin
    uzles .= ∂z(utles[$n])
    vzles .= ∂z(vtles[$n])
    uzi = interior(uzles, 1, 1, 2:Nzles-1)
    vzi = interior(vzles, 1, 1, 2:Nzles-1)
    return @. sqrt(uzi^2 + vzi^2)
end

Tznles = @lift begin
    Tzles .= ∂z(Ttles[$n])
    return Array(interior(Tzles, 1, 1, 2:Nzles-1))
end

Ri⁻¹les = zeros(Nzles, Nt)
Ri⁻¹ = Field{Nothing, Nothing, Face}(lesgrid)

jld2filename = "forcing_and_bcs_and_ics_0N140W.jld2"
file = jldopen(jld2filename)
α = file["thermal_expansion"]
β = file["haline_contraction"]
close(file)
g = 9.81

for n = 1:Nt
    Un = utles[n]
    Vn = vtles[n]
    Tn = Ttles[n]
    Sn = Stles[n]

    N² = g * (α * ∂z(Tn) - β * ∂z(Sn))
    S² = ∂z(Un)^2 + ∂z(Vn)^2
    Ri⁻¹ .= S² / N²
    fill_halo_regions!(Ri⁻¹)
    Ri⁻¹les[:, n] .= interior(Ri⁻¹, 1, 1, :)
end

scmgrid = utscm.grid
uzscm = Field{Nothing, Nothing, Face}(scmgrid)
vzscm = Field{Nothing, Nothing, Face}(scmgrid)
Tzscm = Field{Nothing, Nothing, Face}(scmgrid)
Nzscm = size(uzscm, 3)

Uznscm = @lift begin
    uzscm .= ∂z(utscm[$n])
    vzscm .= ∂z(vtscm[$n])
    uzi = interior(uzscm, 1, 1, 2:Nzscm-1)
    vzi = interior(vzscm, 1, 1, 2:Nzscm-1)
    return @. sqrt(uzi^2 + vzi^2)
end

Tznscm = @lift begin
    Tzscm .= ∂z(Ttscm[$n])
    return Array(interior(Tzscm, 1, 1, 2:Nzscm-1))
end

fig = Figure(size=(1600, 800))

ax_u = Axis(fig[1, 1], xlabel="u (m s⁻¹)", ylabel="z (m)")
ax_T = Axis(fig[1, 2], xlabel="T (ᵒC)", ylabel="z (m)")
ax_S = Axis(fig[2, 1], xlabel="S (g kg⁻¹)", ylabel="z (m)")
ax_e = Axis(fig[2, 2], xlabel="e (m² s⁻²)", ylabel="z (m)")

ax_uz = Axis(fig[1, 3], xlabel="∂z u (s⁻¹)", ylabel="z (m)")
ax_Tz = Axis(fig[2, 3], xlabel="∂z T (ᵒC m⁻¹)", ylabel="z (m)")

xlims!(ax_u, -2, 2)
xlims!(ax_T, 20, 31)
xlims!(ax_S, 34.8, 35.5)
xlims!(ax_e, -1e-5, 5e-4)
xlims!(ax_uz, -0.05, 0.08)
xlims!(ax_Tz, -0.05, 0.3)

zlims = (-108, 0)
tmax = 34

limits = @lift begin
    tn = t[$n] / days
    t1 = max(0, tn - 0.5)
    t2 = min(tn + 6, tmax)
    ((t1, t2), zlims)
end

ax_Rles = Axis(fig[1, 4:6], xlabel="Time (days)", ylabel="z (m)"; limits)
ax_Rscm = Axis(fig[2, 4:6], xlabel="Time (days)", ylabel="z (m)"; limits)

t = utles.times
zles = znodes(utles)
zwhi = znodes(utwhi)
zscm = znodes(utscm)
zles_Ri = znodes(Rtles)
zscm_Ri = znodes(Rtscm)

lines!(ax_u, unscm, zscm, label="CATKE")
lines!(ax_T, Tnscm, zscm, label="CATKE")
lines!(ax_S, Snscm, zscm, label="CATKE")
lines!(ax_e, enscm, zscm, label="CATKE")

lines!(ax_u, unles, zles, label="Oceananigans LES")
lines!(ax_T, Tnles, zles, label="Oceananigans LES")
lines!(ax_S, Snles, zles, label="Oceananigans LES")
lines!(ax_e, enles, zles, label="Oceananigans LES")

lines!(ax_u, unwhi, zwhi, label="Whitt et al 2022")
lines!(ax_T, Tnwhi, zwhi, label="Whitt et al 2022")
lines!(ax_S, Snwhi, zwhi, label="Whitt et al 2022")

zfscm = znodes(uzscm)
lines!(ax_uz, Uznscm, zfscm[2:Nzscm-1], label="CATKE")
lines!(ax_Tz, Tznscm, zfscm[2:Nzscm-1], label="CATKE")

zfles = znodes(uzles)
lines!(ax_uz, Uznles, zfles[2:Nzles-1], label="LES")
lines!(ax_Tz, Tznles, zfles[2:Nzles-1], label="LES")

Legend(fig[1:2, 0], ax_u)

# axislegend(ax_u)
# axislegend(ax_T, position=:rb)
# axislegend(ax_S)
# axislegend(ax_e, position=:rb)

colorrange = (0, 4.5)
colormap = :viridis
hm_R = heatmap!(ax_Rles, t ./ days, zles_Ri, permutedims(Ri⁻¹scm); colorrange, colormap)
hm_R = heatmap!(ax_Rscm, t ./ days, zscm_Ri, permutedims(Ri⁻¹les); colorrange, colormap)

tndays = @lift t[$n] / days
vlines!(ax_Rles, tndays, color=:black)
vlines!(ax_Rscm, tndays, color=:black)

Colorbar(fig[1:2, 7], hm_R, label="1 / Ri")

title = @lift "Tropical turbulence at " * prettytime(t[$n])
Label(fig[0, 0:7], title)

display(fig)

record(fig, "tropical_turbulence_catke_les_comparison.mp4", 1:10:Nt, framerate=12) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end

