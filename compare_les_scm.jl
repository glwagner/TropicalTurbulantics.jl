using Oceananigans
using Oceananigans.Units
using GLMakie

lesfilename = "tropical_turbulence_Nz108_averages.jld2"
#scmfilename = "tropical_turbulence_single_column_model_Nz108.jld2"
#scmfilename = "tropical_turbulence_single_column_model_Nz54.jld2"
scmfilename = "tropical_turbulence_single_column_model_Nz27.jld2"

utles = FieldTimeSeries(lesfilename, "u")
vtles = FieldTimeSeries(lesfilename, "v")
Ttles = FieldTimeSeries(lesfilename, "T")
Stles = FieldTimeSeries(lesfilename, "S")
etles = FieldTimeSeries(lesfilename, "e")
Rtles = FieldTimeSeries(lesfilename, "Ri")

utscm = FieldTimeSeries(scmfilename, "u")
vtscm = FieldTimeSeries(scmfilename, "v")
Ttscm = FieldTimeSeries(scmfilename, "T")
Stscm = FieldTimeSeries(scmfilename, "S")
etscm = FieldTimeSeries(scmfilename, "e")
Rtscm = FieldTimeSeries(scmfilename, "Ri")

n = Observable(1)
unles = @lift interior(utles[$n], 1, 1, :)
Tnles = @lift interior(Ttles[$n], 1, 1, :)
Snles = @lift interior(Stles[$n], 1, 1, :)
enles = @lift interior(etles[$n], 1, 1, :)

unscm = @lift interior(utscm[$n], 1, 1, :)
Tnscm = @lift interior(Ttscm[$n], 1, 1, :)
Snscm = @lift interior(Stscm[$n], 1, 1, :)
enscm = @lift interior(etscm[$n], 1, 1, :)

Ri⁻¹les = 1 ./ interior(Rtles, 1, 1, :, :)
Ri⁻¹scm = 1 ./ interior(Rtscm, 1, 1, :, :)

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
xlims!(ax_T, 21, 31)
xlims!(ax_S, 32, 35.3)
xlims!(ax_e, -1e-5, 5e-4)
xlims!(ax_uz, 0, 0.1)
xlims!(ax_Tz, 0, 0.2)

ax_Rles = Axis(fig[1, 4:6], xlabel="Time (days)", ylabel="z (m)")
ax_Rscm = Axis(fig[2, 4:6], xlabel="Time (days)", ylabel="z (m)")

t = utles.times
zles = znodes(utles)
zscm = znodes(utscm)
zles_Ri = znodes(Rtles)
zscm_Ri = znodes(Rtscm)

lines!(ax_u, unles, zles, label="LES")
lines!(ax_T, Tnles, zles, label="LES")
lines!(ax_S, Snles, zles, label="LES")
lines!(ax_e, enles, zles, label="LES")

lines!(ax_u, unscm, zscm, label="CATKE")
lines!(ax_T, Tnscm, zscm, label="CATKE")
lines!(ax_S, Snscm, zscm, label="CATKE")
lines!(ax_e, enscm, zscm, label="CATKE")

zfles = znodes(uzles)
lines!(ax_uz, Uznles, zfles[2:Nzles-1], label="LES")
lines!(ax_Tz, Tznles, zfles[2:Nzles-1], label="LES")

zfscm = znodes(uzscm)
lines!(ax_uz, Uznscm, zfscm[2:Nzscm-1], label="CATKE")
lines!(ax_Tz, Tznscm, zfscm[2:Nzscm-1], label="CATKE")

axislegend(ax_u)
axislegend(ax_T, position=:rb)
axislegend(ax_S)
axislegend(ax_e, position=:rb)

hm_R = heatmap!(ax_Rles, t ./ days, zles_Ri, permutedims(Ri⁻¹scm), colorrange=(0, 5))
hm_R = heatmap!(ax_Rscm, t ./ days, zscm_Ri, permutedims(Ri⁻¹les), colorrange=(0, 5))

tndays = @lift t[$n] / days
vlines!(ax_Rles, tndays, color=:black)
vlines!(ax_Rscm, tndays, color=:black)

Colorbar(fig[1:2, 7], hm_R, label="1 / Ri")

title = @lift "Tropical turbulence at " * prettytime(t[$n])
Label(fig[0, 1:7], title)

display(fig)

Nt = length(utles)
record(fig, "tropical_turbulence_catke_les_comparison.mp4", 1:10:Nt, framerate=12) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end

#save("tropical_turbulence_large_eddy_simulation.png", fig)

