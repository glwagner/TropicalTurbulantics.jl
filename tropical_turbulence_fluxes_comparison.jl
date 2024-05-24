using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using CairoMakie
using Statistics
using JLD2
using NCDatasets
using Dates
using LaTeXStrings

lesfilename = "tropical_turbulence_whitt2022_0N140W.jld2"
scmfilename = "single_column_tropical_turbulence_tiny_time_step_Nz108.jld2"

glsfilename = "19851002.ocean_hourly_GLS_0.nc"
ds = Dataset(glsfilename)
wTgls = - ds["Tflx_dia_diff"][1, 1, :, :]
wTgls = permutedims(wTgls)
zgls = - ds["zi"][:]

wTgls = reverse(wTgls, dims=2)
zgls = reverse(zgls)

tgls_dates = ds["time"][:]
tgls = datetime2unix.(tgls_dates)
tgls .-= tgls[1]
tglsdays = tgls ./ days

# wTtles = FieldTimeSeries(lesfilename, "wT")
# wTtscm = FieldTimeSeries(scmfilename, "wT_ccf")

tles = wTtles.times
tlesdays = tles ./ days

tscm = wTtscm.times
tscmdays = tscm ./ days

Ntles = length(wTtles)
Ntscm = length(wTtscm)

lesgrid = wTtles.grid
scmgrid = wTtscm.grid

Nzles = size(lesgrid, 3)
Nzscm = size(scmgrid, 3)

zles = znodes(wTtles)
zscm = znodes(wTtscm)

# Uniform spacing.
Δzles = zles[2] - zles[1]
Δzscm = zscm[2] - zscm[1]
Δzgls = zgls[2] - zgls[1]

time_series_to_array(ts) = interior(ts, 1, 1, :, :) |> permutedims |> Array
wTscm = time_series_to_array(wTtscm)
wTles = time_series_to_array(wTtles)

dz_wTscm = (wTscm[:, 2:end] .- wTscm[:, 1:end-1]) ./ Δzscm
dz_wTles = (wTles[:, 2:end] .- wTles[:, 1:end-1]) ./ Δzles
dz_wTgls = (wTgls[:, 2:end] .- wTgls[:, 1:end-1]) ./ Δzgls

z_dzles = (zles[2:end] .+ zles[1:end-1]) / 2
z_dzscm = (zscm[2:end] .+ zscm[1:end-1]) / 2
z_dzgls = (zgls[2:end] .+ zgls[1:end-1]) / 2

set_theme!(Theme(fontsize=16, font="Helvetica"))
fig = Figure(size=(1000, 1000))

zlims = (-100, 0)
tlims = (7.6, 13.8)
limits = (tlims, zlims)

ax_wT_les = Axis(fig[1, 1]; limits, xlabel="Time (days)", ylabel="z (m)")
ax_wT_gls = Axis(fig[3, 1]; limits, xlabel="Time (days)", ylabel="z (m)")
ax_wT_scm = Axis(fig[2, 1]; limits, xlabel="Time (days)", ylabel="z (m)")

ax_dz_wT_les = Axis(fig[4, 1]; limits, xlabel="Time (days)", ylabel="z (m)")
ax_dz_wT_gls = Axis(fig[6, 1]; limits, xlabel="Time (days)", ylabel="z (m)")
ax_dz_wT_scm = Axis(fig[5, 1]; limits, xlabel="Time (days)", ylabel="z (m)")

colorrange = (-1e-4, 1e-4)
colormap = :balance
hm = heatmap!(ax_wT_les, tlesdays, zles, wTles; colormap, colorrange)
hm = heatmap!(ax_wT_gls, tglsdays, zgls, wTgls; colormap, colorrange)
hm = heatmap!(ax_wT_scm, tscmdays, zscm, wTscm; colormap, colorrange)
Colorbar(fig[1:3, 2], hm, label=L"\overline{wT} \mathrm{(m ᵒC s⁻¹)}")

colorrange = (-1e-5, 1e-5)
hm = heatmap!(ax_dz_wT_les, tlesdays, z_dzles, dz_wTles; colormap, colorrange)
hm = heatmap!(ax_dz_wT_gls, tglsdays, z_dzles, dz_wTgls; colormap, colorrange)
hm = heatmap!(ax_dz_wT_scm, tscmdays, z_dzles, dz_wTscm; colormap, colorrange)
Colorbar(fig[4:6, 2], hm, label=L"\partial_z \overline{wT} \mathrm{(ᵒC s⁻¹)}")

tt = tlims[1] + 0.01 * (tlims[end] - tlims[1])
zt = -95
align = (:left, :bottom)
text!(ax_wT_les, tt, zt; align, text="(a) LES (Whitt et al. 2022)")
text!(ax_wT_scm, tt, zt; align, text="(b) CATKE")
text!(ax_wT_gls, tt, zt; align, text="(c) GLS (Reichl et al. 2024)")

text!(ax_dz_wT_les, tt, zt; align, text="(d) LES (Whitt et al. 2022)")
text!(ax_dz_wT_scm, tt, zt; align, text="(e) CATKE")
text!(ax_dz_wT_gls, tt, zt; align, text="(f) GLS (Reichl et al 2024)")

display(fig)

#save("compare_les_scm_fluxes.pdf", fig)
save("compare_les_scm_fluxes.png", fig)

