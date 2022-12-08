using GLMakie
using Oceananigans

filename = "tropical_turbulence_single_column_model_Nz54.jld2"

ut = FieldTimeSeries(filename, "u")
vt = FieldTimeSeries(filename, "v")
Tt = FieldTimeSeries(filename, "T")
St = FieldTimeSeries(filename, "S")
et = FieldTimeSeries(filename, "e")
Rit = FieldTimeSeries(filename, "Ri")

t = ut.times
z = znodes(ut)
z_Ri = znodes(Rit)

u = interior(ut, 1, 1, :, :)
v = interior(vt, 1, 1, :, :)
T = interior(Tt, 1, 1, :, :)
S = interior(St, 1, 1, :, :)
e = interior(et, 1, 1, :, :)
Ri = interior(Rit, 1, 1, :, :)
Ri⁻¹ = 1 ./ interior(Rit, 1, 1, :, :)

umax = max(maximum(abs, u), maximum(abs, v))
ulim = 3umax / 4

fig = Figure(resolution=(2400, 1200))

ax_u = Axis(fig[1, 1], xlabel="Time (days)", ylabel="z (m)")
ax_v = Axis(fig[1, 2], xlabel="Time (days)", ylabel="z (m)")
ax_T = Axis(fig[2, 1], xlabel="Time (days)", ylabel="z (m)")
ax_S = Axis(fig[2, 2], xlabel="Time (days)", ylabel="z (m)")
ax_e = Axis(fig[3, 1], xlabel="Time (days)", ylabel="z (m)")
ax_R = Axis(fig[3, 2], xlabel="Time (days)", ylabel="z (m)")

hm_u = heatmap!(ax_u, t ./ days, z, permutedims(u), colormap=:balance, colorrange=(-1, 1))
hm_v = heatmap!(ax_v, t ./ days, z, permutedims(v), colormap=:balance, colorrange=(-1, 1))
hm_T = heatmap!(ax_T, t ./ days, z, permutedims(T), colormap=:thermal, colorrange=(20, 26.1))
hm_S = heatmap!(ax_S, t ./ days, z, permutedims(S), colormap=:haline, colorrange=(34.9, 35.3))
hm_e = heatmap!(ax_e, t ./ days, z, permutedims(e), colormap=:solar, colorrange=(0, 1e-4))
hm_R = heatmap!(ax_R, t ./ days, z_Ri, permutedims(Ri⁻¹), colorrange=(0, 5))

Colorbar(fig[1, 0], hm_u, flipaxis=false, label="Zonal velocity (m s⁻¹)")
Colorbar(fig[1, 3], hm_v, label="Meridional velocity (m s⁻¹)")
Colorbar(fig[2, 0], hm_T, flipaxis=false, label="Temperature (ᵒC)")
Colorbar(fig[2, 3], hm_S, label="Salinity (g kg⁻¹)")
Colorbar(fig[3, 0], hm_e, flipaxis=false, label="Turbulent kinetic energy (m s⁻¹)")
Colorbar(fig[3, 3], hm_R, label="1 / Ri")

display(fig)

save("tropical_turbulence_single_column_simulation.png", fig)

