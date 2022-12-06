using GLMakie
using Oceananigans

filename = "tropical_turbulence_single_column_model.jld2"

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

fig = Figure()

ax_u = Axis(fig[1, 1], xlabel="Time (days)", ylabel="z (m)")
ax_v = Axis(fig[1, 2], xlabel="Time (days)", ylabel="z (m)")
ax_T = Axis(fig[2, 1], xlabel="Time (days)", ylabel="z (m)")
ax_S = Axis(fig[2, 2], xlabel="Time (days)", ylabel="z (m)")
ax_e = Axis(fig[3, 1], xlabel="Time (days)", ylabel="z (m)")
ax_Ri = Axis(fig[3, 2], xlabel="Time (days)", ylabel="z (m)")

heatmap!(ax_u,  t ./ days, z, permutedims(u))
heatmap!(ax_v,  t ./ days, z, permutedims(v))
heatmap!(ax_T,  t ./ days, z, permutedims(T))
heatmap!(ax_S,  t ./ days, z, permutedims(S))
heatmap!(ax_e,  t ./ days, z, permutedims(e), colorrange=(0, 1e-4))
heatmap!(ax_Ri, t ./ days, z_Ri, permutedims(Ri))

display(fig)

