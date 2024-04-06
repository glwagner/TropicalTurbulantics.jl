using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using GLMakie
using JLD2

lesfilename = "tropical_turbulence_Nz54_averages.jld2"
scmfilename = "tropical_turbulence_single_column_model_Nz108.jld2"

utles = FieldTimeSeries(lesfilename, "u")
vtles = FieldTimeSeries(lesfilename, "v")
Ttles = FieldTimeSeries(lesfilename, "T")

utscm = FieldTimeSeries(scmfilename, "u")
vtscm = FieldTimeSeries(scmfilename, "v")
Ttscm = FieldTimeSeries(scmfilename, "T")

Nt = length(utscm)
lesgrid = utles.grid
Nzles = size(lesgrid, 3)

uscm = Array(interior(utscm, 1, 1, :, 1:Nt))
ules = Array(interior(utles, 1, 1, :, 1:Nt))
Tscm = Array(interior(Ttscm, 1, 1, :, 1:Nt))
Tles = Array(interior(Ttles, 1, 1, :, 1:Nt))

t = utles.times
tdays = t ./ days
zles_u  = znodes(utles)
zscm_u  = znodes(utscm)
zles_Ri = znodes(Rtles)
zscm_Ri = znodes(Rtscm)

fig = Figure(size=(1600, 800))

ax_u1 = Axis(fig[1, 1], xlabel="Time (days)", ylabel="z (m)")
ax_u2 = Axis(fig[2, 1], xlabel="Time (days)", ylabel="z (m)")
ax_du = Axis(fig[3, 1], xlabel="Time (days)", ylabel="z (m)")
ax_T1 = Axis(fig[1, 2], xlabel="Time (days)", ylabel="z (m)")
ax_T2 = Axis(fig[2, 2], xlabel="Time (days)", ylabel="z (m)")
ax_dT = Axis(fig[3, 2], xlabel="Time (days)", ylabel="z (m)")

Nt = length(utles)
Nzles = size(utles, 3)
δu = zeros(Nzles, Nt)
δT = zeros(Nzles, Nt)

ψ = Field{Nothing, Nothing, Center}(lesgrid)

for n = 1:Nt
    regrid!(ψ, utscm[n])
    δu[:, n] .= interior(utles[n], 1, 1, :) .- interior(ψ, 1, 1, :)

    regrid!(ψ, Ttscm[n])
    δT[:, n] .= interior(Ttles[n], 1, 1, :) .- interior(ψ, 1, 1, :)
end

colorrange = (-0.2, 0.2)
colormap = :balance
hm_u = heatmap!(ax_u1, tdays, zscm_u, permutedims(uscm); colormap, colorrange) 
hm_u = heatmap!(ax_u2, tdays, zles_u, permutedims(ules); colormap, colorrange) 

colorrange = (-0.1, 0.1)
hm_u = heatmap!(ax_du, tdays, zles_u, permutedims(δu); colormap, colorrange) 

colormap = :thermal
colorrange = (24.8, 25.0)
hm_u = heatmap!(ax_T1,  tdays, zscm_u, permutedims(Tscm); colormap, colorrange)  
hm_u = heatmap!(ax_T2,  tdays, zles_u, permutedims(Tles); colormap, colorrange)  

colormap = :balance
colorrange = (-0.05, 0.05)
hm_T = heatmap!(ax_dT, tdays, zles_u, permutedims(δT); colormap, colorrange) 

tlims = (6.5, 8.7)
xlims!(ax_u1, tlims...)
xlims!(ax_u2, tlims...)
xlims!(ax_du, tlims...)
xlims!(ax_T1, tlims...)
xlims!(ax_T2, tlims...)
xlims!(ax_dT, tlims...)
=#

display(fig)

