using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using GLMakie
using JLD2

#=
lesfilename = "tropical_turbulence_Nz54_averages.jld2"
#whifilename = "tropical_turbulence_original_les_data.jld2"
scmfilename = "tropical_turbulence_single_column_model_Nz108.jld2"

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
=#

Nt = length(utscm)
lesgrid = utles.grid
Nzles = size(lesgrid, 3)

#=
Ri⁻¹les = zeros(Nzles+1, Nt)
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
=#

#=
uscm = Array(interior(utscm, 1, 1, :, 1:Nt))
ules = Array(interior(utles, 1, 1, :, 1:Nt))
Tscm = Array(interior(Ttscm, 1, 1, :, 1:Nt))
Tles = Array(interior(Ttles, 1, 1, :, 1:Nt))
escm = Array(interior(etscm, 1, 1, :, 1:Nt))
eles = Array(interior(etles, 1, 1, :, 1:Nt))

Ri⁻¹scm = 1 ./ interior(Rtscm, 1, 1, :, 1:Nt)

#Ri⁻¹les = 1 ./ interior(Rtles, 1, 1, :, 1:Nt)
#Ri⁻¹les = 0 * Ri⁻¹scm
=#

fig = Figure(size=(1600, 800))

ax_eles = Axis(fig[1, 1], xlabel="Time (days)", ylabel="z (m)")
ax_escm = Axis(fig[1, 2], xlabel="Time (days)", ylabel="z (m)")

ax_Rles = Axis(fig[2, 1], xlabel="Time (days)", ylabel="z (m)")
ax_Rscm = Axis(fig[2, 2], xlabel="Time (days)", ylabel="z (m)")

t = utles.times
tdays = t ./ days
zles_u  = znodes(utles)
zscm_u  = znodes(utscm)
zles_Ri = znodes(Rtles)
zscm_Ri = znodes(Rtscm)

colorrange = (0, 1e-4)
colormap = :solar
hm_e = heatmap!(ax_eles, tdays, zscm_u, permutedims(escm); colorrange, colormap)

colorrange = (0, 5e-4)
hm_e = heatmap!(ax_escm, tdays, zles_u, permutedims(eles); colorrange, colormap)

colorrange = (0, 5)
colormap = :viridis
levels = 0:0.5:5
hm_R = heatmap!(ax_Rles, tdays, zscm_Ri, permutedims(Ri⁻¹scm); colorrange, colormap) #, levels)
hm_R = heatmap!(ax_Rscm, tdays, zles_Ri, permutedims(Ri⁻¹les); colorrange, colormap) #, levels)

# Colorbar(fig[1:2, 7], hm_R, label="1 / Ri")

tlims = (6.5, 8.7)
xlims!(ax_eles, tlims...)
xlims!(ax_escm, tlims...)
xlims!(ax_Rles, tlims...)
xlims!(ax_Rscm, tlims...)

#=
fig = Figure(size=(1600, 800))

ax_u1 = Axis(fig[1, 1], xlabel="Time (days)", ylabel="z (m)")
ax_u2 = Axis(fig[2, 1], xlabel="Time (days)", ylabel="z (m)")
ax_du = Axis(fig[3, 1], xlabel="Time (days)", ylabel="z (m)")
ax_T1 = Axis(fig[1, 2], xlabel="Time (days)", ylabel="z (m)")
ax_T2 = Axis(fig[2, 2], xlabel="Time (days)", ylabel="z (m)")
ax_dT = Axis(fig[3, 2], xlabel="Time (days)", ylabel="z (m)")

#=
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
=#

#=
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
=#

tlims = (6.5, 8.7)
xlims!(ax_u1, tlims...)
xlims!(ax_u2, tlims...)
xlims!(ax_du, tlims...)
xlims!(ax_T1, tlims...)
xlims!(ax_T2, tlims...)
xlims!(ax_dT, tlims...)
=#

display(fig)

