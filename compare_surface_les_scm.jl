using NCDatasets
using Dates

whittfilepath = joinpath("data_les", "ROMS_PSH_6HRLIN_0N140W_360x360x216_22OCT2020.nc")

ds = Dataset(whittfilepath)
T0whitt = ds["tempme"][215, :]
S0whitt = ds["saltme"][215, :]
u0whitt = ds["ume"][215, :]
twhitt = ds["time"][:]
t0whitt = twhitt[1]
twhitt_seconds = 1e-3 .* map(t -> Millisecond(t - t0whitt).value, twhitt)

using Oceananigans
using Oceananigans.Units
using GLMakie

#lesfilename = "tropical_turbulence_Nz108_averages.jld2"
lesfilename = "tropical_turbulence_Nz54_averages.jld2"
scmfilename = "tropical_turbulence_single_column_model_Nz108.jld2"
#scmfilename = "tropical_turbulence_single_column_model_Nz54.jld2"
#scmfilename = "tropical_turbulence_single_column_model_Nz27.jld2"
#scmfilename = "tropical_turbulence_single_column_model_Nz216.jld2"

utles = FieldTimeSeries(lesfilename, "u")
vtles = FieldTimeSeries(lesfilename, "v")
Ttles = FieldTimeSeries(lesfilename, "T")
Stles = FieldTimeSeries(lesfilename, "S")

utscm = FieldTimeSeries(scmfilename, "u")
vtscm = FieldTimeSeries(scmfilename, "v")
Ttscm = FieldTimeSeries(scmfilename, "T")
Stscm = FieldTimeSeries(scmfilename, "S")

t = utles.times

Nzles = size(Ttles, 3)
Nzscm = size(Ttscm, 3)

T0les = Ttles[1, 1, Nzles, :]
S0les = Stles[1, 1, Nzles, :]
u0les = utles[1, 1, Nzles, :]
v0les = vtles[1, 1, Nzles, :]

T0scm = Ttscm[1, 1, Nzscm, :]
S0scm = Stscm[1, 1, Nzscm, :]
u0scm = utscm[1, 1, Nzscm, :]
v0scm = vtscm[1, 1, Nzscm, :]

fig = Figure(size=(1200, 600))

axT = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Surface temperature (ᵒC)")
axS = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Surface salinity (g kg⁻¹)")
axu = Axis(fig[3, 1], xlabel="Time (days)", ylabel="Surface zonal velocity (m s⁻¹)")

ylims!(axT, 24.6, 26.8)
ylims!(axS, 34.8, 35.1)
ylims!(axu, -1, 0.2)

lines!(axT, t ./ days, T0scm, label="CATKE")
lines!(axS, t ./ days, S0scm, label="CATKE")
lines!(axu, t ./ days, u0scm, label="CATKE")

lines!(axT, t ./ days, T0les, label="Oceananigans LES")
lines!(axS, t ./ days, S0les, label="Oceananigans LES")
lines!(axu, t ./ days, u0les, label="Oceananigans LES")

lines!(axT, twhitt_seconds ./ days, T0whitt)
lines!(axS, twhitt_seconds ./ days, S0whitt)
lines!(axu, twhitt_seconds ./ days, u0whitt)

axislegend(axS, position=:lt)

display(fig)

