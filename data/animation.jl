using Oceananigans
using GLMakie

xyfilename = "tropical_turbulence_Nh128_Nz108_xy.jld2"
yzfilename = "tropical_turbulence_Nh128_Nz108_yz.jld2"
xzfilename = "tropical_turbulence_Nh128_Nz108_xz.jld2"

# wxyt = FieldTimeSeries(xyfilename, "w")
# wyzt = FieldTimeSeries(yzfilename, "w")
#wxzt = FieldTimeSeries(xzfilename, "w")

x, y, z = nodes(wxzt)
grid = wxzt.grid
Lx = grid.Lx
Ly = grid.Ly
Lz = grid.Lz

n = Observable(100)
wxzn = @lift interior(wxzt[$n], :, 1, :)

fig = Figure(size=(800, 400))
wlim = 0.02
aspect = Lx/Lz
colorrange = (-wlim, wlim)
ax = Axis(fig[1, 1]; aspect, xlabel="x (m)", ylabel="z (m)")
heatmap!(ax, x, z, wxzn; colormap=:balance, colorrange)
#heatmap!(ax, wxzn; colormap=:balance, colorrange)

display(fig)

Nt = length(wxzt)
record(fig, "tropical_turbulence.mp4", 1:Nt) do nn
    n[] = nn
end

