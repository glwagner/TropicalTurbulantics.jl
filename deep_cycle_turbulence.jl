using NCDatasets
using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials
using JLD2

setup_filepath = "forcing_and_ics_0N140W.jld2"
file = jldopen(setup_filepath)
U = file["U"]
V = file["V"]
T = file["T"]
S = file["S"]
thermal_expansion = file["thermal_expansion"]
haline_contraction = file["haline_contraction"]
close(file)

Nz = length(z)
Nh = 64
Lh = 64

grid = RectilinearGrid(size = (Nh, Nh, Nz)
                       x = (0, Lh),
                       y = (0, Lh),
                       z = z,
                       topology = (Periodic, Periodic, Bounded))

equation_of_state = LinearEquationOfState(; thermal_expansion, haline_contraction)

model = NonhydrostaticModel(; grid,
                            advection = WENO(),
                            timestepper = :RungeKutta3,
                            tracers = (:T, :S),
                            buoyancy = SeawaterBuoyancy(; equation_of_state),
                            closure = AnisotropicMinimumDissipation())

Uᵢ = U[:, 1]
Vᵢ = U[:, 1]
Tᵢ = U[:, 1]
Sᵢ = U[:, 1]

set!(model; u=Uᵢ, v=Vᵢ, T=Tᵢ, S=Sᵢ)

simulation = Simulation(model, Δt=1e-2, stop_iteration=10)

run!(simulation)

