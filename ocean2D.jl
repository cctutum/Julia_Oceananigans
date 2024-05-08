cd(@__DIR__)
using Pkg
Pkg.activate(".")
#######################################################

using Oceananigans

#=
We are setting out to simulate a 2D (x-z grid) layer of 
fluid in Earth’s gravitational field. The bottom of the layer 
is maintained at a higher temperature than the top. This 
heating from below creates a convective motion.
=#

# The Grid
grid = RectilinearGrid(
            size = (256, 32); # number of grid points in x and z
            topology = (Periodic, Flat, Bounded), # BCs
            extent = (256, 32) # physical lengths in x and z
            )
    
# The Boundary Conditions
bc = FieldBoundaryConditions(
            top = ValueBoundaryCondition(1.0),
            bottom = ValueBoundaryCondition(20.0)
            )

# The Diffusivities
closure = ScalarDiffusivity(ν=0.05, κ=0.01)

# The Equation of State 
buoyancy = SeawaterBuoyancy(
    equation_of_state = LinearEquationOfState(thermal_expansion=0.01,
    haline_contraction=0)
    )

#=
The equation of state is a function that describes how the density 
of the fluid at any point depends on the temperature and salinity 
there (the assumption of incompressibility usually used in 
Oceananigans models means that density has no dependence on 
pressure). Our model is salt free, but our fluid will be lighter 
when it’s hotter. This is what will cause the fluid to move,
as the lighter parts will rise and the heavier parts will sink, 
driven by gravity.

“haline” is essentially a synonym for saline used by oceanographers.
=#

# The Model and Initial Conditions
model = NonhydrostaticModel(;
            grid, buoyancy, closure,
            boundary_conditions=(T=bc,), tracers=(:T, :S)
            )

# The Initial Conditions 

# We will add a small, random perturbation to the temperature field
tper(x, z) = 0.1 * rand()
set!(model; T = tper)

# The Simulation
simulation = Simulation(model; Δt=0.01, stop_time=1800)

simulation.output_writers[:velocities] = 
            JLD2OutputWriter(
                model, model.velocities,
                filename = "conv4.jld2", 
                schedule = TimeInterval(1)
                )

simulation.output_writers[:tracers] =
            JLD2OutputWriter(
                model, model.tracers,
                filename = "conv4T.jld2", 
                schedule = TimeInterval(1)
                )

run!(simulation)

# The Results

using Plots

uF = model.velocities.u;

TF = model.tracers.T;

heatmap(interior(TF, 1:grid.Nx, 1, 1:grid.Nz)';
        aspect_ratio=1, yrange=(0, 1.5grid.Nz))

heatmap(interior(uF, 1:grid.Nx, 1, 1:grid.Nz)';
        aspect_ratio=1, yrange=(0, 1.5grid.Nz))

# Run for additional 10 timesteps
simulation.stop_time+=10;
run!(simulation)

uF = FieldTimeSeries("conv4.jld2", "u")

sizeof(uF)

using Printf

i = 50;
h50 = heatmap(interior(uF[i], 1:grid.Nx, 1, 1:grid.Nz)';
        aspect_ratio=1, yrange=(0, 1.5grid.Nz),
        colorbar=:false, ylabel="z",
        annotations=[
            (0, uF.grid.Nz+15,
            text("Horizontal velocity at timestep $i", 12, :left)),
            (0, uF.grid.Nz+5,
            text((@sprintf "Max = %.3g" maximum(uF[i])), 8, :left)),
            (100, uF.grid.Nz+5,
            text((@sprintf "Min = %.3g" minimum(uF[i])), 8, :left))],
        grid=false, axis=false)

i = 100;
h100 = heatmap(interior(uF[i], 1:grid.Nx, 1, 1:grid.Nz)';
        aspect_ratio=1, yrange=(0, 1.5grid.Nz),
        colorbar=:false, ylabel="z",
        annotations=[
            (0, uF.grid.Nz+15,
            text("Horizontal velocity at timestep $i", 12, :left)),
            (0, uF.grid.Nz+5,
            text((@sprintf "Max = %.3g" maximum(uF[i])), 8, :left)),
            (100, uF.grid.Nz+5,
            text((@sprintf "Min = %.3g" minimum(uF[i])), 8, :left))],
        grid=false, axis=false)

i = 500;
h500 = heatmap(interior(uF[i], 1:grid.Nx, 1, 1:grid.Nz)';
        aspect_ratio=1, yrange=(0, 1.5grid.Nz),
        colorbar=:false, ylabel="z",
        annotations=[
            (0, uF.grid.Nz+15,
            text("Horizontal velocity at timestep $i", 12, :left)),
            (0, uF.grid.Nz+5,
            text((@sprintf "Max = %.3g" maximum(uF[i])), 8, :left)),
            (100, uF.grid.Nz+5,
            text((@sprintf "Min = %.3g" minimum(uF[i])), 8, :left))],
        grid=false, axis=false)

plot(h50, h100, h500; layout=(3,1))



# Creating an animation of an Oceananigans simulation
using Oceananigans, Reel, Plots

function heatmap_at_time(F, time, fmin, fmax, duration)
    ts = F.times
    time = time * ts[end]/duration
    i = indexin(minimum(abs.(ts .- time)), abs.(ts .- time))[1]
    xr = yr = zr = 1
    if F.grid.Nx > 1
        xr = 1:F.grid.Nx
    end
    if F.grid.Ny > 1
        yr = 1:F.grid.Ny
    end
    if F.grid.Nz > 1
        zr = 1:F.grid.Nz
    end
    heatmap(interior(F[i], xr, yr, zr)'; aspect_ratio=1, yrange=(0, 1.5F.grid.Nz),
            clim=(fmin, fmax))
end

uF = FieldTimeSeries("conv4.jld2", "u")
const fmin = 0.5minimum(uF)
const fmax = 0.5maximum(uF)
const duration = 30

function plotframe(t, dt)
    heatmap_at_time(uF, t, fmin, fmax, duration)
end

uMovie = roll(plotframe; fps=30, duration)
write("uMovie.mp4", uMovie)