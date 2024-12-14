using GLMakie
using DifferentialEquations
using GeometryBasics
using Colors

# Pendulum parameters
L = 1.0    # Arm length
m1 = 1.0    # Bob mass

# Initial conditions
θ_0 = pi/4
ω_0 = 0.0

# External parameters
g = 9.81
radius = 0.3
height = 1.5
t_start = 0.0
t_end = 20.0

function generate_cylinder(radius, height, θ_steps=100, z_steps=100)
    θ = range(0, 2π, length=θ_steps)  # Angular steps for the circular cross-section
    z = range(0, height, length=z_steps)  # Steps along the height of the cylinder

    # Parametrize the cylinder
    x = [radius * cos(θi+pi) for θi in θ, zi in z]
    y = [radius * sin(θi+pi) for θi in θ, zi in z]
    z = [zi-0.3 for θi in θ, zi in z]

    return x, y, z
end

function f(θ, ω, radius)
    x = radius .* cos.(θ.+pi)
    y = radius .* sin.(θ.+pi)
    z = ω./10 .+ 0.5
    return Point3f0.(x, y, z)
end

function pendulum(du, u, p, t)
    L, g = p
    θ, ω = u

    du[1] = ω
    du[2] = -(g/L)*sin(θ)
end

function animation()
    x, y, z = generate_cylinder(radius, height)

    tspan = (t_start, t_end)
    p = (L, g)
    
    # Solve the ODE for pendulum

    u0 = [θ_0, ω_0] # Define initial conditions
    prob = ODEProblem(pendulum, u0, tspan, p)
    sol_pendulum = solve(prob, Tsit5(), saveat=0.01)
    

    # Angle data
    θ_data = sol_pendulum[1, :]
    
    # Angular momentum data
    ω_data = sol_pendulum[2, :]

    # Compute pendulum bob positions for each pendulum
    x1_pendulum = L .* sin.(θ_data)
    y1_pendulum = -L .* cos.(θ_data)

    # Animation figure    
    fig = Figure(resolution=(1200, 600))

    # Pendulum plot
    ax1 = Axis(fig[1, 1], aspect=DataAspect())
    xlims!(ax1, -2.5, 2.5)
    ylims!(ax1, -2.5, 2.5)
    hidespines!(ax1)
    
    # Cylinder plot
    ax2 = Axis3(fig[1, 2], aspect=:data, perspectiveness=0.75)
    #hidedecorations!(ax2)
    surface!(ax2, x, y, z, color=:lightblue, transparency=true, alpha=0.5)

    rod1 = Observable([Point2f0(0.0, 0.0), Point2f0(x1_pendulum[1], y1_pendulum[1])])
    bob1 = Observable([Point2f0(x1_pendulum[1], y1_pendulum[1])])
    
    # cylinder trajectory with fading colors
    cylinder_lines = Observable(f(θ_data[1:1], ω_data[1:1], radius))
    cylinder_line_colors = Observable([RGBAf(255, 0, 0, 1.0)])


    # Separatrix initial conditions
    θ_sep_0 = 0.0
    ω_sep_0 = sqrt(4*g/L)

    # Solve the separatrix trajectory
    u0_sep = [θ_sep_0, ω_sep_0]
    prob_sep = ODEProblem(pendulum, u0_sep, tspan, p)
    sol_pendulum_sep = solve(prob_sep, Tsit5(), saveat=0.01)
    θ_sep_data = sol_pendulum_sep[1, :]
    ω_sep_data = sol_pendulum_sep[2, :]

    # Plot the separatrix
    lines!(ax2, f(θ_sep_data, ω_sep_data, radius); color=:black, linewidth=2, label="Separatrix")


    # Plot trajectory with fading colors
    lines!(ax2, cylinder_lines; color=cylinder_line_colors, linewidth=2, label="Pendulum")
    lines!(ax1, rod1; color=:black, linewidth=2)
    scatter!(ax1, bob1; color=RGBAf(255, 0, 0, 1.0), markersize=10)
    

    # Determine number of frames
    num_frames = length(θ_data)

    tail_length = 50  # Number of previous points to display

    # Animation update function
    function update(frame)
        # Update trajectory with fading tail
        start_idx = max(1, frame - tail_length)
        current_trajectory = f(θ_data[start_idx:frame], ω_data[start_idx:frame], radius)
        cylinder_lines[] = current_trajectory

        # Create fading color array
        num_points = length(current_trajectory)
        color_alphas = [((j - start_idx) / num_points)^2 for j in start_idx:frame]
        fading_colors = [RGBAf(255, 0, 0, 1 - alpha) for alpha in color_alphas]
        cylinder_line_colors[] = fading_colors

        # Update double pendulum rods and bobs
        rod1[] = [Point2f0(0.0, 0.0), Point2f0(x1_pendulum[frame], y1_pendulum[frame])]
        
        bob1[] = [Point2f0(x1_pendulum[frame], y1_pendulum[frame])]
    end
     # Record animation of pendulum and cylinder
    filename = "pendulums_sepratix.mp4"
    record(fig, filename, 1:num_frames; framerate=30) do frame
       update(frame)
    end
    return fig
end

GLMakie.activate!()
fig = animation()
fig