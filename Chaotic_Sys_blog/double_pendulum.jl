using GLMakie
using DifferentialEquations
using GeometryBasics
using Colors

# Pendulum parameters
L1 = 1.0    # Arm length
L2 = 1.5
m1 = 1.0    # Bob mass
m2 = 1.0

# HIGH LE: 0.52 LOW LE: 1.5 (computed in LE_plot.jl)
# Initial conditions
θ_0 = pi/2
φ_0 = pi/2
ω1_0 = 0.0
ω2_0 = 0.0

# External parameters
g = 9.81
R = 1.0
r = 0.3
t_start = 0.0
t_end = 100.0

N = 1       # Number of pendulums
perturbation = 1e-3
increments = [(-N + i) * perturbation for i in 1:N]

φ_1 = Float64[]  # Array for θ_data
ω2_1 = Float64[]  # Array for ω1_data

# Torus embedding function
function f(θ, φ, R, r)
    # The 3*pi/2 factor is for aesthetic purposes, if it werent there the poincare map ring on the torus would be on the back
    x = cos.(θ .+ 3*pi/2) .* (R .+ r .* cos.(φ))
    y = sin.(θ .+ 3*pi/2) .* (R .+ r .* cos.(φ))
    z = r .* sin.(φ)
    return Point3f0.(x, y, z)
end

# Generate torus surface
function generate_torus(R, r, θ_steps=100, φ_steps=100)
    θ = range(0, 2π, length=θ_steps)
    φ = range(0, 2π, length=φ_steps)
    x = [(R + r*cos(φi))*cos(θj+3*pi/2) for φi in θ, θj in θ]
    y = [(R + r*cos(φi))*sin(θj+3*pi/2) for φi in θ, θj in θ]
    z = [r*sin(φi) for φi in φ, θj in θ]
    return x, y, z
end

function double_pendulum(du, u, p, t)
    m1, m2, L1, L2, g = p
    θ, φ, ω1, ω2 = u

    sinΔ = sin(φ - θ)
    cosΔ = cos(φ - θ)
 
    denom1 = (m1 + m2)*L1 - m2*L1*cosΔ^2
    denom2 = (L2/L1)*denom1

    du[1] = ω1
    du[2] = ω2

    du[3] = (m2*L1*ω1^2*sinΔ*cosΔ +
             m2*g*sin(φ)*cosΔ +
             m2*L2*ω2^2*sinΔ -
             (m1 + m2)*g*sin(θ)) / denom1

    du[4] = (-m2*L2*ω2^2*sinΔ*cosΔ +
             (m1 + m2)*(g*sin(θ)*cosΔ - L1*ω1^2*sinΔ - g*sin(φ))) / denom2
end

function animation()
    x, y, z = generate_torus(R, r)

    tspan = (t_start, t_end)
    p = (m1, m2, L1, L2, g)
    
    # Solve the ODE for each pendulum
    sol_pendulum = Vector{ODESolution}(undef, N)
    for i in 1:N
        u0 = [θ_0 + increments[i], φ_0 + increments[i], ω1_0, ω2_0] # Define perturbed initial conditions
        prob = ODEProblem(double_pendulum, u0, tspan, p)
        sol_pendulum[i] = solve(prob, Tsit5(), saveat=0.01)
    end

    # Angle data
    θ_data = [mod.(sol_pendulum[i][1, :], 2π) for i in 1:N]
    φ_data = [mod.(sol_pendulum[i][2, :], 2π) for i in 1:N]
    
    # Angular momentum data
    ω1_data = [(m_1+m_2).*sol_pendulum[i][3, :] for i in 1:N]
    ω2_data = [m_2.*sol_pendulum[i][4, :] for i in 1:N]

    # Compute pendulum bob positions for each pendulum
    x1_pendulum = [L1 .* sin.(θ_data[i]) for i in 1:N]
    y1_pendulum = [-L1 .* cos.(θ_data[i]) for i in 1:N]

    x2_pendulum = [x1_pendulum[i] .+ L2 .* sin.(φ_data[i]) for i in 1:N]
    y2_pendulum = [y1_pendulum[i] .- L2 .* cos.(φ_data[i]) for i in 1:N]

    # Animation figure    
    fig = Figure(resolution=(1200, 600))

    # Double pendulum plot
    ax1 = Axis(fig[1, 1], aspect=DataAspect())
    xlims!(ax1, -2.5, 2.5)
    ylims!(ax1, -2.5, 2.5)
    hidespines!(ax1)
    
    # Torus plot
    ax2 = Axis3(fig[1, 2], aspect=:data, perspectiveness=0.75)
    hidedecorations!(ax2)
    surface!(ax2, x, y, z, color=:lightblue, transparency=true, alpha=0.5)

    

    # Observables for N pendulums
    rod1 = Vector{Observable}(undef, N)
    rod2 = Vector{Observable}(undef, N)
    bob1 = Vector{Observable}(undef, N)
    bob2 = Vector{Observable}(undef, N)

    # Observables for N torus trajectories
    torus_lines = Vector{Observable}(undef, N)
    torus_line_colors = Vector{Observable}(undef, N)

    # Observables for Poincare map points
    poincare_points = Vector{Observable}(undef, N)

    # Plot initial positions of all pendulums and torus lines
    colors = [RGB(0.2, 0.5, 1.0) * (1 - i/(N-1)) + RGB(0.1, 0.1, 0.4) * (i/(N-1)) for i in 0:N-1]

    for i in 1:N
        rod1[i] = Observable([Point2f0(0.0, 0.0), Point2f0(x1_pendulum[i][1], y1_pendulum[i][1])])
        rod2[i] = Observable([Point2f0(x1_pendulum[i][1], y1_pendulum[i][1]), Point2f0(x2_pendulum[i][1], y2_pendulum[i][1])])
        bob1[i] = Observable([Point2f0(x1_pendulum[i][1], y1_pendulum[i][1])])
        bob2[i] = Observable([Point2f0(x2_pendulum[i][1], y2_pendulum[i][1])])

        # Torus trajectory with fading colors
        torus_lines[i] = Observable(f(θ_data[i][1:1], φ_data[i][1:1], R, r))
        torus_line_colors[i] = Observable([RGBAf(colors[i].r, colors[i].g, colors[i].b, 1.0)])

        # Observable for Poincare map points
        poincare_points[i] = Observable(Point3f0[])

        # Plot trajectory with fading colors
        lines!(ax2, torus_lines[i]; color=torus_line_colors[i], linewidth=2, label="Pendulum $i")
        lines!(ax1, rod1[i]; color=:black, linewidth=2)
        lines!(ax1, rod2[i]; color=:black, linewidth=2)
        scatter!(ax1, bob1[i]; color=colors[i], markersize=10)
        scatter!(ax1, bob2[i]; color=colors[i], markersize=10)
    end

    # Poincare map points on torus
    scatter!(ax2, poincare_points[1]; color=:red, markersize=15)

    # Determine number of frames
    num_frames = length(θ_data[1])

    tail_length = 50  # Number of previous points to display

    # Animation update function
    function update(frame)
        for i in 1:N
            # Update trajectory with fading tail
            start_idx = max(1, frame - tail_length)
            current_trajectory = f(θ_data[i][start_idx:frame], φ_data[i][start_idx:frame], R, r)
            torus_lines[i][] = current_trajectory

            # Create fading color array
            num_points = length(current_trajectory)
            color_alphas = [((j - start_idx) / num_points)^2 for j in start_idx:frame]
            fading_colors = [RGBAf(colors[i].r, colors[i].g, colors[i].b, 1 - alpha) for alpha in color_alphas]
            torus_line_colors[i][] = fading_colors

            # Update double pendulum rods and bobs for pendulum i
            rod1[i][] = [Point2f0(0.0, 0.0), Point2f0(x1_pendulum[i][frame], y1_pendulum[i][frame])]
            rod2[i][] = [Point2f0(x1_pendulum[i][frame], y1_pendulum[i][frame]), Point2f0(x2_pendulum[i][frame], y2_pendulum[i][frame])]

            bob1[i][] = [Point2f0(x1_pendulum[i][frame], y1_pendulum[i][frame])]
            bob2[i][] = [Point2f0(x2_pendulum[i][frame], y2_pendulum[i][frame])]
            
            # Check for Poincare map points
            if  (frame-1)>1 && θ_data[i][frame-1] > 3*pi/2 && θ_data[i][frame] < pi
                # Map point to torus coordinates
                poincare_point = f([floor(θ_data[i][frame])], [φ_data[i][frame]], R, r)[1]

                push!(φ_1, φ_data[i][frame])
                push!(ω2_1, ω2_data[i][frame])

                # Append the point to the Poincare map points
                current_points = poincare_points[i][]
                push!(current_points, poincare_point)
                poincare_points[i][] = current_points
            end
        end
    end

    
    
     # Record animation of pendulum and torus
    filename = "pendulums_poincare_fading_L2-$(L2)_L1-$(L1)_N-$(N)_sep-$(perturbation).mp4"
    record(fig, filename, 1:num_frames; framerate=30) do frame
       update(frame)
    end

    # Plot and save poincare points 
    fig2 = Figure(resolution=(600, 600))
    ax = Axis(fig2[1, 1], title="phi vs Omega2", xlabel="Phi (rad)", ylabel="Omega2 (rad/s)")
    scatter!(ax, φ_1, ω2_1; color=:blue, markersize=5, label="Data Points")
    save("phi_vs_omega2.png", fig2)

    return fig
end

GLMakie.activate!()
fig = animation()
fig