using DynamicalSystems
using DifferentialEquations
using ChaosTools
using GLMakie
using Statistics
using GeometryBasics

# Parameters
m1 = 1.0
m2 = 1.0
L1 = 1.0
g = 9.81

# Initial conditions
θ1_0 = pi / 2
θ2_0 = pi / 2
ω1_0 = 0.0
ω2_0 = 0.0
u0 = [θ1_0, θ2_0, ω1_0, ω2_0]

# Time span for the simulation
tmax = 100.0
dt = 0.01
tspan = (0.0, tmax)

# Function to define the double pendulum dynamics
function double_pendulum!(du, u, p, t)
    m1, m2, L1, L2, g = p
    θ1, θ2, ω1, ω2 = u

    Δθ = θ2 - θ1
    sinΔθ = sin(Δθ)
    cosΔθ = cos(Δθ)

    denom1 = (m1 + m2) * L1 - m2 * L1 * cosΔθ^2
    denom2 = (L2 / L1) * denom1

    du[1] = ω1
    du[2] = ω2

    du[3] = (m2 * L1 * ω1^2 * sinΔθ * cosΔθ +
             m2 * g * sin(θ2) * cosΔθ +
             m2 * L2 * ω2^2 * sinΔθ -
             (m1 + m2) * g * sin(θ1)) / denom1

    du[4] = (-m2 * L2 * ω2^2 * sinΔθ * cosΔθ +
             (m1 + m2) * (g * sin(θ1) * cosΔθ -
             L1 * ω1^2 * sinΔθ -
             g * sin(θ2))) / denom2
end

# Function to compute and plot Lyapunov exponent with standard deviation band
function compute_and_plot_lyapunov_band(ratios, filename)
    num_calculations = 100  # Number of LE calculations per ratio
    num_ratios = length(ratios)
    
    avg_LE = zeros(num_ratios)
    std_LE = zeros(num_ratios)
    
    for (i, ratio) in enumerate(ratios)
        L2 = ratio * L1
        p = (m1, m2, L1, L2, g)
        
        LE_values = zeros(num_calculations)
        
        for j in 1:num_calculations
            # Perturb initial conditions slightly
            δu0 = 0.01 * randn(length(u0))  # Small random perturbation
            u0_perturbed = u0 .+ δu0
            
            # Create the ContinuousDynamicalSystem
            ds = ContinuousDynamicalSystem(double_pendulum!, u0_perturbed, p)
            
            # Compute the largest Lyapunov exponent
            λ = lyapunov(ds, tspan[2]; Δt = dt)
            
            LE_values[j] = λ
        end
        
        avg_LE[i] = mean(LE_values)
        std_LE[i] = std(LE_values)
        
        println("Ratio L2/L1 = $(ratio): LE mean = $(avg_LE[i]), std = $(std_LE[i])")
    end
    
    # Identify the indices of the highest and lowest average LE
    idx_max = argmax(avg_LE)
    idx_min = argmin(avg_LE)
    
    # Get the corresponding ratios and round to 3 decimal places
    ratio_max = round(ratios[idx_max], digits=3)
    ratio_min = round(ratios[idx_min], digits=3)
    
    println("\nHighest average LE at ratio L2/L1 = $(ratio_max)")
    println("Lowest average LE at ratio L2/L1 = $(ratio_min)")
    
    # Convert ratios to a standard array for plotting
    x = collect(ratios)
    
    # Compute upper and lower bands (avg ± std)
    upper_band = avg_LE .+ std_LE
    lower_band = avg_LE .- std_LE
    
    # Plotting the results with Makie
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1],
              title = "Lyapunov Exponent vs Arm Length Ratio",
              xlabel = "Arm Length Ratio (L₂ / L₁)",
              ylabel = "Largest Lyapunov Exponent (λ₁)")
    
    # Plot average LE as a line
    lines!(ax, x, avg_LE, linewidth = 2, color = :blue, label="Average LE")
    
    # Plot band between avg-std and avg+std
    band!(ax, x, lower_band, upper_band, color = (:blue, 0.3), label="±1 Std Dev")
    
    # Highlight the highest and lowest average LE points
    scatter!(ax, [x[idx_max]], [avg_LE[idx_max]], color = :red, markersize = 10, label="Highest Avg LE")
    scatter!(ax, [x[idx_min]], [avg_LE[idx_min]], color = :green, markersize = 10, label="Lowest Avg LE")
    
    # Add a legend
    axislegend(ax, position = :rt)
    
    save(filename, fig)
    println("Plot saved as: $filename")
end

# Generate and save the plot with 50 arm length ratios
ratios_50 = LinRange(0.5, 1.5, 50)
compute_and_plot_lyapunov_band(ratios_50, "lyapunov_plot.png")
