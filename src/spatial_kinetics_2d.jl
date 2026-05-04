# ==============================================================================
#     2D SPATIAL NEUTRON DIFFUSION MODEL FOR CHERNOBYL RBMK-1000
#     Part III: Extended Spatially-Resolved Kinetics
# ==============================================================================
#
# This module implements the full 2D neutron diffusion equation coupled to
# thermal-hydraulic feedback effects. Instead of treating the entire core as
# a single "point," this model resolves spatial variations of neutron flux,
# fuel temperature, coolant void fraction, and xenon concentration across
# a 2D cross-section of the RBMK core.
#
# Governing equations:
#   1/v_g * ∂φ/∂t = D∇²φ - Σ_a*φ + (ν*Σ_f/k)*φ + S_delayed  [FLUX]
#   ∂C_i/∂t = (β_i/Λ)*φ - λ_i*C_i                             [PRECURSORS]
#   ρ_cf * ∂T_f/∂t = κ_f*∇²T_f + Q_fission - Q_coolant       [FUEL TEMP]
#
# Numerical method:
#   - Finite difference discretization (2nd order central differences)
#   - Sparse matrix storage (critical for 50x100 grid = 55,000 unknowns)
#   - Implicit time stepping (CVODE_BDF from Sundials.jl)
# ==============================================================================

using DifferentialEquations, Sundials, SparseArrays, LinearAlgebra, Plots, Printf

# === CONSTANTS ===
const GRID_R = 50          # Radial grid points
const GRID_Z = 100         # Axial grid points
const N_NODES = GRID_R * GRID_Z  # Total spatial nodes

const RADIUS_CORE = 12.0   # Core radius (m)
const HEIGHT_CORE = 14.0   # Core height (m)
const DR = RADIUS_CORE / GRID_R
const DZ = HEIGHT_CORE / GRID_Z

# Physics constants
const V_G = 1e7            # Neutron group velocity (cm/s)
const D = 1.5              # Diffusion coefficient (cm) - averaged
const SIGMA_A_BASE = 0.1   # Base absorption cross-section
const SIGMA_F = 0.095      # Fission cross-section
const NU = 2.42            # Neutrons per fission
const K_EFF = 1.0          # Effective multiplication factor (target)

const LAMBDA = 0.0005      # Prompt neutron generation time
const BETA = 0.0065        # Delayed neutron fraction
const BETA_FRACTIONS = [0.033, 0.219, 0.196, 0.395, 0.115, 0.042] .* BETA
const LAMBDA_PRECURS = [0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01]

const ALPHA_DOPPLER = -3e-3  # Temperature coefficient (dk/k/K)
const ALPHA_VOID = +0.005    # Void coefficient (dk/k per % void)

# Material properties
const RHO_FUEL = 10960.0   # UO₂ density (kg/m³)
const CP_FUEL = 300.0      # Heat capacity (J/kg/K)
const KAPPA_FUEL = 3.0     # Thermal conductivity (W/m/K)

# === GRID INITIALIZATION ===
struct SpatialGrid
    r::Vector{Float64}         # Radial positions
    z::Vector{Float64}         # Axial positions
    nodes::Int                 # Total nodes
    node_to_idx::Dict{Tuple{Int,Int}, Int}
end

function create_grid()
    r = range(0, RADIUS_CORE, length=GRID_R+1)[1:end-1] .+ DR/2
    z = range(-HEIGHT_CORE/2, HEIGHT_CORE/2, length=GRID_Z+1)[1:end-1] .+ DZ/2
    
    node_to_idx = Dict{Tuple{Int,Int}, Int}()
    for i in 1:GRID_R
        for j in 1:GRID_Z
            idx = (i-1) * GRID_Z + j
            node_to_idx[(i, j)] = idx
        end
    end
    
    return SpatialGrid(r, z, GRID_R * GRID_Z, node_to_idx)
end

# === MATERIAL PROPERTY MAP ===
@enum MaterialType FUEL COOLANT GRAPHITE CONTROL_ROD REFLECTOR

function get_material_type(r::Float64, z::Float64)::MaterialType
    # Simple material map: core is mostly fuel+coolant mix
    # z < -6m: reflector (bottom)
    # -6 to 7m: active core
    # z > 7m: reflector (top)
    
    if z < -6.0 || z > 7.0
        return REFLECTOR
    elseif r > 8.0
        return GRAPHITE  # Outer graphite moderator
    else
        return FUEL  # Simplified: treat core as fuel region
    end
end

# === SPARSE LAPLACIAN MATRIX (2D finite difference) ===
function build_laplacian_2d(grid::SpatialGrid)::SparseMatrixCSC{Float64, Int64}
    """
    Build the 2D Laplacian matrix: ∇²φ = ∂²φ/∂r² + 1/r*∂φ/∂r + ∂²φ/∂z².
    Uses 5-point stencil. Returns sparse matrix for N_NODES × N_NODES.
    """
    N = grid.nodes
    I, J, V = Int[], Int[], Float64[]
    
    # Coefficients for second-order central differences
    coeff_rr = 1.0 / (DR^2)
    coeff_r = 1.0 / (2 * DR^2)
    coeff_zz = 1.0 / (DZ^2)
    
    for i in 1:GRID_R
        for j in 1:GRID_Z
            idx = grid.node_to_idx[(i, j)]
            
            # Center coefficient (depends on stencil type)
            center_coeff = -2.0 * coeff_rr - 2.0 * coeff_zz
            
            # Radial second derivative: ∂²φ/∂r²
            if i > 1
                idx_left = grid.node_to_idx[(i-1, j)]
                push!(I, idx); push!(J, idx_left); push!(V, coeff_rr)
            end
            if i < GRID_R
                idx_right = grid.node_to_idx[(i+1, j)]
                push!(I, idx); push!(J, idx_right); push!(V, coeff_rr)
            else
                # Boundary condition: vacuum BC at r = R_core
                center_coeff += coeff_rr
            end
            
            # Radial first derivative: 1/r * ∂φ/∂r (only if r ≠ 0)
            if i > 1 && grid.r[i] > 0.1
                idx_left = grid.node_to_idx[(i-1, j)]
                idx_right = grid.node_to_idx[(i+1, j)]
                r_val = grid.r[i]
                coeff_r_val = coeff_r / r_val
                push!(I, idx); push!(J, idx_right); push!(V, coeff_r_val)
                push!(I, idx); push!(J, idx_left); push!(V, -coeff_r_val)
            end
            
            # Axial second derivative: ∂²φ/∂z²
            if j > 1
                idx_bottom = grid.node_to_idx[(i, j-1)]
                push!(I, idx); push!(J, idx_bottom); push!(V, coeff_zz)
            else
                center_coeff += coeff_zz  # Boundary condition
            end
            if j < GRID_Z
                idx_top = grid.node_to_idx[(i, j+1)]
                push!(I, idx); push!(J, idx_top); push!(V, coeff_zz)
            else
                center_coeff += coeff_zz  # Boundary condition
            end
            
            # Center coefficient
            push!(I, idx); push!(J, idx); push!(V, center_coeff)
        end
    end
    
    return sparse(I, J, V, N, N)
end

# === FULL 2D SPATIAL KINETICS SYSTEM ===
"""
State vector u has 11*N_NODES components:
  u[1:N_NODES]                 → φ(r,z) — neutron flux
  u[N_NODES+1:N_NODES*7]       → C_i(r,z) for i=1..6 — precursors
  u[N_NODES*7+1:N_NODES*8]     → I(r,z) — Iodine-135
  u[N_NODES*8+1:N_NODES*9]     → X(r,z) — Xenon-135
  u[N_NODES*9+1:N_NODES*10]    → T_f(r,z) — fuel temperature
  u[N_NODES*10+1:N_NODES*11]   → T_c(r,z) — coolant temperature
"""

mutable struct SpatialKineticsParams
    laplacian::SparseMatrixCSC{Float64, Int64}   # Pre-computed Laplacian
    grid::SpatialGrid
    void_coefficient::Float64
    graphite_tip_enabled::Bool
    scram_active::Bool
    scram_time::Float64
    t_scram_trigger::Float64
    initial_power::Float64
end

function spatial_kinetics!(du, u, p::SpatialKineticsParams, t)
    grid = p.grid
    N = grid.nodes
    lap = p.laplacian
    
    # Unpack state (views for efficiency)
    φ    = @view u[1:N]
    C    = [@view u[N + (i-1)*N + 1:N + i*N] for i in 1:6]
    I    = @view u[N*7 + 1:N*8]
    X    = @view u[N*8 + 1:N*9]
    T_f  = @view u[N*9 + 1:N*10]
    T_c  = @view u[N*10 + 1:N*11]
    
    # Unpack derivatives
    dφ   = @view du[1:N]
    dC   = [@view du[N + (i-1)*N + 1:N + i*N] for i in 1:6]
    dI   = @view du[N*7 + 1:N*8]
    dX   = @view du[N*8 + 1:N*9]
    dT_f = @view du[N*9 + 1:N*10]
    dT_c = @view du[N*10 + 1:N*11]
    
    # Compute feedback effects per node
    ρ_doppler = α_doppler_field(T_f)
    ρ_void = α_void_field(T_c, p.void_coefficient)
    ρ_xenon = xenon_reactivity_field(X)
    
    # Control rod insertion profile (SCRAM simulation)
    if !p.scram_active && t > p.t_scram_trigger
        p.scram_active = true
        p.scram_time = t
    end
    
    ρ_control = control_reactivity_field(t, p.scram_active, p.scram_time)
    
    # Total reactivity field
    ρ_total = ρ_doppler .+ ρ_void .+ ρ_xenon .+ ρ_control
    
    # === 2D NEUTRON DIFFUSION EQUATION ===
    # Point kinetics core: ((ρ-β)/Λ)*φ + sum(λ*C) = (1/v_g) * dφ/dt  [rearranged for implicit solve]
    
    # Diffusion term (apply Laplacian)
    diff_term = -D .* (lap * φ) .+ (SIGMA_A_BASE .+ ρ_total) .* φ
    
    # Fission term
    fission_term = (NU * SIGMA_F / K_EFF) .* φ
    
    # Precursor coupling
    precursor_sum = zeros(N)
    for i in 1:6
        precursor_sum .+= LAMBDA_PRECURS[i] .* C[i]
    end
    
    dφ .= ((fission_term .- diff_term .- (BETA .- ρ_total) .* φ) ./ LAMBDA) .+ precursor_sum .* φ ./ LAMBDA
    
    # === DELAYED PRECURSOR EQUATIONS ===
    for i in 1:6
        dC[i] .= (BETA_FRACTIONS[i] / LAMBDA) .* φ .- LAMBDA_PRECURS[i] .* C[i]
    end
    
    # === IODINE/XENON KINETICS ===
    dI .= 0.0639 .* (SIGMA_F .* φ) .- lambda_I .* I
    dX .= 0.0023 .* (SIGMA_F .* φ) .+ lambda_I .* I .- lambda_Xe .* X .- sigma_Xe .* φ .* X
    
    # === FUEL TEMPERATURE (Heat conduction equation) ===
    # ρ_f * cp_f * ∂T_f/∂t = κ_f * ∇²T_f + Q_fission - Q_coolant
    
    # Heat generation from fission: Q ~ ν*Σ_f*φ*E_fission (scaled)
    Q_fission = SIGMA_F .* φ  # Simplified fission heat
    
    # Heat loss to coolant (simple model)
    Q_coolant = 0.1 .* (T_f .- T_c)
    
    # Heat conduction in fuel
    heat_diffusion = KAPPA_FUEL .* (lap * T_f)
    
    dT_f .= (heat_diffusion .+ Q_fission .- Q_coolant) ./ (RHO_FUEL * CP_FUEL)
    
    # === COOLANT TEMPERATURE ===
    # Simplified: coolant equilibrates with fuel
    dT_c .= 0.5 .* (T_f .- T_c)
end

# === FEEDBACK FUNCTIONS ===
const lambda_I = 2.92e-5
const lambda_Xe = 2.09e-5
const sigma_Xe = 2.6e-18

function α_doppler_field(T_f::Vector)::Vector
    """Doppler coefficient reactivity field"""
    return ALPHA_DOPPLER .* (T_f .- 450.0)
end

function α_void_field(T_c::Vector, α_void::Float64)::Vector
    """Void coefficient reactivity field"""
    void_frac = clamp.((T_c .- 520.0) ./ 100.0, 0.0, 1.0)
    return α_void .* void_frac
end

function xenon_reactivity_field(X::Vector)::Vector
    """Xenon poisoning reactivity"""
    return -sigma_Xe .* X ./ max(SIGMA_F, eps())
end

function control_reactivity_field(t::Float64, scram_active::Bool, scram_time::Float64)::Vector
    """Control rod insertion profile"""
    if !scram_active
        return zeros(N_NODES)
    end
    
    t_rel = t - scram_time
    if t_rel < 1.0
        # Graphite tip transient
        return fill(+0.003 * (t_rel / 1.0), N_NODES)
    elseif t_rel < 6.0
        # Boron insertion (negative reactivity)
        frac = (t_rel - 1.0) / 5.0
        return fill(-0.04 * frac, N_NODES)
    else
        return fill(-0.04, N_NODES)
    end
end

# === INITIAL CONDITIONS ===
function build_u0_spatial(grid::SpatialGrid, initial_power::Float64)::Vector
    u0 = zeros(N_NODES * 11)
    
    # Initial flux: parabolic profile (thermal reactor steady-state)
    φ0 = @view u0[1:N_NODES]
    for i in 1:grid.nodes
        r_idx = div(i-1, GRID_Z) + 1
        z_idx = mod(i-1, GRID_Z) + 1
        
        r = grid.r[r_idx]
        z = grid.z[z_idx]
        
        # Parabolic profile for steady-state
        profile_r = max(0.0, 1.0 - (r / RADIUS_CORE)^2)
        profile_z = max(0.0, 1.0 - (z / (HEIGHT_CORE/2))^2)
        
        φ0[i] = initial_power * 1e13 * profile_r * profile_z
    end
    
    # Initial precursors (equilibrium)
    for k in 1:6
        C = @view u0[N_NODES + (k-1)*N_NODES + 1:N_NODES + k*N_NODES]
        C .= (BETA_FRACTIONS[k] / (LAMBDA * LAMBDA_PRECURS[k])) .* φ0
    end
    
    # Initial iodine/xenon
    I0 = @view u0[N_NODES*7 + 1:N_NODES*8]
    X0 = @view u0[N_NODES*8 + 1:N_NODES*9]
    I0 .= 0.0639 .* (SIGMA_F .* φ0) ./ lambda_I
    X0 .= (0.0023 .* SIGMA_F .* φ0 .+ lambda_I .* I0) ./ (lambda_Xe + sigma_Xe .* φ0)
    
    # Initial temperatures
    T_f0 = @view u0[N_NODES*9 + 1:N_NODES*10]
    T_c0 = @view u0[N_NODES*10 + 1:N_NODES*11]
    T_f0 .= 450.0
    T_c0 .= 520.0
    
    return u0
end

# === SOLVER ===
function solve_spatial_kinetics(initial_power::Float64, tspan::Tuple, scram_trigger_time::Float64)
    println("\n" * "="^80)
    println("2D SPATIAL NEUTRON DIFFUSION MODEL")
    println("Grid: $(GRID_R) × $(GRID_Z) = $(N_NODES) spatial nodes")
    println("State vector size: $(N_NODES * 11) = $(11 * N_NODES) coupled ODEs")
    println("="^80)
    
    grid = create_grid()
    lap = build_laplacian_2d(grid)
    
    u0 = build_u0_spatial(grid, initial_power)
    
    params = SpatialKineticsParams(
        laplacian = lap,
        grid = grid,
        void_coefficient = +0.005,
        graphite_tip_enabled = true,
        scram_active = false,
        scram_time = 0.0,
        t_scram_trigger = scram_trigger_time,
        initial_power = initial_power
    )
    
    prob = ODEProblem(spatial_kinetics!, u0, tspan, params)
    
    println("Solving system (this may take several minutes)...")
    println("Using CVODE_BDF with GMRES Krylov solver...")
    
    # Use implicit BDF solver suitable for stiff systems
    @time sol = solve(
        prob,
        CVODE_BDF(linear_solver = :GMRES),
        reltol = 1e-4,
        abstol = 1e-6,
        saveat = 0.1,
        progress = true
    )
    
    return sol, grid
end

# === VISUALIZATION ===
function visualize_spatial_results(sol, grid, time_indices)
    """Plot heatmaps of flux distribution at selected times"""
    N = grid.nodes
    
    # Time points to visualize
    times = sol.t[time_indices]
    
    fig = nothing
    for (idx, t_idx) in enumerate(time_indices)
        u = sol.u[t_idx]
        φ = u[1:N]
        
        # Reshape to 2D grid
        flux_2d = reshape(φ, GRID_Z, GRID_R)'
        
        # Create heatmap
        p = heatmap(
            grid.z,
            grid.r,
            log10.(max.(flux_2d, 1e10)),
            title = "Neutron Flux (log scale) at t = $(round(times[idx], digits=3)) s",
            xlabel = "z (axial, m)",
            ylabel = "r (radial, m)",
            color = :hot,
            clims = (10, 15)  # Log scale bounds
        )
        
        if fig === nothing
            fig = p
        else
            fig = plot(fig, p, layout = length(time_indices))
        end
    end
    
    display(fig)
    return fig
end

# === ANALYSIS: Bottom-to-Top Prompt Criticality Wave ===
function analyze_criticality_wavefront(sol, grid)
    """
    Analyze how prompt criticality propagates spatially.
    Track which regions go critical first (bottom) vs last (top).
    """
    N = grid.nodes
    n_times = length(sol.t)
    
    # Track when each spatial region reaches prompt criticality (φ > 1e15)
    criticality_time = fill(Inf, N)
    
    for t_idx in 1:n_times
        u = sol.u[t_idx]
        φ = u[1:N]
        
        for node_idx in 1:N
            if φ[node_idx] > 1e15 && isinf(criticality_time[node_idx])
                criticality_time[node_idx] = sol.t[t_idx]
            end
        end
    end
    
    # Analyze spatial pattern
    criticality_2d = reshape(criticality_time, GRID_Z, GRID_R)'
    
    # Print report
    println("\n" * "="^80)
    println("PROMPT CRITICALITY WAVE ANALYSIS")
    println("="^80)
    
    finite_times = filter(x -> isfinite(x), criticality_time)
    if !isempty(finite_times)
        min_time = minimum(finite_times)
        max_time = maximum(finite_times)
        
        println("First region to go critical: t = $(round(min_time, digits=3)) s")
        println("Last region to go critical:  t = $(round(max_time, digits=3)) s")
        println("Wavefront propagation time:  $(round(max_time - min_time, digits=3)) s")
        
        # Bottom vs top
        bottom_zones = criticality_time[1:10]  # First 10 axial zones
        top_zones = criticality_time[end-9:end]  # Last 10 axial zones
        
        bottom_critical = minimum(filter(isfinite, bottom_zones))
        top_critical = minimum(filter(isfinite, top_zones))
        
        if isfinite(bottom_critical) && isfinite(top_critical)
            println("Bottom goes critical at: t = $(round(bottom_critical, digits=3)) s")
            println("Top goes critical at:    t = $(round(top_critical, digits=3)) s")
            println("Upward propagation time: $(round(top_critical - bottom_critical, digits=3)) s")
        end
    else
        println("No Prompt Criticality Achieved (Safe Shutdown)")
    end
    println("="^80 * "\n")
    
    return criticality_2d
end

# === MAIN EXECUTION ===
if @isdefined(run_spatial_model) && run_spatial_model
    println("\nStarting 2D Spatial Kinetics Model...")
    
    # Run baseline scenario with spatial resolution
    sol, grid = solve_spatial_kinetics(
        initial_power = 0.07,  # 7% initial power
        tspan = (0.0, 10.0),   # 10 seconds
        scram_trigger_time = 1.0  # SCRAM at 1.0s
    )
    
    # Visualize results
    time_indices = [1, div(length(sol.t), 4), div(length(sol.t), 2), length(sol.t)]
    visualize_spatial_results(sol, grid, time_indices)
    
    # Analyze criticality wave
    criticality_map = analyze_criticality_wavefront(sol, grid)
end

# Export functions for external use
export solve_spatial_kinetics, visualize_spatial_results, analyze_criticality_wavefront, SpatialGrid
