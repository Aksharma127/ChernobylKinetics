#!/usr/bin/env julia
"""
Comparison tests for key design scenarios (short runs).
"""

using ChernobylKinetics

function peak_power(sol)
    return maximum(sol[1, :])
end

println("Running scenario comparison tests...")

solA = run_scenario("A: Baseline", baseline_params(), (0.0, 3.0))
solF3 = run_scenario("F3: 15-rod margin", scenario_f3_params(), (0.0, 3.0))
solG = run_scenario("G: TMI-like", scenario_g_tmi_params(), (0.0, 3.0))

peakA = peak_power(solA)
peakF3 = peak_power(solF3)
peakG = peak_power(solG)

@assert peakA > peakF3
@assert peakA > peakG
@assert peakF3 > 0.5
@assert peakG > 0.5

println("Scenario comparison tests passed.")
