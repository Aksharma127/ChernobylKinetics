#!/usr/bin/env julia
"""
Quick smoke tests for core scenarios (short runs).
"""

using ChernobylKinetics

function peak_power(sol)
    return maximum(sol[1, :])
end

println("Running scenario smoke tests...")

solA = run_scenario("A: Baseline", baseline_params(), (0.0, 3.0))
solF1 = run_scenario("F1: Shutdown at 50%", scenario_f1_params(), (0.0, 3.0))
solI = run_scenario("I: Safe RBMK", scenario_i_safe_rbmk_params(), (0.0, 3.0))

peakA = peak_power(solA)
peakF1 = peak_power(solF1)
peakI = peak_power(solI)

@assert peakA > peakF1
@assert peakA > peakI
@assert peakF1 < 2.0
@assert peakI < 2.0

println("Scenario smoke tests passed.")
