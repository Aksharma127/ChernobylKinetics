#!/usr/bin/env julia
using ChernobylKinetics

run_spatial_model = true
sol, grid = solve_spatial_kinetics(0.07, (0.0, 10.0), 1.0)
visualize_spatial_results(sol, grid, [1, div(length(sol.t), 2), length(sol.t)])
