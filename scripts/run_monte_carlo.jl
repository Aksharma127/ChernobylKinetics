#!/usr/bin/env julia
using ChernobylKinetics

peak_powers, peak_temps, outcomes, stats = run_monte_carlo_uq(500)
plot_monte_carlo_results(peak_powers, stats)
