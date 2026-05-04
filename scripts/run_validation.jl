#!/usr/bin/env julia
using ChernobylKinetics

sol, power_mw, residuals = validate_against_skala(baseline_validation_params())
plot_validation_comparison(sol, power_mw)
