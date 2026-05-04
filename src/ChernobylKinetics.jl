module ChernobylKinetics

include("scenarios_F_through_I.jl")
include("spatial_kinetics_2d.jl")
include("validation_and_uq.jl")
include("results_summary_and_comparison.jl")

export ScenarioParams,
       baseline_params,
       scenario_f1_params,
       scenario_f2_params,
       scenario_f3_params,
       scenario_g_tmi_params,
       scenario_i_safe_rbmk_params,
       run_scenario,
       run_scenario_f2,
       plot_scenario_f_comparison,
       plot_scenario_g_comparison,
       plot_scenario_h_24h_xenon,
       plot_scenario_i_safe_rbmk,
       run_monte_carlo_uq,
       plot_monte_carlo_results

export SpatialGrid,
       solve_spatial_kinetics,
       visualize_spatial_results,
       analyze_criticality_wavefront

export baseline_validation_params,
       validate_against_skala,
       plot_validation_comparison,
       run_comprehensive_monte_carlo

export print_results_summary_table,
       create_comparison_curves_plot,
       create_void_coefficient_comparison_plot,
       create_monte_carlo_histogram,
       generate_complete_results_package

end
