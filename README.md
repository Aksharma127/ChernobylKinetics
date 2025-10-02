# RBMK-1000 Chernobyl Accident Kinetics Simulator

[![Language](https://img.shields.io/badge/Language-Julia-9558B2.svg)](https://julialang.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains a Julia-based point-kinetics simulator developed to model the final seconds of the Chernobyl disaster on April 26, 1986. The goal of this project is to provide an educational tool for understanding the complex interplay of physics—specifically the positive void coefficient, xenon poisoning, and control rod design—that led to the catastrophic power excursion.

The simulation runs several "what-if" scenarios to explore the reactor's instability and saves the detailed results to both plots and CSV files for analysis.

## The Physics Behind the Simulation

The simulation is built on the **point-kinetics equations**, a set of coupled ordinary differential equations (ODEs) that describe the time-dependent behavior of a reactor's neutron population without considering its spatial distribution.

### 1. Point-Kinetics Equations

The core of the engine models the neutron population and the concentration of delayed neutron precursors.

**Neutron Population ($n$):**
$$
\frac{dn(t)}{dt} = \frac{\rho(t) - \beta}{\Lambda} n(t) + \sum_{i=1}^{6} \lambda_i C_i(t)
$$

**Delayed Neutron Precursors ($C_i$):**
$$
\frac{dC_i(t)}{dt} = \frac{\beta_i}{\Lambda} n(t) - \lambda_i C_i(t)
$$

Where:
* $n(t)$ is the relative neutron population, proportional to reactor power.
* $C_i(t)$ is the concentration of the i-th group of delayed neutron precursors.
* $\rho(t)$ is the total reactivity of the core.
* $\beta$ is the total fraction of delayed neutrons.
* $\Lambda$ is the prompt neutron lifetime.
* $\lambda_i$ and $\beta_i$ are the decay constant and fraction for the i-th precursor group, respectively.

### 2. Reactivity Feedback

The crucial part of the simulation is the calculation of total reactivity, $\rho(t)$, which is the sum of external control and internal feedback mechanisms.

$$
\rho_{total}(t) = \rho_{control}(t) + \rho_{feedback}(T) + \rho_{xenon}(t)
$$

* **$\rho_{control}(t)$:** This represents the reactivity from the control rods. Our model simulates the **AZ-5 (SCRAM)** event, which includes the fatal, brief positive reactivity insertion from the graphite tips followed by the main negative insertion from the boron carbide absorbers.

* **$\rho_{feedback}(T)$:** This models how the reactor's power changes as its temperature changes. It has two key components in the RBMK:
    * **Fuel Temperature (Negative Feedback):** As fuel heats up, it becomes less efficient at causing fission. This is a natural, stabilizing effect, represented by a negative coefficient, $\alpha_f$.
    * **Coolant Temperature / Void (Positive Feedback):** As the water coolant boils into steam (voids), fewer neutrons are absorbed by the coolant. In the graphite-moderated RBMK, this means more neutrons are available for fission, increasing the reaction rate. This dangerous, destabilizing effect is represented by a positive coefficient, $\alpha_c$, and is the central flaw of the RBMK design.

* **$\rho_{xenon}(t)$:** Models the effect of Xenon-135, a powerful neutron-absorbing fission byproduct. The buildup of Xenon (the "Xenon Pit") played a critical role in forcing the operators into an unsafe reactor configuration.

## Code Structure

The simulation is contained in a single script/notebook and is organized into several key functions:
* `ReactorParams`: A `mutable struct` that holds all the physical and scenario-specific parameters of the reactor.
* `reactor_kinetics!`: The main function that defines the 11-equation ODE system for the solver.
* `build_u0`: A helper function that constructs the correct initial conditions (power, precursor levels, xenon concentration, etc.) for the start of the simulation.
* `run_scenario`: A high-level function that takes a set of parameters, runs a full simulation, generates plots, and saves the output data to a CSV file.
* `run_all_examples`: The main execution block that defines and runs all the different experimental scenarios.

## Scenarios Analyzed

The script runs four primary scenarios to investigate the reactor's physics:

* **Scenario A: Baseline:** Simulates the accident with the best-estimate parameters for the reactor's state on April 26, 1986.
* **Scenario B: Delayed SCRAM:** Investigates the effect of a slower control rod insertion on the power excursion.
* **Scenario C: High Void Coef:** Models a stronger positive void coefficient to see its effect on the thermal transient.
* **Scenario D: Graphite Tip Strong:** Doubles the positive reactivity from the graphite tips to analyze its sensitivity.

## How to Run the Simulation

### Prerequisites
1.  **Install Julia:** Download and install the latest version of Julia from [julialang.org](https://julialang.org/).
2.  **Clone this repository:**
    ```bash
    git clone [https://github.com/Aksharma127/ChernobylKinetics.git](https://github.com/Aksharma127/ChernobylKinetics.git)
    cd ChernobylKinetics
    ```

### Running the Code
1.  **Start Julia:** Open your terminal and type `julia`.
2.  **Enter Package Manager:** Press `]` to enter the `pkg>` prompt.
3.  **Install Dependencies:** Run the following command to install the necessary packages.
    ```julia
    pkg> add DifferentialEquations Plots DelimitedFiles
    ```
4.  **Run the Simulation:** You can run the script from inside the Julia REPL or execute the notebook in a compatible environment like VS Code.

## Expected Output
Running the script will:
1.  Generate and display a series of plots for each scenario, showing the relative power, fuel/coolant temperatures, and isotope concentrations over time.
2.  Create four CSV files in the same directory (`scenarioA_Baseline.csv`, `scenarioB_Delayed_SCRAM.csv`, etc.), containing the raw time-series data for all 11 state variables for each scenario.

## Disclaimer
This project is an educational tool intended to demonstrate the principles of reactor kinetics and the specific flaws of the RBMK design. It is a simplified **point-kinetics model** and does not account for spatial effects, which were critical in the real accident. The parameters used are based on public reports but are simplified for this model. This simulation is not a substitute for professional, high-fidelity reactor analysis codes.

## Acknowledgments
This project was developed by **Aksharma127** in collaboration with Google's AI assistant, **Gemini**. Gemini provided guidance on the physics, helped structure the code, debug errors, and analyze the results.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.
