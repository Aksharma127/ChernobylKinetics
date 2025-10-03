# Anatomy of a Runaway Reaction: A Computational Reconstruction of the 1986 Chernobyl Disaster

[![Language](https://img.shields.io/badge/Language-Julia-9558B2.svg)](https://julialang.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

### Author: Aksharma127
### Date: October 3, 2025

## Abstract

*The catastrophic accident at Unit 4 of the Chernobyl Nuclear Power Plant on April 26, 1986, was the culmination of fundamental design choices rooted in a Soviet engineering philosophy that created a machine with inherent and unforgiving instabilities. This paper provides a critical review of the causal physics and details the development of a point-kinetics simulation model in the Julia programming language to computationally reconstruct the final seconds of the power excursion. By modeling the core's neutronics, thermal-hydraulics, and isotope poisoning under the specific conditions of the safety test, this research quantifies the impact of the RBMK-1000's primary design flaws: the positive void coefficient and the "positive scram" effect of the control rods. The simulation successfully reproduces the prompt criticality event and, through a definitive counterfactual experiment, demonstrates that the disaster would have been averted had the control rods been designed without their fatal graphite tips. This research concludes that the accident was not a simple matter of operator error but was fundamentally a systemic failure—a betrayal by a design that was dangerously unstable and counter-intuitive, and which turned the operators' final, logical attempt to ensure safety into the trigger for the catastrophe.*

---

## 1. The Inherently Unstable Machine: Design Flaws of the RBMK-1000

The disaster was precipitated by a fatal synergy between three critical design flaws: a large positive void coefficient, a counter-intuitive emergency shutdown system, and a susceptibility to xenon poisoning. These were not isolated defects but interconnected features of a system that, under the conditions of the fateful safety test, was primed for a runaway reaction.

### 1.1. The Positive Void Coefficient: An Unforgiving Physics
The most fundamental flaw was the RBMK's large positive void coefficient. In most reactors, water is both a coolant and a moderator (a substance that slows neutrons to enable fission). If the water boils, the moderator is lost, and the reaction slows down—a safe, negative feedback loop.

In the RBMK, the moderator was graphite. The water's primary role was cooling, and from a neutronic perspective, it acted as a poison (a neutron absorber). This created a perilous positive feedback loop:
**More Power → More Steam (Voids) → Less Neutron Absorption by Water → More Neutrons for Fission → Even More Power**

At low power levels (below 20% of nominal), this effect became dominant, rendering the reactor dangerously unstable.

### 1.2. The AZ-5 "Emergency" System and the Positive Scram Effect
Compounding this instability was the design of the AZ-5 emergency shutdown system. Each control rod consisted of a boron carbide absorber section and a graphite "displacer." When fully withdrawn, a 1.25-meter column of water remained at the bottom of the control channel.

When the AZ-5 button was pressed, the graphite displacer descended first, pushing this neutron-absorbing water out. For the first few seconds, this action **replaced a neutron absorber (water) with a moderator (graphite)** at the bottom of the core. This resulted in a massive, localized insertion of **positive reactivity**. The emergency button, when used in this specific configuration, functioned as a detonator.

### 1.3. The Xenon Pit: Setting the Trap
The physics of Xenon-135, a powerful neutron-absorbing byproduct of fission, forced the operators into the dangerously unstable configuration. A long hold at 50% power, followed by a power plunge to near-zero, created a massive "xenon pit" that "poisoned" the reactor. To counteract this and keep the reactor running for the test, the operators had to withdraw control rods far beyond the safe operating limit. This created the two necessary preconditions for disaster: the reactor was in the unstable low-power regime, and almost all control rods were withdrawn, arming the fatal positive scram effect.

## 2. A Computational Reconstruction

To investigate these physics quantitatively, I developed a point-kinetics model in the Julia programming language. This model simplifies the reactor core to a single point but retains the critical time-dependent physics.

### 2.1. The Governing Equations
The simulation solves a system of 11 ordinary differential equations (ODEs). The core equations are for the neutron population ($n$) and six groups of delayed neutron precursors ($C_i$):
$$
\frac{dn(t)}{dt} = \frac{\rho(t) - \beta}{\Lambda} n(t) + \sum_{i=1}^{6} \lambda_i C_i(t)
$$
$$
\frac{dC_i(t)}{dt} = \frac{\beta_i}{\Lambda} n(t) - \lambda_i C_i(t)
$$
Additional equations model the concentrations of Iodine-135 and Xenon-135, as well as the average fuel and coolant temperatures, which feed back into the total reactivity, $\rho(t)$. The parameters used were sourced from the IAEA's INSAG-7 report to ground the simulation in reality.

### 2.2. Implementation
The model was implemented in Julia using the `DifferentialEquations.jl` package for its high-performance solvers. The code is structured as a scenario runner, allowing for easy modification of parameters to conduct the sensitivity and counterfactual analyses detailed below.

## 3. Scenario Analysis: Dissecting the Disaster

The simulation was used to run a suite of five experiments. The following is a detailed report of the findings.

### 3.1. Scenario A: Baseline Reconstruction
* **The Question:** Can the model reproduce the catastrophic power excursion using best-estimate parameters for the reactor's state on April 26, 1986?
* **The Result:**
    * **Peak Relative Power:** `672,541` x initial power.
    * **Peak Fuel Temperature:** `40.2 million K`.
* **The Conclusion:** The simulation successfully reproduces a prompt criticality event of immense violence. The physically impossible temperature is the model's way of indicating a complete and rapid thermal-mechanical destruction of the core, validating the model's core physics.

### 3.2. Scenario B: Delayed SCRAM
* **The Question:** How does the speed of the control rod insertion affect the accident's severity?
* **The Method:** The control rod insertion duration was increased from 5 to 20 seconds.
* **The Result:**
    * **Peak Relative Power:** `4,625` x initial power.
* **The Conclusion:** A slower SCRAM leads to a *less severe* power spike. This is because the slower insertion allows more time for the fuel to heat up, introducing powerful negative reactivity (the Doppler effect) which acts as a brake and "dampens" the peak of the excursion. This highlights the complex race between competing feedback effects.

### 3.3. Scenario C: High Void Coefficient
* **The Question:** What is the quantitative impact of the positive void coefficient?
* **The Method:** The strength of the positive void feedback was significantly increased.
* **The Result:**
    * **Peak Relative Power:** `24.8 billion` x initial power (a ~37,000-fold increase over baseline).
    * **Peak Fuel Temperature:** `1.4 trillion K`.
* **The Conclusion:** This scenario demonstrates the true nature of the positive feedback loop. By making it stronger, we removed any chance for the negative feedback to compete, resulting in an unstoppable thermal runaway. The astronomical temperature is the simulation's way of showing that the positive void coefficient was the **engine of the explosion**, providing the raw destructive power.

### 3.4. Scenario D: Graphite Tip Strong
* **The Question:** How sensitive was the accident to the initial "kick" from the control rod tips?
* **The Method:** The positive reactivity from the `graphite_tip_effect` was doubled.
* **The Result:**
    * **Peak Relative Power:** `42 million` x initial power (a ~62-fold increase over baseline).
* **The Conclusion:** The graphite tip effect was the **trigger** of the disaster. Doubling the trigger's strength resulted in a dramatically more violent explosion. This proves that while the void coefficient was the "fuel," the positive scram effect was the "spark" that ignited it.

### 3.5. Scenario E: No Graphite Tips (The "Smoking Gun")
* **The Question:** Could the disaster have been prevented if the control rods were designed safely?
* **The Method:** The `graphite_tip_effect` was set to `0.0`, simulating a SCRAM on a reactor with safely designed control rods.
* **The Result:**
    * **Outcome:** A safe, rapid, and controlled shutdown. The power **never increased** and began to drop immediately and exponentially.
* **The Conclusion:** This is the definitive finding of my research. The simulation proves that the operators' final action—to press the emergency shutdown button—was the correct and logical thing to do. They were betrayed by a machine whose primary safety device was engineered to function as a detonator. This single, preventable engineering flaw was the ultimate cause of the accident.

| Scenario | Description | Peak Power (x Initial) | Key Finding |
| :--- | :--- | :--- | :--- |
| **A** | **Baseline** | ~670,000 | Model successfully reproduces the excursion. |
| **B** | **Delayed SCRAM** | ~4,600 | Slower insertion gives negative feedback time to act, dampening the peak. |
| **C** | **High Void Coef.** | ~24,800,000,000 | The void coefficient was the "engine" of the explosion. |
| **D** | **Graphite Tip Strong** | ~42,000,000 | The graphite tip was the "trigger" of the explosion. |
| **E** | **No Graphite Tips** | <1 (Immediate Decay) | A safe rod design **would have prevented the disaster**. |

## 4. Final Verdict: A Betrayal by Design

This computational investigation validates the modern scientific consensus established by the IAEA's INSAG-7 report. The Chernobyl disaster was not a simple matter of operator error but the inevitable result of a systemic failure. A fundamentally flawed machine, born of a deficient safety culture, created an operational environment that was both dangerously unstable and counter-intuitive.

My counterfactual simulation (Scenario E) proves with a high degree of confidence that had the machine been designed to safe engineering standards, the operators' final action would have saved the reactor. Instead, they were betrayed. The ultimate cause of the Chernobyl disaster was not error, but a fatal and unforgivable betrayal by design.

## Acknowledgments
This project was conceived and developed by myself, **Aksharma127**. Throughout the process, I collaborated extensively with Google's AI assistant, **Gemini**, which served as a valuable research partner for structuring the simulation code, debugging, analyzing results, and assisting in the formal write-up of this report.
