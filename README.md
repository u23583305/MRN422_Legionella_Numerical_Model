
# MRN 422 – Legionella Thermal Mitigation Numerical Model

## Project Overview

This repository contains the Python numerical model developed for the University of Pretoria MRN 422 research project, *Thermal Mitigation of Legionella Risk in Hot Water Distribution Systems*.

The model simulates transient temperature behaviour in a copper hot-water distribution network, including the main piping and a dead leg. It accounts for heat losses to the surroundings, insulation, axial conduction, recirculation and simplified cavity-flow heat transfer within the dead leg.

The numerical model was validated against averaged experimental temperature measurements for insulated and uninsulated stagnation and recirculation conditions.

## Repository Contents

- `network.py`, `pipe.py`, `dead_leg.py` – Core numerical model components.
- `controller.py` – Recirculation control logic.
- `Validation_stagnation_I.py` and `Validation_stagnation_NI.py` – Insulated and uninsulated stagnation validation.
- `Validation_recirculation_I.py` and `Validation_recirculation_NI.py` – Insulated and uninsulated recirculation validation.
- `24_hour_comparison_pipe_losses_only.py` – 24-hour simulation of temperature-controlled recirculation, pipe heat losses and pump electricity consumption.
- `mitigation_optimisation.py` – Evaluation of mitigation operating conditions.
- `mitigation_optimisation_csv.py` – Recirculation optimisation and CSV output.
- `calculate_48C_performance.py` – Thermal performance analysis.
- `Experimental_data/` – Six averaged experimental Excel datasets used for validation.
- `Results/` – Numerical simulation results and energy cost comparisons.
- `requirements.txt` – Required Python packages.

## Installation

Download the repository using **Code → Download ZIP** and extract it to a folder on your computer.

Open the extracted folder in Visual Studio Code and install the required Python libraries using:

```bash
python -m pip install -r requirements.txt
```

## Running the Validation Scripts

From the project folder, execute the required validation script.

For example, to run the uninsulated stagnation validation:

```bash
python Validation_stagnation_NI.py
```

The other validation scripts can be executed in the same way.

## Running the Extended Simulations

The 24-hour simulation and mitigation analysis scripts can be executed individually.

To run the 24-hour simulation:

```bash
python 24_hour_comparison_pipe_losses_only.py
```

To run the mitigation optimisation:

```bash
python mitigation_optimisation_csv.py
```

## Experimental Data

The `Experimental_data` folder contains the averaged temperature measurements from the six experimental datasets. These files are used as inputs for the numerical model validation.

The experimental data files must remain within this folder, and the scripts should be executed from the repository's main directory.

## Numerical Simulation Results

The `Results/` directory contains the numerical outputs generated using the Python model.

The following files are provided:

- `24h_no_recirc.csv` – Simulated 24-hour temperature history for the non-recirculating reference case.
- `24h_uninsulated_recirc.csv` – Simulated 24-hour temperature history with uninsulated recirculation.
- `24h_insulated_recirc.csv` – Simulated 24-hour temperature history with 10 mm pipe insulation.
- `Table_10_24h_Risk_Comparison.csv` – Comparison of thermal-risk exposure and recirculation pump operation.
- `Table_11_48C_Thermal_Performance.csv` – Thermal performance comparison relative to the 48°C threshold.
- `Table_11_24h_Energy_Cost_Comparison.csv` – Simulated pipe heat losses, pump electricity consumption and estimated additional annual operating costs.
- `Optimisation_1_to_10_min_results.csv` – Results of the recirculation duration optimisation.

These files allow the numerical results to be inspected without rerunning the full simulations.

### Recirculation Optimisation

The `mitigation_optimisation_csv.py` script evaluates recirculation durations from 1 to 10 minutes for insulated and uninsulated piping, using a 47°C pump activation threshold at thermocouple 12A.

Thermal risk exposure is evaluated at thermocouple 9A over the 32–47°C temperature range and normalised to minutes per hour. An 8-minute recirculation duration was selected to balance thermal risk reduction and pump operation.

The results are provided in `Results/Optimisation_1_to_10_min_results.csv`.

### 24-Hour Energy and Operating Cost Analysis

The 24-hour simulation compares the heat losses and electricity requirements of non-recirculating, uninsulated recirculating and insulated recirculating systems.

The simulation considers two five-minute shower periods per day and an 8-minute temperature-controlled recirculation cycle. Heat losses from the main pipe and dead leg are calculated using the numerical model.

The additional geyser heating energy is estimated from the increase in pipe heat losses relative to the non-recirculating system. An 85 W recirculation pump and an electricity tariff of R3.83/kWh are assumed.

The estimated additional annual operating costs are R2,436.43 for uninsulated recirculation and R794.38 for insulated recirculation. These costs are calculated using the rounded daily heat-loss values reported in Appendix F4, assuming normal geyser heating requirements remain constant across configurations.

The results are provided in `Results/Table_11_24h_Energy_Cost_Comparison.csv`.

## Further Information

The mathematical formulation, engineering assumptions, numerical solution procedure, validation results and mitigation analysis are documented in the accompanying MRN 422 research report.

### Computational Requirements

The computational time required depends on the simulation duration, numerical time step and number of operating conditions evaluated. The extended 24-hour simulations and mitigation optimisation may require several hours to complete, depending on the computer used.
