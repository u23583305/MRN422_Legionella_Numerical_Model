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
- `24_hour_comparison.py` and `24_hour_comparison_all_showers.py` – Extended simulation and comparison scripts.
- `mitigation_optimisation.py` – Evaluation of mitigation operating conditions.
- `calculate_48C_performance.py` – Thermal performance analysis.
- `Experimental_data/` – Six averaged experimental Excel datasets used for validation.
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

The 24-hour simulation and mitigation analysis scripts can be executed individually, for example:

```bash
python 24_hour_comparison.py
```

```bash
python mitigation_optimisation.py
```

## Experimental Data

The `Experimental_data` folder contains the averaged temperature measurements from the six experimental datasets. These files are used as inputs for the numerical model validation.

The experimental data files must remain within this folder, and the scripts should be executed from the repository's main directory.

## Further Information

The mathematical formulation, engineering assumptions, numerical solution procedure, validation results and mitigation analysis are documented in the accompanying MRN 422 research report.

### Computational Requirements

The computational time required depends on the simulation duration, numerical time step and number of operating conditions evaluated. The extended 24-hour simulations and mitigation optimisation may require several hours to complete, depending on the computer used.

## Numerical Simulation Results

The `Results/` directory contains the numerical outputs generated using the Python model developed for this investigation.

The following files are provided:

- `24h_no_recirc.csv` — Simulated 24-hour temperature history for the no-recirculation reference case.
- `24h_uninsulated_recirc.csv` — Simulated 24-hour temperature history with recirculation and uninsulated piping.
- `24h_insulated_recirc.csv` — Simulated 24-hour temperature history with recirculation and 10 mm pipe insulation.
- `Table_10_24h_Risk_Comparison.csv` — Comparison of thermal-risk exposure and recirculation pump operation for the three scenarios.
- `Table_11_48C_Thermal_Performance.csv` — Comparison of thermal performance relative to the 48°C threshold.

These CSV files allow the reported numerical results to be inspected without rerunning the full simulations.

### Recirculation Optimisation

The `mitigation_optimisation_csv.py` script evaluates recirculation durations from 1 to 10 minutes for insulated and uninsulated piping, using a 47 °C pump activation threshold at thermocouple 12A.

Thermal risk exposure is evaluated at thermocouple 9A over the 32–47 °C temperature range and normalised to minutes per hour. An 8-minute recirculation duration was selected to balance thermal risk reduction and pump operation.

The results are provided in `Results/Optimisation_1_to_10_min_results.csv`.

To run the optimisation:

```bash
python mitigation_optimisation_csv.py
```
