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
