import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from network import network


# ==========================================================
# NON-INSULATED STAGNATION
# COUPLED MODEL VALIDATION
# ==========================================================

FILE = "Experimental_data/Non_Insulated_Stagnation_AVERAGED.xlsx"
SHEET = "Average Data"

T_INF = 21.322
DT = 1.0

MAIN_TCS = [
    "15B",
    "7B",
    "5B",
    "6B",
    "4B",
    "8A",
    "11A",
    "12A"
]

DEAD_LEG_TCS = [
    "7A",
    "9A"
]

ALL_TCS = (
    MAIN_TCS
    + DEAD_LEG_TCS
)


# ==========================================================
# LOAD EXPERIMENTAL DATA
# ==========================================================

data = pd.read_excel(
    FILE,
    sheet_name=SHEET
)

time_exp = (
    data["Time (s)"].to_numpy()
)

duration = float(
    time_exp[-1]
)


# ==========================================================
# INITIAL CONDITIONS
# ==========================================================

T_init_main = {}

for tc in MAIN_TCS:

    T_init_main[tc] = (
        data[
            f"{tc} Mean (°C)"
        ].iloc[0]
    )


T_init_dead_leg = {}

for tc in DEAD_LEG_TCS:

    T_init_dead_leg[tc] = (
        data[
            f"{tc} Mean (°C)"
        ].iloc[0]
    )


# ==========================================================
# SET EXPERIMENTAL AMBIENT
# ==========================================================

network.T_inf = T_INF


for pipe in network.pipes:

    pipe.T_inf = T_INF


for dead_legs in (
    network.dead_leg_connections.values()
):

    for dead_leg in dead_legs:

        dead_leg.T_inf = T_INF


# ==========================================================
# GET STAGNANT h_i
# ==========================================================

H_I_STAG = (
    network.pipes[0].h_i_stag
)


# ==========================================================
# PRINT CONDITIONS
# ==========================================================

print(
    "\nNON-INSULATED STAGNATION "
    "- COUPLED MODEL"
)

print(
    "\nCALIBRATION CONDITION"
)

print(
    f"Stagnant h_i = "
    f"{H_I_STAG:.1f} W/m²K"
)

print(
    "\nEXPERIMENTAL CONDITIONS"
)

print(
    f"Ambient temperature = "
    f"{T_INF:.3f} °C"
)

print(
    f"Duration = "
    f"{duration:.0f} s"
)

print(
    f"Time step = "
    f"{DT:.1f} s"
)


print(
    "\nMAIN-LINE INITIAL TEMPERATURES"
)

for tc in MAIN_TCS:

    print(
        f"{tc} = "
        f"{T_init_main[tc]:.3f} °C"
    )


print(
    "\nDEAD-LEG INITIAL TEMPERATURES"
)

for tc in DEAD_LEG_TCS:

    print(
        f"{tc} = "
        f"{T_init_dead_leg[tc]:.3f} °C"
    )


# ==========================================================
# RUN COUPLED MODEL
# ==========================================================

(
    time_model,
    main_model,
    dead_model
) = network.coupled_stagnation(

    T_init_main=T_init_main,

    T_init_dead_leg=T_init_dead_leg,

    duration=duration,

    dt=DT
)


time_model = np.array(
    time_model
)


# ==========================================================
# INTERPOLATE MODEL ONTO EXPERIMENTAL TIME GRID
# ==========================================================

model = {}


for tc in MAIN_TCS:

    model[tc] = np.interp(
        time_exp,
        time_model,
        np.array(
            main_model[tc]
        )
    )


for tc in DEAD_LEG_TCS:

    model[tc] = np.interp(
        time_exp,
        time_model,
        np.array(
            dead_model[tc]
        )
    )


# ==========================================================
# EXPERIMENTAL DATA
# ==========================================================

experimental = {}
experimental_sd = {}


for tc in ALL_TCS:

    experimental[tc] = (
        data[
            f"{tc} Mean (°C)"
        ].to_numpy()
    )

    experimental_sd[tc] = (
        data[
            f"{tc} SD (°C)"
        ].to_numpy()
    )


# ==========================================================
# ERROR METRICS
# ==========================================================

def error_metrics(
    experimental_values,
    numerical_values
):

    error = (
        numerical_values
        - experimental_values
    )

    rmse = np.sqrt(
        np.mean(
            error ** 2
        )
    )

    mae = np.mean(
        np.abs(error)
    )

    max_error = np.max(
        np.abs(error)
    )

    final_error = (
        error[-1]
    )

    mean_error = (
        np.mean(error)
    )

    return (
        rmse,
        mae,
        max_error,
        final_error,
        mean_error
    )


# ==========================================================
# CALCULATE RESULTS
# ==========================================================

results = []


for tc in ALL_TCS:

    (
        rmse,
        mae,
        max_error,
        final_error,
        mean_error
    ) = error_metrics(
        experimental[tc],
        model[tc]
    )


    results.append({

        "TC": tc,

        "Experimental Initial":
            experimental[tc][0],

        "Model Initial":
            model[tc][0],

        "Experimental Final":
            experimental[tc][-1],

        "Model Final":
            model[tc][-1],

        "RMSE":
            rmse,

        "MAE":
            mae,

        "Max Error":
            max_error,

        "Final Error":
            final_error,

        "Mean Error":
            mean_error
    })


results_df = pd.DataFrame(
    results
)


# ==========================================================
# PRINT FULL RESULTS
# ==========================================================

print(
    "\nFULL COUPLED VALIDATION RESULTS"
)

print(
    "\n"
    "TC     RMSE     MAE      MAX      "
    "FINAL     MEAN"
)

print(
    "       (°C)     (°C)     (°C)     "
    "(°C)      (°C)"
)


for _, row in (
    results_df.iterrows()
):

    print(
        f"{row['TC']:<5}  "
        f"{row['RMSE']:>6.3f}   "
        f"{row['MAE']:>6.3f}   "
        f"{row['Max Error']:>6.3f}   "
        f"{row['Final Error']:>7.3f}   "
        f"{row['Mean Error']:>7.3f}"
    )


# ==========================================================
# OVERALL PERFORMANCE
# ==========================================================

all_exp = np.concatenate(
    [
        experimental[tc]
        for tc in ALL_TCS
    ]
)

all_model = np.concatenate(
    [
        model[tc]
        for tc in ALL_TCS
    ]
)


overall_error = (
    all_model
    - all_exp
)


overall_rmse = np.sqrt(
    np.mean(
        overall_error ** 2
    )
)

overall_mae = np.mean(
    np.abs(
        overall_error
    )
)

overall_max = np.max(
    np.abs(
        overall_error
    )
)


print(
    "\nOVERALL 10-TC "
    "COUPLED PERFORMANCE"
)

print(
    f"Stagnant h_i = "
    f"{H_I_STAG:.1f} W/m²K"
)

print(
    f"Overall RMSE = "
    f"{overall_rmse:.3f} °C"
)

print(
    f"Overall MAE = "
    f"{overall_mae:.3f} °C"
)

print(
    f"Overall maximum absolute error = "
    f"{overall_max:.3f} °C"
)


# ==========================================================
# INDIVIDUAL VALIDATION PLOTS
# ==========================================================

for tc in ALL_TCS:

    plt.figure()

    plt.fill_between(
        time_exp,
        experimental[tc]
        - experimental_sd[tc],
        experimental[tc]
        + experimental_sd[tc],
        alpha=0.2,
        label="Experimental ±1 SD"
    )

    plt.plot(
        time_exp,
        experimental[tc],
        label="Experimental mean"
    )

    plt.plot(
        time_exp,
        model[tc],
        "--",
        label="Coupled numerical model"
    )

    plt.xlabel(
        "Time (s)"
    )

    plt.ylabel(
        "Temperature (°C)"
    )

    plt.title(
        f"Uninsulated Stagnation "
        f"Validation - {tc}"
    )

    plt.grid(True)

    plt.legend()

    plt.tight_layout()

    plt.show()

    plt.close()


# ==========================================================
# RMSE SUMMARY
# ==========================================================

plt.figure()

plt.bar(
    results_df["TC"],
    results_df["RMSE"]
)

plt.xlabel(
    "Thermocouple"
)

plt.ylabel(
    "RMSE (°C)"
)

plt.title(
    "Coupled Non-Insulated Stagnation "
    "- RMSE"
)

plt.grid(
    True,
    axis="y"
)

plt.tight_layout()

plt.show()

plt.close()


# ==========================================================
# MEAN ERROR SUMMARY
# ==========================================================

plt.figure()

plt.bar(
    results_df["TC"],
    results_df["Mean Error"]
)

plt.axhline(
    0,
    linestyle="--"
)

plt.xlabel(
    "Thermocouple"
)

plt.ylabel(
    "Mean Model - Experiment (°C)"
)

plt.title(
    "Coupled Non-Insulated Stagnation "
    "- Mean Error"
)

plt.grid(
    True,
    axis="y"
)

plt.tight_layout()

plt.show()

plt.close()


# ==========================================================
# FINAL ERROR SUMMARY
# ==========================================================

plt.figure()

plt.bar(
    results_df["TC"],
    results_df["Final Error"]
)

plt.axhline(
    0,
    linestyle="--"
)

plt.xlabel(
    "Thermocouple"
)

plt.ylabel(
    "Final Model - Experiment (°C)"
)

plt.title(
    "Coupled Non-Insulated Stagnation "
    "- Final Error"
)

plt.grid(
    True,
    axis="y"
)

plt.tight_layout()

plt.ylim(40, 65)
plt.show()

# ==========================================================
# MAXIMUM PERCENTAGE ERROR
# ==========================================================

print(
    "\nMAXIMUM PERCENTAGE ERROR"
)

print(
    "\nTC       MAX ERROR (%)"
)

max_percentage_errors = {}

for tc in ALL_TCS:

    percentage_error = (
        np.abs(
            model[tc]
            - experimental[tc]
        )
        / experimental[tc]
    ) * 100

    max_percentage_error = (
        np.max(
            percentage_error
        )
    )

    max_percentage_errors[tc] = (
        max_percentage_error
    )

    print(
        f"{tc:<6}"
        f"{max_percentage_error:>12.3f}"
    )


overall_max_percentage_error = max(
    max_percentage_errors.values()
)

print(
    "\nOverall maximum percentage error = "
    f"{overall_max_percentage_error:.3f} %"
)

plt.close()