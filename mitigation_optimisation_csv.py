import numpy as np
import csv
from pathlib import Path
import matplotlib.pyplot as plt

from network import network


# ==========================================================
# SETTINGS
# ==========================================================

T_SOURCE = 65.0
T_AMBIENT = 20.0

FLOW_RATE = 26.7
U_CAV = 0.0050

RISK_LOW = 32.0
RISK_HIGH = 47.0

TRIGGER_TEMPERATURES = [
    47.0
]

RUN_TIMES_MIN = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

DT_STAG = 1.0
DT_RECIRC = 0.05

# Maximum stagnation period allowed
MAX_STAGNATION = 10800.0


# Save to the same directory as this script, independent of terminal location.
OUTPUT_CSV = Path(__file__).resolve().parent / "Optimisation_1_to_10_min_results.csv"
CSV_FIELDS = [
    "condition", "T_on", "run_time_min", "T_9A_stop",
    "risk_recirc_min", "time_to_next_trigger_min",
    "cycle_length_min", "risk_cooling_min", "total_risk_min",
    "normalised_risk_min_per_hour", "pump_duty_percent"
]

# Start a new results file for this run; flush after each completed case.
with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(file, fieldnames=CSV_FIELDS)
    writer.writeheader()
print(f"Saving optimisation results to: {OUTPUT_CSV}")


# ==========================================================
# INITIAL TEMPERATURES
# ==========================================================

INITIAL_MAIN = {
    "15B": T_SOURCE,
    "7B": T_SOURCE,
    "5B": T_SOURCE,
    "6B": T_SOURCE,
    "4B": T_SOURCE,
    "8A": T_SOURCE,
    "11A": T_SOURCE,
    "12A": T_SOURCE
}

INITIAL_DEAD = {
    "7A": T_SOURCE,
    "9A": T_SOURCE
}


# ==========================================================
# RISK TIME
# ==========================================================

def calculate_risk_time(
    times,
    temperatures
):

    risk_seconds = 0.0

    for i in range(1, len(times)):

        if (
            RISK_LOW
            <= temperatures[i]
            <= RISK_HIGH
        ):

            risk_seconds += (
                times[i]
                - times[i - 1]
            )

    return risk_seconds / 60.0


# ==========================================================
# FIND TRIGGER
# ==========================================================

def find_trigger_index(
    curve,
    T_on
):

    for i, T in enumerate(curve):

        if T <= T_on:
            return i

    return None


# ==========================================================
# GET MAIN-LINE STATE
# ==========================================================

def get_main_state(
    curves,
    index
):

    names = [
        "15B",
        "7B",
        "5B",
        "6B",
        "4B",
        "8A",
        "11A",
        "12A"
    ]

    return {
        name: curves[name][index]
        for name in names
    }


# ==========================================================
# GET DEAD-LEG STATE
# ==========================================================

def get_dead_state(
    curves,
    index
):

    return {
        "7A": curves["7A"][index],
        "9A": curves["9A"][index]
    }


# ==========================================================
# RUN ONE INSULATION CONDITION
# ==========================================================

def run_condition(
    condition_name,
    t_ins
):

    print("\n")
    print("=" * 95)
    print(condition_name.upper())
    print("=" * 95)

    network.T_source = T_SOURCE
    network.T_inf = T_AMBIENT

    # ======================================================
    # INITIAL STAGNATION
    # ======================================================

    print(
        "\nRunning initial stagnation..."
    )

    (
        stag_times,
        stag_main,
        stag_dead
    ) = network.coupled_stagnation(
        T_init_main=INITIAL_MAIN,
        T_init_dead_leg=INITIAL_DEAD,
        duration=MAX_STAGNATION,
        dt=DT_STAG,
        t_ins=t_ins,
        k_ins=0.036
    )

    all_results = []

    # ======================================================
    # TEST TRIGGER TEMPERATURES
    # ======================================================

    for T_on in TRIGGER_TEMPERATURES:

        trigger_index = (
            find_trigger_index(
                stag_main["12A"],
                T_on
            )
        )

        if trigger_index is None:

            print(
                f"\nTrigger {T_on:.0f} °C "
                f"was not reached."
            )

            continue

        initial_trigger_time = (
            stag_times[trigger_index]
        )

        main_at_trigger = (
            get_main_state(
                stag_main,
                trigger_index
            )
        )

        dead_at_trigger = (
            get_dead_state(
                stag_dead,
                trigger_index
            )
        )

        print("\n")
        print("-" * 95)

        print(
            f"T_on = {T_on:.0f} °C"
        )

        print(
            f"Initial trigger time = "
            f"{initial_trigger_time / 60:.2f} min"
        )

        print(
            f"9A at trigger = "
            f"{dead_at_trigger['9A']:.2f} °C"
        )

        # ==================================================
        # TEST EACH PUMP RUN TIME
        # ==================================================

        for run_time_min in RUN_TIMES_MIN:

            run_time_s = (
                run_time_min * 60.0
            )

            # ----------------------------------------------
            # RECIRCULATION
            # ----------------------------------------------

            inlet_time = [
                0.0,
                run_time_s
            ]

            inlet_temperature = [
                T_SOURCE,
                T_SOURCE
            ]

            (
                recirc_times,
                recirc_main,
                recirc_dead,
                velocity
            ) = network.coupled_recirculation(
                T_init_main=main_at_trigger,
                T_init_dead_leg=dead_at_trigger,
                inlet_time=inlet_time,
                inlet_temperature=inlet_temperature,
                flow_rate_L_min=FLOW_RATE,
                duration=run_time_s,
                dt=DT_RECIRC,
                t_ins=t_ins,
                k_ins=0.036,
                dead_leg_cavity_velocity=U_CAV
            )

            # ----------------------------------------------
            # FINAL STATE AFTER RECIRCULATION
            # ----------------------------------------------

            main_after_recirc = (
                get_main_state(
                    recirc_main,
                    -1
                )
            )

            dead_after_recirc = (
                get_dead_state(
                    recirc_dead,
                    -1
                )
            )

            T_9A_stop = (
                dead_after_recirc["9A"]
            )

            # ----------------------------------------------
            # RISK DURING RECIRCULATION
            # ----------------------------------------------

            risk_recirc = (
                calculate_risk_time(
                    recirc_times,
                    recirc_dead["9A"]
                )
            )

            # ----------------------------------------------
            # SECOND STAGNATION
            # ----------------------------------------------

            (
                cool_times,
                cool_main,
                cool_dead
            ) = network.coupled_stagnation(
                T_init_main=main_after_recirc,
                T_init_dead_leg=dead_after_recirc,
                duration=MAX_STAGNATION,
                dt=DT_STAG,
                t_ins=t_ins,
                k_ins=0.036
            )

            next_trigger_index = (
                find_trigger_index(
                    cool_main["12A"],
                    T_on
                )
            )

            if next_trigger_index is None:

                time_to_next_trigger = None
                risk_cooling = None

            else:

                time_to_next_trigger = (
                    cool_times[
                        next_trigger_index
                    ]
                )

                risk_cooling = (
                    calculate_risk_time(
                        cool_times[
                            :next_trigger_index + 1
                        ],
                        cool_dead["9A"][
                            :next_trigger_index + 1
                        ]
                    )
                )

            # ----------------------------------------------
            # TOTAL RISK PER COMPLETE REPEATING CYCLE
            # ----------------------------------------------

            if risk_cooling is None:

                total_risk = None

            else:

                total_risk = (
                    risk_recirc
                    + risk_cooling
                )

            # ----------------------------------------------
            # CYCLE LENGTH
            # ----------------------------------------------

            if time_to_next_trigger is None:

                cycle_length_min = None

            else:

                cycle_length_min = (
                    (
                        run_time_s
                        + time_to_next_trigger
                    )
                    / 60.0
                )

            # ----------------------------------------------
            # STORE
            # ----------------------------------------------

            result = {
                "condition":
                    condition_name,

                "T_on":
                    T_on,

                "run_time_min":
                    run_time_min,

                "T_9A_stop":
                    T_9A_stop,

                "risk_recirc_min":
                    risk_recirc,

                "time_to_next_trigger_min":
                    (
                        None
                        if time_to_next_trigger
                        is None
                        else
                        time_to_next_trigger
                        / 60.0
                    ),

                "cycle_length_min":
                    cycle_length_min,

                "risk_cooling_min":
                    risk_cooling,

                "total_risk_min":
                    total_risk
            }

            if cycle_length_min is not None and cycle_length_min > 0:
                result["normalised_risk_min_per_hour"] = (
                    None if total_risk is None
                    else total_risk / cycle_length_min * 60.0
                )
                result["pump_duty_percent"] = (
                    run_time_min / cycle_length_min * 100.0
                )
            else:
                result["normalised_risk_min_per_hour"] = None
                result["pump_duty_percent"] = None

            all_results.append(
                result
            )

            # Save this case immediately so a later interruption does
            # not discard the completed results.
            with OUTPUT_CSV.open("a", newline="", encoding="utf-8") as file:
                writer = csv.DictWriter(file, fieldnames=CSV_FIELDS)
                writer.writerow(result)
                file.flush()

            # ----------------------------------------------
            # PRINT
            # ----------------------------------------------

            print(
                f"\nRun time = "
                f"{run_time_min:.1f} min"
            )

            print(
                f"9A at pump stop = "
                f"{T_9A_stop:.2f} °C"
            )

            print(
                f"9A risk during recirc = "
                f"{risk_recirc:.2f} min"
            )

            if (
                time_to_next_trigger
                is not None
            ):

                print(
                    f"Time to next trigger = "
                    f"{time_to_next_trigger / 60:.2f} min"
                )

                print(
                    f"Complete cycle length = "
                    f"{cycle_length_min:.2f} min"
                )

                print(
                    f"9A risk during cooling = "
                    f"{risk_cooling:.2f} min"
                )

                print(
                    f"Total 9A risk/cycle = "
                    f"{total_risk:.2f} min"
                )

            else:

                print(
                    "Next trigger not reached."
                )

    return all_results


# ==========================================================
# NON-INSULATED
# ==========================================================

results_NI = run_condition(
    condition_name="Non-insulated",
    t_ins=0.0
)


# ==========================================================
# INSULATED
# ==========================================================

results_I = run_condition(
    condition_name="10 mm insulated",
    t_ins=0.010
)


# ==========================================================
# COMBINE RESULTS
# ==========================================================

all_results = (
    results_NI
    + results_I
)


# ==========================================================
# FINAL SUMMARY
# ==========================================================

print("\n")
print("=" * 115)
print("RECIRCULATION RUN-TIME OPTIMISATION SUMMARY")
print("=" * 115)

print(
    "Condition | T_on | Run min | "
    "9A stop | Next trigger min | "
    "Cycle min | Risk/cycle min"
)

print("-" * 115)

for result in all_results:

    next_trigger = (
        "-"
        if result[
            "time_to_next_trigger_min"
        ] is None
        else
        f"{result['time_to_next_trigger_min']:.2f}"
    )

    cycle_length = (
        "-"
        if result[
            "cycle_length_min"
        ] is None
        else
        f"{result['cycle_length_min']:.2f}"
    )

    risk_cycle = (
        "-"
        if result[
            "total_risk_min"
        ] is None
        else
        f"{result['total_risk_min']:.2f}"
    )

    print(
        f"{result['condition']:16s} | "
        f"{result['T_on']:4.0f} | "
        f"{result['run_time_min']:7.1f} | "
        f"{result['T_9A_stop']:7.2f} | "
        f"{next_trigger:16s} | "
        f"{cycle_length:9s} | "
        f"{risk_cycle:14s}"
    )


# ==========================================================
# PLOT: RISK EXPOSURE VS RUN TIME
# ==========================================================

plt.figure(
    figsize=(9, 5)
)

for condition in [
    "Non-insulated",
    "10 mm insulated"
]:

    for T_on in TRIGGER_TEMPERATURES:

        selected = [
            result
            for result in all_results
            if (
                result["condition"]
                == condition
                and result["T_on"]
                == T_on
            )
        ]

        run_times = [
            result["run_time_min"]
            for result in selected
        ]

        risks = [
            result["total_risk_min"]
            for result in selected
        ]

        plt.plot(
            run_times,
            risks,
            marker="o",
            label=(
                f"{condition}, "
                f"T_on={T_on:.0f} °C"
            )
        )

plt.xlabel(
    "Recirculation run time (min)"
)

plt.ylabel(
    "9A thermal-risk exposure per cycle (min)"
)

plt.title(
    "Effect of recirculation run time "
    "on dead-leg thermal-risk exposure"
)

plt.grid(
    alpha=0.25
)

plt.legend()

plt.tight_layout()

plt.show()


# ==========================================================
# PLOT: CYCLE LENGTH VS RUN TIME
# ==========================================================

plt.figure(
    figsize=(9, 5)
)

for condition in [
    "Non-insulated",
    "10 mm insulated"
]:

    for T_on in TRIGGER_TEMPERATURES:

        selected = [
            result
            for result in all_results
            if (
                result["condition"]
                == condition
                and result["T_on"]
                == T_on
            )
        ]

        run_times = [
            result["run_time_min"]
            for result in selected
        ]

        cycle_lengths = [
            result["cycle_length_min"]
            for result in selected
        ]

        plt.plot(
            run_times,
            cycle_lengths,
            marker="o",
            label=(
                f"{condition}, "
                f"T_on={T_on:.0f} °C"
            )
        )

plt.xlabel(
    "Recirculation run time (min)"
)

plt.ylabel(
    "Complete cycle length (min)"
)

plt.title(
    "Effect of recirculation run time "
    "on cycle duration"
)

plt.grid(
    alpha=0.25
)

plt.legend()

plt.tight_layout()

plt.show()