
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from network import network
from pipe import Pipe

# ==========================================================
# SETTINGS
# ==========================================================

T_SOURCE = 65.0
T_AMBIENT = 20.0

FLOW_RATE = 26.7
U_CAV_RECIRC = 0.0050

TRIGGER_12A = 47.0
RUN_SECONDS = 8 * 60

DAY_SECONDS = 24 * 3600

SHOWER_STARTS = [7 * 3600, 19 * 3600]
SHOWER_SECONDS = 5 * 60

DT_STAG = 5.0
DT_FLOW = 0.05
CHECK_SECONDS = 60.0

MAIN_NAMES = [
    "15B", "7B", "5B", "6B",
    "4B", "8A", "11A", "12A"
]

DEAD_NAMES = ["7A", "9A"]


# ==========================================================
# INITIAL TEMPERATURES
# ==========================================================

def initial_states():

    main = {
        name: T_SOURCE
        for name in MAIN_NAMES
    }

    dead = {
        name: T_SOURCE
        for name in DEAD_NAMES
    }

    return main, dead


def last_states(main_curves, dead_curves):

    main = {
        name: main_curves[name][-1]
        for name in MAIN_NAMES
    }

    dead = {
        name: dead_curves[name][-1]
        for name in DEAD_NAMES
    }

    return main, dead


# ==========================================================
# SHOWER FLOW THROUGH OPEN DEAD LEG
# ==========================================================

def shower_flow(main, dead, duration, t_ins):

    # During showering, the dead leg becomes an open
    # branch carrying water towards the showerhead.
    #
    # Flow direction:
    # Main pipe -> T2 -> 7A -> 9A -> showerhead

    times, main_curves, _, _ = (
        network.coupled_recirculation(
            T_init_main=main,
            T_init_dead_leg=dead,
            inlet_time=[0.0, duration],
            inlet_temperature=[T_SOURCE, T_SOURCE],
            flow_rate_L_min=FLOW_RATE,
            duration=duration,
            dt=DT_FLOW,
            t_ins=t_ins,
            k_ins=0.036,
            dead_leg_cavity_velocity=0.0
        )
    )

    times = np.asarray(times)

    # Pipe dimensions
    diameter = 0.020

    # Convert flow rate from L/min to m3/s
    flow_m3_s = FLOW_RATE / 1000 / 60

    area = np.pi * diameter**2 / 4
    velocity = flow_m3_s / area

    # Main pipe temperature at the branch connection
    T_T2 = np.asarray(main_curves["T2"])

    dead_curves = {}

    # Thermocouple positions along the shower branch
    positions = {
        "7A": 0.07,
        "9A": 0.39
    }

    for name, x in positions.items():

        # Time required for water to reach this position
        delay = x / velocity

        # Temperature arriving from T-piece 2
        inlet = np.interp(
            np.maximum(times - delay, 0),
            times,
            T_T2
        )

        # Pipe heat loss during shower flow
        pipe = Pipe(
            name=name,
            L=x,
            D_i=diameter,
            D_o=0.022,
            t_ins=t_ins,
            k_ins=0.036,
            T_inf=T_AMBIENT
        )

        T_reference = (
            pipe.steady_flow_temperature(
                T_SOURCE,
                velocity
            )[0]
        )

        attenuation = (
            (T_reference - T_AMBIENT) /
            (T_SOURCE - T_AMBIENT)
        )

        # Initial stagnant water is displaced by
        # incoming hot water during the shower.
        T_branch = np.where(
            times < delay,
            dead[name],
            T_AMBIENT +
            (inlet - T_AMBIENT) * attenuation
        )

        dead_curves[name] = T_branch.tolist()

    return times, main_curves, dead_curves


# ==========================================================
# RUN ONE OPERATING SEGMENT
# ==========================================================

def run_segment(main, dead, duration, t_ins,
                mode="stagnation"):

    if mode == "shower":

        return shower_flow(
            main, dead, duration, t_ins
        )

    elif mode == "recirculation":

        (
            times,
            main_curves,
            dead_curves,
            velocity
        ) = network.coupled_recirculation(
            T_init_main=main,
            T_init_dead_leg=dead,
            inlet_time=[0.0, duration],
            inlet_temperature=[T_SOURCE, T_SOURCE],
            flow_rate_L_min=FLOW_RATE,
            duration=duration,
            dt=DT_FLOW,
            t_ins=t_ins,
            k_ins=0.036,
            dead_leg_cavity_velocity=U_CAV_RECIRC
        )

    else:

        (
            times,
            main_curves,
            dead_curves
        ) = network.coupled_stagnation(
            T_init_main=main,
            T_init_dead_leg=dead,
            duration=duration,
            dt=DT_STAG,
            t_ins=t_ins,
            k_ins=0.036
        )

    return np.asarray(times), main_curves, dead_curves


# ==========================================================
# SIMULATE ONE 24-HOUR CONDITION
# ==========================================================

def simulate_case(label, t_ins, recirculation):

    network.T_source = T_SOURCE
    network.T_inf = T_AMBIENT

    main, dead = initial_states()

    time_s = 0.0

    times = [0.0]
    temperatures = [dead["9A"]]

    pump_seconds = 0.0
    pump_starts = 0

    while time_s < DAY_SECONDS - 1e-7:

        remaining = DAY_SECONDS - time_s

        # ----------------------------------------------
        # RECIRCULATION CONDITION
        # ----------------------------------------------

        if (
            recirculation
            and main["12A"] <= TRIGGER_12A
        ):

            duration = min(
                RUN_SECONDS,
                remaining
            )

            mode = "recirculation"

            pump_starts += 1
            pump_seconds += duration

        # ----------------------------------------------
        # SHOWER DRAW CONDITION
        # ----------------------------------------------

        elif (
            not recirculation
            and any(
                start <= time_s < start + SHOWER_SECONDS
                for start in SHOWER_STARTS
            )
        ):

            end = min(
                start + SHOWER_SECONDS
                for start in SHOWER_STARTS
                if start <= time_s < start + SHOWER_SECONDS
            )

            duration = min(
                end - time_s,
                remaining
            )

            mode = "shower"

        # ----------------------------------------------
        # STAGNATION CONDITION
        # ----------------------------------------------

        else:

            duration = min(
                CHECK_SECONDS,
                remaining
            )

            mode = "stagnation"

            # Stop stagnation exactly when the next
            # shower is scheduled to begin.
            if not recirculation:

                future = [
                    start
                    for start in SHOWER_STARTS
                    if start > time_s
                ]

                if future:
                    duration = min(
                        duration,
                        min(future) - time_s
                    )

        # ----------------------------------------------
        # RUN MODEL
        # ----------------------------------------------

        (
            seg_t,
            seg_main,
            seg_dead
        ) = run_segment(
            main,
            dead,
            duration,
            t_ins,
            mode
        )

        # Store continuous temperature history
        times.extend(
            (time_s + seg_t[1:]).tolist()
        )

        temperatures.extend(
            seg_dead["9A"][1:]
        )

        # Carry final temperatures into next segment
        main, dead = last_states(
            seg_main,
            seg_dead
        )

        time_s += float(seg_t[-1])

    result = pd.DataFrame({
        "time_min": np.array(times) / 60,
        "T_9A": temperatures
    })

    return (
        result,
        pump_seconds / 3600,
        pump_starts
    )


# ==========================================================
# CALCULATE THERMAL-RISK EXPOSURE
# ==========================================================

def band_minutes(df, low, high):

    time = df["time_min"].to_numpy()
    temperature = df["T_9A"].to_numpy()

    total = 0.0

    for i in range(len(time) - 1):

        T1 = temperature[i]
        T2 = temperature[i + 1]

        cuts = [0.0, 1.0]

        # Interpolate temperature threshold crossings
        if T2 != T1:

            for boundary in [low, high]:

                fraction = (
                    (boundary - T1) /
                    (T2 - T1)
                )

                if 0 < fraction < 1:
                    cuts.append(fraction)

        cuts.sort()

        for a, b in zip(cuts[:-1], cuts[1:]):

            midpoint = (
                T1 + (T2 - T1) * (a + b) / 2
            )

            if low <= midpoint < high:

                total += (
                    (b - a) *
                    (time[i + 1] - time[i])
                )

    return total


# ==========================================================
# FIGURE 30: 24-HOUR TEMPERATURE COMPARISON
# ==========================================================

def make_plot(cases):

    fig, ax = plt.subplots(figsize=(12, 6))

    # ----------------------------------------------
    # LEGIONELLA GROWTH-POTENTIAL REGIONS
    # ----------------------------------------------

    # Green: none / very low
    ax.axhspan(
        15, 32,
        color="green",
        alpha=0.12
    )

    # Orange: moderate
    ax.axhspan(
        32, 38,
        color="orange",
        alpha=0.23
    )

    # Red: very high
    ax.axhspan(
        38, 47,
        color="red",
        alpha=0.23
    )

    # Green: lower growth potential above 47 °C
    ax.axhspan(
        47, 70,
        color="green",
        alpha=0.12
    )

    # Temperature boundaries
    for boundary in [32, 38, 47]:

        ax.axhline(
            boundary,
            color="gray",
            linestyle="--",
            linewidth=0.7,
            alpha=0.6
        )

    # ----------------------------------------------
    # PLOT TEMPERATURE CURVES
    # ----------------------------------------------

    styles = [
        (
            "Uninsulated recirculation (8 min)",
            "tab:blue",
            "-"
        ),
        (
            "10 mm insulated recirculation (8 min)",
            "tab:orange",
            "-"
        ),
        (
            "No recirculation (07:00 and 19:00 draws)",
            "dimgray",
            "--"
        )
    ]

    for name, color, linestyle in styles:

        df = cases[name]

        ax.plot(
            df["time_min"] / 60,
            df["T_9A"],
            color=color,
            linestyle=linestyle,
            linewidth=1.5,
            label=name
        )

    # ----------------------------------------------
    # INDICATE SHOWER TIMES
    # ----------------------------------------------

    for start in SHOWER_STARTS:

        ax.axvspan(
            start / 3600,
            (start + SHOWER_SECONDS) / 3600,
            color="gray",
            alpha=0.13
        )

    # ----------------------------------------------
    # AXES AND FORMATTING
    # ----------------------------------------------

    ax.set_xlim(0, 24)
    ax.set_ylim(15, 70)

    ax.set_xticks(
        np.arange(0, 25, 2)
    )

    ax.set_xlabel(
        "Time (hours)"
    )

    ax.set_ylabel(
        "Temperature at thermocouple 9A (°C)"
    )

    ax.set_title(
        "24-hour dead-leg temperature comparison"
    )

    ax.grid(alpha=0.15)

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=1
    )

    plt.tight_layout()

    plt.savefig(
        "Figure_30_24h_Comparison.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


# ==========================================================
# RUN ALL THREE CONDITIONS
# ==========================================================

if __name__ == "__main__":

    definitions = [
        (
            "Uninsulated recirculation (8 min)",
            0.0,
            True,
            "24h_uninsulated_recirc.csv"
        ),
        (
            "10 mm insulated recirculation (8 min)",
            0.010,
            True,
            "24h_insulated_recirc.csv"
        ),
        (
            "No recirculation (07:00 and 19:00 draws)",
            0.0,
            False,
            "24h_no_recirc.csv"
        )
    ]

    cases = []
    plot_cases = {}
    summary = []

    for name, thickness, recirc, filename in definitions:

        print("\nStarting:", name, flush=True)

        df, pump_h, pump_count = simulate_case(
            name,
            thickness,
            recirc
        )

        df.to_csv(filename, index=False)

        plot_cases[name] = df

        # Moderate growth-potential exposure
        moderate = band_minutes(
            df, 32, 38
        )

        # Very-high growth-potential exposure
        very_high = band_minutes(
            df, 38, 47
        )

        total_risk = moderate + very_high

        summary.append({
            "Configuration": name,
            "32-38 C (min)": moderate,
            "38-47 C (min)": very_high,
            "32-47 C (min)": total_risk,
            "Pump hours": pump_h,
            "Pump starts": pump_count
        })

    # ----------------------------------------------
    # RESULTS TABLE
    # ----------------------------------------------

    table = pd.DataFrame(summary)

    # No recirculation is the reference condition
    baseline = table.loc[
        2, "32-47 C (min)"
    ]

    if baseline > 0:

        table["Reduction vs no recirc (%)"] = (
            (
                baseline -
                table["32-47 C (min)"]
            ) / baseline
        ) * 100

    else:

        table["Reduction vs no recirc (%)"] = np.nan

    print("\n")
    print("=" * 90)
    print("24-HOUR THERMAL-RISK COMPARISON")
    print("=" * 90)

    print(
        table.round(2).to_string(index=False)
    )

    table.to_csv(
        "Table_10_24h_Risk_Comparison.csv",
        index=False
    )

    # ----------------------------------------------
    # GENERATE FIGURE 30
    # ----------------------------------------------

    make_plot(plot_cases)
