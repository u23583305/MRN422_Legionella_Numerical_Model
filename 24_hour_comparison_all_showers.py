import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from network import network
from pipe import Pipe

# 24-hour comparison using existing network solver; shower hydraulics are approximate.
T_SOURCE = 65.0
T_AMBIENT = 20.0
FLOW_RATE = 26.7  # L/min; same as optimisation. Shower flow is an assumption.
U_CAV_RECIRC = 0.0050
TRIGGER_12A = 47.0  # IMPORTANT: original optimisation uses 12A to trigger
RUN_SECONDS = 8 * 60
DAY_SECONDS = 24 * 3600
SHOWER_STARTS = (7 * 3600, 19 * 3600)
SHOWER_SECONDS = 5 * 60
DT_STAG = 5.0      # Coarser than original 1 s; convergence-check required
DT_FLOW = 0.05     # Same as original recirculation solver
CHECK_SECONDS = 60.0  # Maximum trigger detection delay ~60 s
MAIN_NAMES = ['15B', '7B', '5B', '6B', '4B', '8A', '11A', '12A']
DEAD_NAMES = ['7A', '9A']


def initial_states():
    return ({name: T_SOURCE for name in MAIN_NAMES},
            {name: T_SOURCE for name in DEAD_NAMES})


def last_states(main_curves, dead_curves):
    return ({name: main_curves[name][-1] for name in MAIN_NAMES},
            {name: dead_curves[name][-1] for name in DEAD_NAMES})


def shower_flow(main, dead, duration, t_ins):
    """Open shower branch: T2 -> 7A -> 9A -> outlet.

    Main-line flow uses the existing network solver. The open branch uses
    plug-flow displacement plus the Pipe steady-flow heat-loss relation.
    This branch calculation is an assumption requiring validation.
    """
    times, main_curves, _, _ = network.coupled_recirculation(
        T_init_main=main, T_init_dead_leg=dead,
        inlet_time=[0.0, duration],
        inlet_temperature=[T_SOURCE, T_SOURCE],
        flow_rate_L_min=FLOW_RATE, duration=duration, dt=DT_FLOW,
        t_ins=t_ins, k_ins=0.036, dead_leg_cavity_velocity=0.0)
    times = np.asarray(times)
    diameter = 0.020
    velocity = (FLOW_RATE / 1000 / 60) / (np.pi * diameter**2 / 4)
    T_T2 = np.asarray(main_curves['T2'])
    dead_curves = {}
    for name, x in [('7A', 0.07), ('9A', 0.39)]:
        delay = x / velocity
        inlet = np.interp(np.maximum(times - delay, 0), times, T_T2)
        pipe = Pipe(name, L=x, D_i=diameter, D_o=0.022,
                    t_ins=t_ins, k_ins=0.036, T_inf=T_AMBIENT)
        reference = pipe.steady_flow_temperature(T_SOURCE, velocity)[0]
        attenuation = (reference - T_AMBIENT) / (T_SOURCE - T_AMBIENT)
        values = np.where(times < delay, dead[name],
                          T_AMBIENT + (inlet - T_AMBIENT) * attenuation)
        dead_curves[name] = values.tolist()
    return times, main_curves, dead_curves


def run_segment(main, dead, duration, t_ins, mode='stagnation'):
    if mode == 'shower':
        return shower_flow(main, dead, duration, t_ins)
    if mode == 'recirculation':
        times, main_curves, dead_curves, _ = network.coupled_recirculation(
            T_init_main=main, T_init_dead_leg=dead,
            inlet_time=[0.0, duration],
            inlet_temperature=[T_SOURCE, T_SOURCE],
            flow_rate_L_min=FLOW_RATE, duration=duration, dt=DT_FLOW,
            t_ins=t_ins, k_ins=0.036,
            dead_leg_cavity_velocity=U_CAV_RECIRC)
    else:
        times, main_curves, dead_curves = network.coupled_stagnation(
            T_init_main=main, T_init_dead_leg=dead,
            duration=duration, dt=DT_STAG,
            t_ins=t_ins, k_ins=0.036)
    return np.asarray(times), main_curves, dead_curves


def simulate_case(label, t_ins, recirculation):
    """Use identical shower schedules for all cases; shower overrides the pump.

    If a shower interrupts an eight-minute pump cycle, the unused pump time
    resumes after the shower. A new cycle is only triggered once the previous
    one has completed and 12A is at or below the trigger temperature.
    """
    network.T_source = T_SOURCE
    network.T_inf = T_AMBIENT
    main, dead = initial_states()
    time_s = 0.0
    times = [0.0]
    temperatures = [dead['9A']]
    pump_seconds = 0.0
    pump_starts = 0
    pump_remaining = 0.0
    next_progress = 3600.0

    while time_s < DAY_SECONDS - 1e-6:
        remaining = DAY_SECONDS - time_s
        active_shower = next((start for start in SHOWER_STARTS
                              if start - 1e-6 <= time_s < start + SHOWER_SECONDS - 1e-6), None)
        if active_shower is not None:
            mode = 'shower'
            duration = min(active_shower + SHOWER_SECONDS - time_s, remaining)
        else:
            if recirculation and pump_remaining <= 1e-6 and main['12A'] <= TRIGGER_12A:
                pump_remaining = float(RUN_SECONDS)
                pump_starts += 1
            if recirculation and pump_remaining > 1e-6:
                mode = 'recirculation'
                duration = min(pump_remaining, remaining)
            else:
                mode = 'stagnation'
                duration = min(CHECK_SECONDS, remaining)
            future_showers = [start for start in SHOWER_STARTS if start > time_s + 1e-6]
            if future_showers:
                duration = min(duration, min(future_showers) - time_s)

        if duration <= 0:
            raise RuntimeError(f'Non-positive segment duration at t={time_s} s')
        seg_t, seg_main, seg_dead = run_segment(main, dead, duration, t_ins, mode)
        elapsed = float(seg_t[-1])
        if elapsed <= 0 or abs(elapsed - duration) > max(DT_STAG, DT_FLOW) + 1e-3:
            raise RuntimeError(f'Unexpected segment duration: requested {duration}, got {elapsed}')
        times.extend((time_s + seg_t[1:]).tolist())
        temperatures.extend(seg_dead['9A'][1:])
        main, dead = last_states(seg_main, seg_dead)
        time_s += elapsed
        if mode == 'recirculation':
            pump_seconds += elapsed
            pump_remaining = max(0.0, pump_remaining - elapsed)
        if time_s >= next_progress:
            print(f'{label}: {time_s / 3600:.1f} / {DAY_SECONDS / 3600:.1f} h', flush=True)
            next_progress += 3600.0

    result = pd.DataFrame({'time_min': np.array(times) / 60.0,
                           'T_9A': temperatures})
    if not np.all(np.diff(result['time_min']) > 0):
        raise RuntimeError('Non-increasing output times')
    if not np.isfinite(result['T_9A']).all():
        raise RuntimeError('Non-finite temperatures')
    return result, pump_seconds / 3600.0, pump_starts


def band_minutes(df, low, high):
    # Linear interpolation within each saved interval for threshold crossings.
    t = df['time_min'].to_numpy()
    y = df['T_9A'].to_numpy()
    total = 0.0
    for i in range(len(t) - 1):
        a, b = y[i], y[i + 1]
        cuts = [0.0, 1.0]
        if b != a:
            for bound in (low, high):
                u = (bound - a) / (b - a)
                if 0.0 < u < 1.0:
                    cuts.append(u)
        cuts.sort()
        for u0, u1 in zip(cuts[:-1], cuts[1:]):
            midpoint = a + (b - a) * (u0 + u1) / 2.0
            if low <= midpoint < high:
                total += (u1 - u0) * (t[i + 1] - t[i])
    return total


def make_plot(cases):
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axhspan(15, 32, color='green', alpha=0.12)
    ax.axhspan(32, 38, color='orange', alpha=0.23)
    ax.axhspan(38, 47, color='red', alpha=0.23)
    ax.axhspan(47, 70, color='green', alpha=0.12)
    for bound in (32, 38, 47):
        ax.axhline(bound, color='gray', ls='--', lw=0.7, alpha=0.6)
    styles = [
        ('Uninsulated recirculation (8 min)', 'tab:blue', '-'),
        ('10 mm insulated recirculation (8 min)', 'tab:orange', '-'),
        ('No recirculation (shower draws only)', 'dimgray', '--')]
    for name, color, linestyle in styles:
        df = cases[name]
        ax.plot(df['time_min'] / 60.0, df['T_9A'],
                color=color, ls=linestyle, lw=1.4, label=name)
    for start in SHOWER_STARTS:
        ax.axvspan(start / 3600, (start + SHOWER_SECONDS) / 3600,
                   color='gray', alpha=0.13)
    ax.set(xlim=(0, 24), ylim=(15, 70),
           xlabel='Time (hours)', ylabel='Temperature at 9A (°C)',
           title='24-hour comparison at 9A — showers at 07:00 and 19:00 in all cases')
    ax.set_xticks(np.arange(0, 25, 2))
    ax.grid(alpha=0.15)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=1)
    fig.tight_layout()
    fig.savefig('Figure_30_24h_Comparison.png', dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    definitions = [
        ('Uninsulated recirculation (8 min)', 0.0, True,
         '24h_uninsulated_recirc.csv'),
        ('10 mm insulated recirculation (8 min)', 0.010, True,
         '24h_insulated_recirc.csv'),
        ('No recirculation (shower draws only)', 0.0, False,
         '24h_no_recirc.csv')]
    cases = {}
    summary = []
    for name, thickness, recirc, filename in definitions:
        print('Starting:', name, flush=True)
        df, pump_h, pump_count = simulate_case(name, thickness, recirc)
        df.to_csv(filename, index=False)
        cases[name] = df
        moderate = band_minutes(df, 32, 38)
        very_high = band_minutes(df, 38, 47)
        summary.append({'Configuration': name, '32-38 C (min)': moderate,
                        '38-47 C (min)': very_high,
                        '32-47 C (min)': moderate + very_high,
                        'Pump hours': pump_h, 'Pump starts': pump_count})
    table = pd.DataFrame(summary)
    baseline = table.loc[2, '32-47 C (min)']
    table['Reduction vs no recirc (%)'] = (
        100 * (baseline - table['32-47 C (min)']) / baseline
        if baseline > 0 else np.nan)
    table.to_csv('Table_10_24h_Risk_Comparison.csv', index=False)
    print('\n24-HOUR RESULTS\n', table.round(2).to_string(index=False))
    make_plot(cases)
