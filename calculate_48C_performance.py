"""Post-process existing 24-hour CSV files. Does NOT rerun thermal simulations."""
from pathlib import Path
import numpy as np
import pandas as pd

THRESHOLD = 48.0
FILES = [
    ('Uninsulated recirculation (8 min)', '24h_uninsulated_recirc.csv'),
    ('10 mm insulated recirculation (8 min)', '24h_insulated_recirc.csv'),
    ('No recirculation (shower draws only)', '24h_no_recirc.csv'),
]

def minutes_at_or_above(df, threshold):
    """Integrate time at/above threshold, interpolating crossing times linearly."""
    t = df['time_min'].to_numpy(dtype=float)
    temp = df['T_9A'].to_numpy(dtype=float)
    if len(t) < 2 or not np.isfinite(t).all() or not np.isfinite(temp).all():
        raise ValueError('Missing or non-finite simulation data')
    if not np.all(np.diff(t) > 0):
        raise ValueError('Simulation times must increase strictly')
    duration = np.diff(t)
    a, b = temp[:-1], temp[1:]
    fractions = np.zeros_like(duration)
    both_above = (a >= threshold) & (b >= threshold)
    crossing_up = (a < threshold) & (b >= threshold)
    crossing_down = (a >= threshold) & (b < threshold)
    fractions[both_above] = 1.0
    fractions[crossing_up] = (b[crossing_up] - threshold) / (b[crossing_up] - a[crossing_up])
    fractions[crossing_down] = (a[crossing_down] - threshold) / (a[crossing_down] - b[crossing_down])
    return float(np.sum(duration * fractions))

rows = []
for label, filename in FILES:
    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError(f'Missing {filename}. Save this script alongside the three CSV files.')
    df = pd.read_csv(path)
    minutes = minutes_at_or_above(df, THRESHOLD)
    total = float(df['time_min'].iloc[-1] - df['time_min'].iloc[0])
    if abs(total - 1440) > 0.2:
        raise ValueError(f'{filename} spans {total:.2f} min rather than 1440 min')
    rows.append({
        'Configuration': label,
        'Time at or above 48 C (min)': round(minutes, 2),
        'Time at or above 48 C (%)': round(100 * minutes / total, 2),
        'Time below 48 C (min)': round(total - minutes, 2),
    })

results = pd.DataFrame(rows)
results.to_csv('Table_11_48C_Thermal_Performance.csv', index=False)
print('\n48 C THERMAL PERFORMANCE (24 HOURS)\n')
print(results.to_string(index=False))
print('\nSaved: Table_11_48C_Thermal_Performance.csv')
print('Note: This is thermal performance, NOT measured Legionella risk reduction.')
