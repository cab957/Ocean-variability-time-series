# -*- coding: utf-8 -*-
"""
Fourier analysis of irregular monthly temperature change rate (°C/s)
Robust version for date,value CSV files
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal

# =====================
# CONFIGURACIÓN
# =====================
CSV_FILE = r"C:\Users\Asus\Desktop\TFG\Python\dTint_dt_results(9524).csv"
DATE_COL = "date"
VALUE_COL = "value"

TOP_K = 5
APPLY_DETREND = True
APPLY_WINDOW = True

# =====================
# SAFE READING
# =====================
df = pd.read_csv(CSV_FILE)

# date parsing
df[DATE_COL] = pd.to_datetime(df[DATE_COL], errors="coerce")
df[VALUE_COL] = pd.to_numeric(df[VALUE_COL], errors="coerce")

# removing rows with NaN in date or value columns
df = df.dropna(subset=[DATE_COL, VALUE_COL])

# time sorting
df = df.sort_values(DATE_COL).set_index(DATE_COL)

# =====================
# CONVERSION TO MONTHLY REGULAR DATA
# =====================
# monthly average
df_monthly = df.resample("MS").mean()

# time gap interpolation
df_monthly[VALUE_COL] = df_monthly[VALUE_COL].interpolate("time")

# removing possible residue NaNs
df_monthly = df_monthly.dropna()

# =====================
# NUMERIC SERIES
# =====================
x = df_monthly[VALUE_COL].to_numpy()
n = x.size

if n < 24:
    raise ValueError("Series is too short for spectral analysis.")

# =====================
# PREPROCESSING
# =====================
x_proc = x.copy()

# removing linear trend and mean
if APPLY_DETREND:
    x_proc = signal.detrend(x_proc, type="linear")

# Hann window
if APPLY_WINDOW:
    x_proc *= signal.windows.hann(n)

# =====================
# FFT
# =====================
fs = 1.0  # 1 muestra por mes
fft_vals = np.fft.rfft(x_proc)

freq_cpm = np.fft.rfftfreq(n, d=1/fs)   # ciclos/mes
freq_cpy = freq_cpm * 12.0              # ciclos/año

amplitude = 2 * np.abs(fft_vals) / n
power = np.abs(fft_vals)**2 / n

# removing 0 component
freq = freq_cpy[1:]
amp = amplitude[1:]

# =====================
# PEAK DETECTION
# =====================
peaks, props = signal.find_peaks(
    amp,
    height=np.percentile(amp, 75)
)

order = np.argsort(props["peak_heights"])[::-1]
peaks = peaks[order][:TOP_K]

dom_freqs = freq[peaks]
dom_amps = amp[peaks]

# =====================
# RESULTS
# =====================
print("\n=== DOMINANT FREQUENCIES ===")
for f, a in zip(dom_freqs, dom_amps):
    print(f"f = {f:6.3f} cycles/year  | period ≈ {1/f:5.2f} years")

# =====================
# GRAPHICS
# =====================
plt.style.use("seaborn-v0_8-whitegrid")
fig, ax = plt.subplots(2, 1, figsize=(11, 8), constrained_layout=True)

# Month series
ax[0].plot(df_monthly.index, x, lw=1.2)
ax[0].set_title("Ta change rate – monthly average (°C/s)")
ax[0].set_ylabel("dT/dt")

# Spectrum
ax[1].plot(freq, amp, lw=1.2)
ax[1].scatter(dom_freqs, dom_amps, color="crimson", zorder=3)
ax[1].set_xlim(0, 6)
ax[1].set_xlabel("Frequency (cycle/year)")
ax[1].set_ylabel("Amplitude")
ax[1].set_title("Fourier Spectrum")

plt.show()
