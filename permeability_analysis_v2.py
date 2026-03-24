import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

# Parameters
d_pipe = 0.050 # m (50mm)
r_pipe = d_pipe / 2
D_bore = 0.075 # m (75mm)
R_bore = D_bore / 2
L_screen = 0.5 # m
A_pipe = np.pi * r_pipe**2
LR = L_screen / R_bore
F = (2 * np.pi * L_screen) / np.log(LR + np.sqrt(1 + LR**2))
H_total = 1.4 # Assumed borehole depth matching DCP range

# BH1 Data
time_s_bh1 = np.array([0, 10, 20, 30, 60, 120, 180, 240, 300, 360, 420, 480, 540, 600, 900, 1200, 1500, 1800, 2100, 2400, 2700])
readings_bh1 = np.array([0, 0.01, 0.02, 0.04, 0.09, 0.18, 0.24, 0.30, 0.37, 0.41, 0.47, 0.51, 0.55, 0.57, 0.67, 0.77, 0.82, 0.86, 0.90, 0.92, 0.93])
dcp_depths_bh1 = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])
dcp_blows_bh1 = np.array([1, 1, 1, 2, 2, 3, 5, 5, 8, 5, 5, 10, 20])

# BH2 Data
time_s_bh2 = np.array([0, 10, 20, 30, 60, 120, 180, 240, 300, 360, 420, 480, 540, 600, 900, 1200, 1500, 1800, 2100, 2400, 2700])
readings_bh2 = np.array([0, 0.01, 0.01, 0.02, 0.06, 0.07, 0.09, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.23, 0.31, 0.35, 0.42, 0.44, 0.48, 0.52, 0.55])
dcp_depths_bh2 = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4])
dcp_blows_bh2 = np.array([2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 3, 5, 5, 20])

def analyze_bh(time, readings, dcp_d, dcp_b, label):
    H_t = H_total - readings
    dt = np.diff(time)
    dh = np.diff(readings)
    t_mid = (time[:-1] + time[1:]) / 2
    r_mid = (readings[:-1] + readings[1:]) / 2
    H_mid = H_total - r_mid
    k_inst = (A_pipe / (F * H_mid)) * (dh / dt)
    k_inst_m_day = k_inst * 3600 * 24
    dcp_interp = np.interp(r_mid, dcp_d, dcp_b)
    return r_mid, k_inst_m_day, dcp_interp

r1, k1, d1 = analyze_bh(time_s_bh1, readings_bh1, dcp_depths_bh1, dcp_blows_bh1, "BH1")
r2, k2, d2 = analyze_bh(time_s_bh2, readings_bh2, dcp_depths_bh2, dcp_blows_bh2, "BH2")

# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# DCP Profiles
ax1.plot(dcp_blows_bh1, dcp_depths_bh1, 'r-o', label='DCP BH1')
ax1.plot(dcp_blows_bh2, dcp_depths_bh2, 'b-o', label='DCP BH2')
ax1.set_xlabel('DCP Blows')
ax1.set_ylabel('Depth (m)')
ax1.invert_yaxis()
ax1.set_title('DCP Profiles')
ax1.grid(True)
ax1.legend()

# Correlation
ax2.scatter(d1, k1, color='red', alpha=0.6, label='BH1 Data')
ax2.scatter(d2, k2, color='blue', alpha=0.6, label='BH2 Data')

all_dcp = np.concatenate([d1, d2])
all_k = np.concatenate([k1, k2])

def power_law(x, a, b):
    return a * np.power(x, b)

try:
    popt, _ = curve_fit(power_law, all_dcp, all_k)
    x_range = np.linspace(1, 10, 100)
    y_range = power_law(x_range, *popt)
    ax2.plot(x_range, y_range, 'k--', label=f'Combined Fit: k = {popt[0]:.3f} * DCP^{popt[1]:.2f}')
except:
    pass

ax2.set_xlabel('DCP Blows / 0.1m')
ax2.set_ylabel('Permeability k (m/day)')
ax2.set_title('Correlation: Permeability vs DCP (BH1 & BH2)')
ax2.legend()
ax2.grid(True)
ax2.set_yscale('log')

plt.tight_layout()
plt.savefig('permeability_vs_dcp_combined.png')
print("Plot saved as permeability_vs_dcp_combined.png")

# Combined Table
results = pd.DataFrame({
    'DCP_Blows': all_dcp,
    'k_m_day': all_k
})
print("\nCorrelation Summary (Combined):")
mean_k = results.groupby('DCP_Blows')['k_m_day'].mean()
print(mean_k)
