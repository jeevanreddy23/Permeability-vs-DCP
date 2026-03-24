import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

# Data from Table 2 (BH1)
# Time (sec) and Readings from Top of Standpipe (m)
time_s = np.array([0, 10, 20, 30, 60, 120, 180, 240, 300, 360, 420, 480, 540, 600, 900, 1200, 1500, 1800, 2100, 2400, 2700])
readings_m = np.array([0, 0.01, 0.02, 0.04, 0.09, 0.18, 0.24, 0.30, 0.37, 0.41, 0.47, 0.51, 0.55, 0.57, 0.67, 0.77, 0.82, 0.86, 0.90, 0.92, 0.93])

# DCP Data (BH1) - blows per 0.1m
dcp_depths = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])
dcp_blows = np.array([1, 1, 1, 2, 2, 3, 5, 5, 8, 5, 5, 10, 20])

# Constants for Hvorslev Case F (from report)
d_pipe = 0.050 # m (50mm)
r_pipe = d_pipe / 2
D_bore = 0.075 # m (75mm)
R_bore = D_bore / 2
L_screen = 0.5 # m
A_pipe = np.pi * r_pipe**2

# Shape factor F for Hvorslev Case F
LR = L_screen / R_bore
F = (2 * np.pi * L_screen) / np.log(LR + np.sqrt(1 + LR**2))

# Assuming total borehole depth H_total = 1.3m (matches DCP end)
H_total = 1.3 
H_t = H_total - readings_m

# Calculate instantaneous infiltration rate i = dh/dt (m/s)
dt = np.diff(time_s)
dh = np.diff(readings_m)
i_rate = dh / dt

# Calculate instantaneous k: k = (A / (F * H)) * (dh/dt)
times_mid = (time_s[:-1] + time_s[1:]) / 2
depths_mid = (readings_m[:-1] + readings_m[1:]) / 2
H_mid = H_total - depths_mid
k_inst = (A_pipe / (F * H_mid)) * i_rate
k_inst_m_day = k_inst * 3600 * 24

# Interpolate DCP at the depths of the water level
dcp_interp = np.interp(depths_mid, dcp_depths, dcp_blows)

# Create a summary DataFrame
df_results = pd.DataFrame({
    'Time (s)': times_mid,
    'Depth (m)': depths_mid,
    'DCP (blows)': dcp_interp,
    'k_inst (m/day)': k_inst_m_day
})

print(df_results.to_string(index=False))

# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# Left Plot: DCP and Permeability vs Depth
ax1.plot(dcp_blows, dcp_depths, 'r-o', label='DCP (Blows/0.1m)')
ax1.set_xlabel('DCP Blows', color='r')
ax1.set_ylabel('Depth (m)')
ax1.invert_yaxis()
ax1.set_title('BH1: DCP and Permeability Profile')
ax1.grid(True, linestyle='--', alpha=0.7)

ax1b = ax1.twiny()
ax1b.plot(k_inst_m_day, depths_mid, 'b-s', label='k (m/day)')
ax1b.set_xlabel('Permeability k (m/day)', color='b')
ax1.legend(loc='upper right')
ax1b.legend(loc='lower right')

# Right Plot: Correlation k vs DCP
ax2.scatter(dcp_interp, k_inst_m_day, color='blue', edgecolor='k', s=50, label='Measured Data')

# Fit power law: k = a * DCP^b
def power_law(x, a, b):
    return a * np.power(x, b)

try:
    popt, _ = curve_fit(power_law, dcp_interp, k_inst_m_day)
    x_range = np.linspace(min(dcp_interp), max(dcp_interp), 100)
    y_fit = power_law(x_range, *popt)
    ax2.plot(x_range, y_fit, 'r--', label=f'Fit: k = {popt[0]:.3f} * DCP^{popt[1]:.2f}')
except Exception as e:
    print("Curve fit failed:", e)

ax2.set_xlabel('DCP Blows / 0.1m')
ax2.set_ylabel('Permeability k (m/day)')
ax2.set_title('Correlation: Permeability vs DCP (BH1)')
# ax2.set_yscale('log')
# ax2.set_xscale('log')
ax2.grid(True, which="both", ls="-", alpha=0.5)
ax2.legend()

plt.tight_layout()
plt.savefig('permeability_vs_dcp_bh1.png')
print("\nPlot saved as permeability_vs_dcp_bh1.png")
