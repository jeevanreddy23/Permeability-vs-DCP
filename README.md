# Permeability vs DCP Analysis – Correlating Dynamic Cone Penetrometer Blow Counts with Soil Permeability for Mining & Site Investigation

This project provides a methodology and Python scripts to correlate in-situ **Dynamic Cone Penetrometer (DCP)** blow counts with **saturated hydraulic conductivity ($k_{sat}$)** determined from infiltration tests.

## Overview

The analysis uses field infiltration test data (falling head / slug test) and calculates instantaneous permeability using the **Hvorslev (Case F)** method. This is then correlated with the soil strength/density profile measured by the DCP at corresponding depths.
Built using real field data from standpipe piezometer tests (Hvorslev Case F method) and DCP testing.
Developed a power-law correlation: k = a · DCP^b showing permeability drops from ~0.18 m/day (soft soil, DCP=1) to ~0.04 m/day (stiff soil, DCP≥8).
Automates data processing, statistical fitting (scipy), and professional visualisation with Matplotlib.
Directly applicable to groundwater monitoring, slope stability, and dewatering design in open-pit mining.
Technologies: Python, Pandas, NumPy, SciPy, Matplotlib.


Key findings of this project (from silty clay site data):
- Permeability decreases logarithmically as soil density (DCP blows) increases.
- In soft layers ($DCP = 1$), $k_{sat} \approx 0.18$ m/day.
- In stiff layers ($DCP \ge 8$), $k_{sat} \approx 0.04$ m/day (matching typical geotechnical report values).

## Project Structure

- `permeability_analysis_v2.py`: The main Python analysis script using `numpy`, `matplotlib`, and `scipy`.
- `permeability_vs_dcp_combined.png`: Final visualization of the correlation for both boreholes.
- `permeability_vs_dcp_bh1.png`: Individual analysis for BH1.

## Installation

Ensure you have Python installed with the following dependencies:
```bash
pip install numpy matplotlib scipy pandas
```

## Usage

Run the analysis script to generate the correlation plots:
```bash
python permeability_analysis_v2.py
```

## Methodology

Permeability is modeled using:
$$k = \frac{A}{F \cdot H(t)} \cdot \frac{dh}{dt}$$
where $F$ is the shape factor for a screened borehole section ($L/R = 13.3$).
The correlation is fitted using a power-law relationship: $k = a \cdot DCP^b$.

---
*Created as part of the Minor Geotech Project.*
