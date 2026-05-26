# LISA Mission: Gravitational Wave Timing Residuals

## Overview
This repository contains numerical simulations developed to model the timing residuals of gravitational waves ($\tau_{GW}$) interacting with the Laser Interferometer Space Antenna (LISA) geometry. The project focuses on the temporal evolution of these residuals and their analytical validation.

## Key Features
- **Dynamic Geometry Modeling:** Simulates the propagation of gravitational waves through LISA's orbital configuration and specific geometry.
- **Stationary Phase Approximation (SPA):** Implements SPA to analytically validate the numerical results.
- **Cosmological Sensitivity:** Analyzes the impact and sensitivity of cosmological parameters, specifically the Hubble constant, on the timing residual signatures.
- **Multi-language Implementation:** Integrates high-performance numerical routines written in Modern Fortran with Python for data analysis, manipulation, and visual representation.

## Core Technologies
- **Python & Fortran:** Fast execution of complex numerical derivations.
- **Python:** Data manipulation (NumPy, SciPy) and data visualization (Matplotlib).
- **LaTeX:** Mathematical documentation and report generation.

## Project Structure & Results

* 📂 **`theory/`**: Contains the full mathematical derivation and the Final Degree Project document. Read the [Mathematical Framework](./theory/THEORY.md) for the analytical SPA approach.
* 📂 **`figures/`**: Includes all the graphical representations of the $\tau_{GW}$ evolution and parameter sensitivity plots.
* 📂 **`code/`**: Source scripts for the numerical simulations written in Python and Modern Fortran.
