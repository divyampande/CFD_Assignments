# Assignment 1: 1D Transient Heat Conduction

This module solves the 1D parabolic heat equation, analyzing the effects of varying time steps ($\Delta t$) on the stability and accuracy of the numerical solution. 

## Physical Problem
The solver evaluates the temperature distribution along a 1D rod over time, capturing the transient thermal response before steady-state is reached. The core focus is on the numerical error introduced by the time integration scheme.

## Code Architecture
* **`1D_parabolic_solver.py`**: The core Python script containing the numerical scheme and discretization logic.
* **`AE25M021.ipynb`**: A Jupyter Notebook used for interactive execution, analysis, and visualization of the results.

## Results & Analysis
The solver outputs the temperature distributions and corresponding error plots at three specific time step resolutions:
* $\Delta t = 0.001$
* $\Delta t = 0.01$
* $\Delta t = 0.1$

These plots (saved as `.png` files) demonstrate how larger time steps introduce truncation errors or numerical instability, heavily dependent on the grid Fourier number.

## Usage
To reproduce the analysis, execute the Jupyter notebook:
```bash
jupyter notebook AE25M021.ipynb
```
Alternatively, run the standalone Python solver:
```bash
python 1D_parabolic_solver.py
```
