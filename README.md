# AM5630: Foundations of Computational Fluid Dynamics
**Author:** Divyam Pandey (AE25M021)

This repository contains numerical solvers and computational analyses developed for the AM5630 course. The projects demonstrate a progression from 1D transient problems using Python to 2D steady-state boundary value problems utilizing high-performance Fortran 90 computing cores.

## Repository Structure

### [Assignment 1: 1D Parabolic Solver](./Assignment_1/)
Focuses on the 1D transient heat conduction equation (parabolic). Evaluates the stability, accuracy, and error propagation of temporal discretization schemes across varying time steps ($\Delta t$). 
* **Language:** Python / Jupyter Notebook
* **Key Concepts:** Time-marching schemes, truncation error, stability analysis.

### [Assignment 2: 2D Elliptic Solvers](./Assignment_2/)
Focuses on 2D steady-state heat conduction (Laplace equation) across a rectangular domain. Features a suite of point-explicit and line-implicit iterative solvers, including a custom Tridiagonal Matrix Algorithm (TDMA) and relaxation factor optimization.
* **Core Logic:** Fortran 90
* **Post-Processing:** Python (`pandas`, `matplotlib`)
* **Key Concepts:** Point/Line Gauss-Seidel, PSOR, LSOR, ADI, Symmetry Boundaries (Neumann conditions).

## Environment Requirements
* **Fortran:** `gfortran` (or equivalent compiler for F90)
* **Python 3.x:** `numpy`, `pandas`, `matplotlib`, `jupyter`
