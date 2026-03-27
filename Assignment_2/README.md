# Assignment 2: 2D Steady-State Heat Conduction

This module implements a suite of numerical solvers for the 2D Laplace equation. It utilizes a decoupled architecture: a high-performance Fortran 90 core for the mathematical heavy lifting, and a Python script for automated data handling and visualization.

## Solvers Implemented
1.  **Point Gauss-Seidel (PGS):** Baseline point-by-point explicit solver.
2.  **Line Gauss-Seidel (LGS):** Line-by-line implicit solver utilizing a custom TDMA inversion.
3.  **Point Successive Over-Relaxation (PSOR):** PGS accelerated by an optimized relaxation factor ($\omega$).
4.  **Line Successive Over-Relaxation (LSOR):** LGS accelerated by an optimized relaxation factor ($\omega$).
5.  **Alternating Direction Implicit (ADI):** Two-step implicit sweeps (X then Y direction) per iteration.

## Physical Domain & Cases
* **Problem 1:** $0.3\text{m} \times 0.4\text{m}$ plate ($31 \times 41$ grid) with Dirichlet boundary conditions (Bottom=40°C, Top=10°C, Left=0°C, Right=0°C).
* **Problem 2a:** Same domain with modified Top BC to 40°C.
* **Problem 2b (Symmetry):** Solved on a quarter-domain ($16 \times 21$ grid) utilizing Neumann symmetry boundaries ($\frac{\partial T}{\partial x} = 0$, $\frac{\partial T}{\partial y} = 0$) handled via 2nd-order ghost nodes to dramatically reduce computational cost.

## Code Architecture
* `cfd_solvers.f90`: Fortran module containing the numerical schemes and TDMA algorithm.
* `main.f90`: The main execution script. Handles grid initialization, boundary application, $\omega$ optimization loops, and CSV I/O.
* `plotter.py`: Python post-processing pipeline. Dynamically reads Fortran CSV outputs, mirrors the quarter-domain symmetry back to full scale, and generates 2D contours and line plots.

## Build and Execute
1. **Compile the Fortran core:**
   ```bash
   gfortran -O3 cfd_solvers.f90 main.f90 -o heat_solver
   ```
2. **Run the computational engine:**
   ```bash
   ./heat_solver
   ```
   *(Note: This will populate the `results/` directory with CSV data).*
3. **Run the post-processing pipeline:**
   ```bash
   python plotter.py
   ```
