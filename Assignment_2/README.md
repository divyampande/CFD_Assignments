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

## Repository Structure
```text
├── results/             # Directory for storing results
├── cfd_solvers.f90      # Core Fortran module (Solvers & TDMA)
├── main.f90             # Main execution script
├── plotter.py           # Python visualization pipeline
├── sample_input.in      # Template for grid and boundary parameters
├── .gitignore           # Ignores compiled binaries and custom .in files
├── .gitattributes
├── LICENSE
└── README.md
```

## Dependencies

This project requires both a Fortran compiler for the computational engine and a Python environment for data visualization.

### 1. Fortran Environment
The core numerical solvers are written in Fortran 90. **gfortran** (GNU Fortran) is the recommended compiler.

* **Linux:**
  ```bash
  sudo apt install gfortran
  ```
* **macOS:**
  ```bash
  brew install gcc
  ```
* **Windows:** Available via [MinGW-w64](https://www.mingw-w64.org/) or [MSYS2](https://www.msys2.org/).

### 2. Python Environment
The `plotter.py` script requires **Python 3.8+** and the following libraries for data manipulation and visualization:

* **Pandas:** Used for reading CSV data.
* **NumPy:** Used for array mirroring (symmetry handling) and data parsing.
* **Matplotlib:** Used for generating 2D contour maps and convergence line plots.

You can easily install the required Python dependencies using pip:
```bash
pip install pandas numpy matplotlib
```

## Input File Configuration
The Fortran solver requires a formatted input file to run. This prevents the need to recompile the source code every time you want to test a new grid size or boundary condition. 

A sample input file (`sample_input.in`) is included in this repository. The solver expects the variables to be provided in the exact order shown below. Inline comments starting with `!` are ignored by the parser.

**Format:**
```text
! CFD Solver Input File
! Ideally, do not give very long names to the input file. Keep it 10 characters or less.
! Also, since the name of the output file depends on the input file, you have to manually change the names in plotting script.
! ==========================================
! Grid Parameters
31          ! imax (Number of nodes in X)
41          ! jmax (Number of nodes in Y)

! Physical Dimensions (meters)
0.3         ! L (Length in X)
0.4         ! W (Width in Y)

! Dirichlet Boundary Conditions (Celsius)
40.0        ! T_Bottom 
10.0        ! T_Top    
0.0         ! T_Left   
0.0         ! T_Right

! Performance Checks
.TRUE.      ! Omega optimization
.FALSE.      ! Benchmarking
.FALSE.     ! Quarter Domain

! Omega, these will be used if the omega optimization is set to FALSE
1.835       ! PSOR Omega 
1.777       ! LSOR Omega
```

## Build and Execute
*Open the project root in a terminal and follow these steps:*

1. **Prepare the Output Directory:**
   Ensure there is a folder to catch the CSV data.
   ```bash
   mkdir -p results
   ```
2. **Compile the Fortran core:**
   ```bash
   gfortran -O3 cfd_solvers.f90 main.f90 -o heat_solver
   ```
3. **Run the computational engine:**
   Pass your desired input file as an argument.
   ```bash
   ./heat_solver sample_input.in
   ```
   *(Note: This will populate the `results/` directory with CSV data).*
5. **Run the post-processing pipeline:**
   Generate the contours and line plots from the output data.
   ```bash
   python plotter.py
   ```
   *(Note: The generated .png plots will be saved inside the results/ directory).*
   
## Expected Outputs
Upon successful execution of both the solver and the plotter, the `results/` directory will contain:
* `*_results.csv`: Raw temperature grid data for each problem.
* `*_contour.png`: 2D color maps of the steady-state temperature distribution.
* `*_line_plots.png`: Temperature profiles extracted at 0.05m intervals along the X and Y axes.
* `performance.csv` *(Optional)*: If benchmarking is enabled, contains wall-time and iteration counts for the different solver methods.
