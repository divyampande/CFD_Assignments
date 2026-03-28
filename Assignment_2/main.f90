! ======================================================================
! File: main.f90
! Author: Divyam Pandey (AE25M021)
! Description: Main execution script for Computer Assignment 2. 
!              Handles grid initialization, boundary conditions, and I/O.
! ======================================================================

program main
    use kinds_mod
    use cfd_solvers
    implicit none

    ! Grid Parameters (Problem 1)
    integer, parameter :: imax = 31, jmax = 41
    real(wp), parameter :: L = 0.3_wp, W = 0.4_wp
    real(wp), parameter :: BC(4) = [40.0_wp, 10.0_wp, 0.0_wp, 0.0_wp]  ! Bottom, Top, Left, Right
    real(wp), parameter :: BC_P2(4) = [40.0_wp, 40.0_wp, 0.0_wp, 0.0_wp]  ! Bottom, Top, Left, Right
    real(wp) :: dx, dy
    
    ! Data Arrays
    real(wp) :: T(imax, jmax)

    ! Data Arrays for Quarter Domain
    integer, parameter :: imax_q = 16, jmax_q = 21
    real(wp) :: T_q(imax_q, jmax_q)
    
    ! Tracking variables
    integer :: iters, csv_id, w_int
    real(wp) :: c_time, omega
    integer :: min_iters_psor, min_iters_lsor
    real(wp) :: best_omega_psor, best_omega_lsor
    
    ! Timing variables
    integer, parameter :: N_REPEAT = 200

    ! Pre-compute grid spacing
    dx = L / real(imax - 1, wp)
    dy = W / real(jmax - 1, wp)

    ! INITIALIZATION & BOUNDARY CONDITIONS (Problem 1)
    call reset_grid(T, imax, jmax, BC)
    
    print *, "--- 2D Heat Conduction Solver Initialized ---"
    print *, "Grid: ", imax, "x", jmax
    print *, "dx = ", dx, "m, dy = ", dy, "m"
    print *, "---------------------------------------------"

    ! RELAXATION FACTOR OPTIMIZERS (Sweeping omega 1.0 to 1.9)
    
    ! Optimize PSOR
    min_iters_psor = 999999
    best_omega_psor = 1.0_wp
    
    ! Since I already ran between 1 to 1.9 with a step of 0.001,
    ! I already know the optimal omega is around 1.835, so I will just sweep,
    ! between 1.8 and 1.85 to save time. If I had to do the entire range, I would just change the limits of this loop.
    do w_int = 1800, 1850
        omega = real(w_int, wp) / 1000.0_wp
        call reset_grid(T, imax, jmax, BC)
        call solve_PSOR(T, imax, jmax, dx, dy, omega, iters, c_time)
        
        ! Print every test to the table
        write(*, '(A10, I20, F20.6, F15.3)') "PSOR", iters, c_time, omega
        
        ! Save the best one
        if (iters < min_iters_psor) then
            min_iters_psor = iters
            best_omega_psor = omega
        end if
    end do
    print *, ">> OPTIMAL PSOR: Omega = ", best_omega_psor, min_iters_psor, " iters"

    ! Optimize LSOR
    min_iters_lsor = 999999
    best_omega_lsor = 1.0_wp
    
    ! Since I already ran between 1 to 1.9 with a step of 0.001,
    ! I already know the optimal omega is around 1.777, so I will just sweep,
    ! between 1.75 and 1.8 to save time. If I had to do the entire range, I would just change the limits of this loop.
    do w_int = 1750, 1800
        omega = real(w_int, wp) / 1000.0_wp
        call reset_grid(T, imax, jmax, BC)
        call solve_LSOR(T, imax, jmax, dx, dy, omega, iters, c_time)
        
        write(*, '(A10, I20, F20.6, F15.3)') "LSOR", iters, c_time, omega
        
        if (iters < min_iters_lsor) then
            min_iters_lsor = iters
            best_omega_lsor = omega
        end if
    end do
    print *, ">> OPTIMAL LSOR: Omega = ", best_omega_lsor, min_iters_lsor, " iters"

    ! DATA EXPORT
    ! We will use LSOR to generate the final steady-state field since all 
    ! methods converge to the exact same physical answer.
    print *, "--- Running Problem 1 ---"
    call reset_grid(T, imax, jmax, BC)
    call solve_LSOR(T, imax, jmax, dx, dy, best_omega_lsor, iters, c_time)
    call export_to_csv("results/prob1_results.csv", T, imax, jmax, dx, dy)

    ! Problem 2 Case A: Change top BC to 40.0 and solve again
    ! We will just use the best LSOR from Problem 1 to generate the final field for Problem 2a as well.

    print *, "--- Running Problem 2 Cases ---"
    call reset_grid(T, imax, jmax, BC_P2)
    call solve_LSOR(T, imax, jmax, dx, dy, best_omega_lsor, iters, c_time)
    call export_to_csv("results/prob2_caseA_results.csv", T, imax, jmax, dx, dy)

    ! Problem 2 Case B (Quarter Domain Symmetry)
    T_q = 0.0_wp
    
    ! Apply Dirichlet Boundaries
    T_q(:, 1) = 40.0_wp     ! Bottom Wall
    T_q(1, :) = 0.0_wp      ! Left Wall
    ! T_q(1, 1) = (40.0_wp + 0.0_wp) / 2.0_wp  ! Corner Averaging
    
    ! We will use the best_omega_lsor
    call solve_LSOR_sym(T_q, imax_q, jmax_q, dx, dy, best_omega_lsor, iters, c_time)
    call export_to_csv("results/prob2_caseB_results.csv", T_q, imax_q, jmax_q, dx, dy)
    
    ! ! TIMING ANALYSIS 
    ! ! For a more rigorous timing analysis, we could repeat the solver multiple times and average the time.
    ! ! This is especially useful for very fast methods where timing can be noisy.
    ! ! Comment these out if you just want to run the solvers without benchmarking. 
    ! ! Recommended to run the benchmarks separately since they will take a long time to execute.
    ! open(newunit=csv_id, file='results/performance.csv', status='replace')
    ! write(csv_id, '(A10, A, A20, A, A20, A, A15)') "Method", ",", "Iterations", ",", "Comp. Time (ms)", ",", "Omega"
    ! close(csv_id)

    ! ! --- Benchmarking Problem 1 ---
    ! print *, "--- Running Problem 1 Benchmarks ---"
    ! call benchmark_solver("PGS")
    ! call benchmark_solver("LGS")
    ! call benchmark_solver("ADI")
    ! call benchmark_solver("PSOR", best_omega_psor)
    ! call benchmark_solver("LSOR", best_omega_lsor)
    
    ! ! --- Benchmarking Problem 2 ---
    ! print *, "--- Running Problem 2 Benchmarks ---"
    ! call benchmark_solver("LSOR_P2a", best_omega_lsor)
    ! call benchmark_solver("LSOR_P2b", best_omega_lsor)

! INTERNAL SUBROUTINES
contains

    ! Instantly resets the grid to 0.0 and re-applies Problem 1 Boundaries
    subroutine reset_grid(T_array, max_i, max_j, BC_list)
        real(wp), intent(inout) :: T_array(:,:)
        integer, intent(in) :: max_i, max_j
        real(wp), intent(in) :: BC_list(4)
        
        T_array = 0.0_wp 
        ! Apply BCs
        T_array(:, 1)     = BC_list(1)  ! Bottom
        T_array(:, max_j) = BC_list(2)  ! Top
        T_array(1, :)     = BC_list(3)  ! Left
        T_array(max_i, :) = BC_list(4)  ! Right

        ! ! Corners (Average of adjacent walls)
        ! T_array(1, 1)         = (BC_list(1) +  BC_list(3)) / 2.0_wp  ! Bottom-Left
        ! T_array(max_i, 1)     = (BC_list(1) +  BC_list(4)) / 2.0_wp  ! Bottom-Right
        ! T_array(1, max_j)     = (BC_list(2) +  BC_list(3)) / 2.0_wp  ! Top-Left
        ! T_array(max_i, max_j) = (BC_list(2) +  BC_list(4)) / 2.0_wp  ! Top-Right
    end subroutine reset_grid

    ! Writes the 2D grid into a 1D X, Y, T format for Python
    subroutine export_to_csv(filename, T_array, max_i, max_j, d_x, d_y)
        character(len=*), intent(in) :: filename
        real(wp), intent(in)         :: T_array(:,:), d_x, d_y
        integer, intent(in)          :: max_i, max_j
        
        integer  :: i, j
        real(wp) :: x_pos, y_pos

        open(newunit=csv_id, file=filename, status='replace')
        write(csv_id, '(A)') "X,Y,T"
        
        do j = 1, max_j
            y_pos = real(j - 1, wp) * d_y
            do i = 1, max_i
                x_pos = real(i - 1, wp) * d_x
                ! Dump coordinates and temperature
                write(csv_id, '(F10.6, A, F10.6, A, F10.6)') x_pos, ",", y_pos, ",", T_array(i,j)
            end do
        end do
        
        close(csv_id)
        print *, ">> Steady-state data exported to: ", filename
    end subroutine export_to_csv

    ! Benchmarker
    subroutine benchmark_solver(solver_name, omega_val)
        character(len=*), intent(in)   :: solver_name
        real(wp), intent(in), optional :: omega_val
        
        real(wp) :: t_start, t_end, avg_time
        integer  :: rep
        real(wp) :: dummy_c_time

        if (present(omega_val)) then
            omega = omega_val
        else
            omega = 1.0_wp
        end if

        avg_time = 0.0_wp
        call cpu_time(t_start)
        
        do rep = 1, N_REPEAT
            select case (solver_name)
                case ("PGS")
                    call reset_grid(T, imax, jmax, BC)
                    call solve_PGS(T, imax, jmax, dx, dy, iters, dummy_c_time)
                case ("LGS")
                    call reset_grid(T, imax, jmax, BC)
                    call solve_LGS(T, imax, jmax, dx, dy, iters, dummy_c_time)
                case ("ADI")
                    call reset_grid(T, imax, jmax, BC)
                    call solve_ADI(T, imax, jmax, dx, dy, iters, dummy_c_time)
                case ("PSOR")
                    call reset_grid(T, imax, jmax, BC)
                    call solve_PSOR(T, imax, jmax, dx, dy, omega, iters, dummy_c_time)
                case ("LSOR")
                    call reset_grid(T, imax, jmax, BC)
                    call solve_LSOR(T, imax, jmax, dx, dy, omega, iters, dummy_c_time)
                
                ! PROBLEM 2 CASES
                case ("LSOR_P2a")
                    call reset_grid(T, imax, jmax, BC_P2) ! Use Problem 2 Boundaries
                    call solve_LSOR(T, imax, jmax, dx, dy, omega, iters, dummy_c_time)
                case ("LSOR_P2b")
                    T_q = 0.0_wp
                    T_q(:, 1) = 40.0_wp     ! Bottom Wall
                    T_q(1, :) = 0.0_wp      ! Left Wall
                    call solve_LSOR_sym(T_q, imax_q, jmax_q, dx, dy, omega, iters, dummy_c_time)
                
                case default
                    print *, "ERROR: Unknown solver passed to benchmark routine!"
                    stop
            end select
        end do
        
        call cpu_time(t_end)
        avg_time = (t_end - t_start) / real(N_REPEAT, wp) * 1000.0_wp  ! Convert to milliseconds

        open(newunit=csv_id, file='results/performance.csv', status='old', position='append', action='write')
        if (present(omega_val)) then
            write(csv_id, '(A10, A, I20, A, F20.3, A, F15.3)') solver_name, ",", iters, ",", avg_time, ",", omega
        else
            write(csv_id, '(A10, A, I20, A, F20.3, A, A15)') solver_name, ",", iters, ",", avg_time, ",", "N/A"
        end if
        close(csv_id)
        
        print *, "Benchmarked: ", solver_name, " | Avg Time: ", avg_time, "ms | Iters: ", iters
    end subroutine benchmark_solver
end program main