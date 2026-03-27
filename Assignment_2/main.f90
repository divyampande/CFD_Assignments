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
    real(wp) :: dx, dy
    
    ! Data Arrays
    real(wp) :: T(imax, jmax)
    
    ! Tracking variables
    integer :: iters, csv_id, w_int
    real(wp) :: c_time, omega
    integer :: min_iters_psor, min_iters_lsor
    real(wp) :: best_omega_psor, best_omega_lsor
    real(wp) :: best_c_time_psor, best_c_time_lsor
    
    ! Pre-compute grid spacing
    dx = L / real(imax - 1, wp)
    dy = W / real(jmax - 1, wp)

    ! INITIALIZATION & BOUNDARY CONDITIONS (Problem 1)
    
    ! Initialize all interior points to 0.0 C
    T = 0.0_wp 
    
    ! Apply Boundary Conditions
    T(:, 1)    = 40.0_wp    ! Bottom wall (T1)
    T(1, :)    = 0.0_wp     ! Left wall (T2)
    T(:, jmax) = 10.0_wp    ! Top wall (T3)
    T(imax, :) = 0.0_wp     ! Right wall (T4)
    
    print *, "--- 2D Heat Conduction Solver Initialized ---"
    print *, "Grid: ", imax, "x", jmax
    print *, "dx = ", dx, "m, dy = ", dy, "m"
    print *, "---------------------------------------------"

    ! Problem 1: Solve using all 5 methods and record performance
    open(newunit=csv_id, file='results/performance.csv', status='replace')
    write(csv_id, '(A10, A20, A20, A15)') "Method", "Iterations", "Comp. Time (s)", "Omega"

    ! PGS
    call reset_grid(T, imax, jmax)
    call solve_PGS(T, imax, jmax, dx, dy, iters, c_time)
    write(csv_id, '(A10, I20, F20.6, A15)') "PGS", iters, c_time, "N/A"

    ! LGS
    call reset_grid(T, imax, jmax)
    call solve_LGS(T, imax, jmax, dx, dy, iters, c_time)
    write(csv_id, '(A10, I20, F20.6, A15)') "LGS", iters, c_time, "N/A"

    ! ADI
    call reset_grid(T, imax, jmax)
    call solve_ADI(T, imax, jmax, dx, dy, iters, c_time)
    write(csv_id, '(A10, I20, F20.6, A15)') "ADI", iters, c_time, "N/A"

    ! RELAXATION FACTOR OPTIMIZERS (Sweeping omega 1.0 to 1.9)
    
    ! Optimize PSOR
    min_iters_psor = 999999
    best_omega_psor = 1.0_wp
    best_c_time_psor = 0.0_wp
    
    do w_int = 1, 19
        omega = real(w_int, wp) / 10.0_wp
        call reset_grid(T, imax, jmax)
        call solve_PSOR(T, imax, jmax, dx, dy, omega, iters, c_time)
        
        ! Print every test to the table
        write(*, '(A10, I20, F20.6, F15.2)') "PSOR", iters, c_time, omega
        
        ! Save the best one
        if (iters < min_iters_psor) then
            min_iters_psor = iters
            best_omega_psor = omega
            best_c_time_psor = c_time
        end if
    end do
    print *, ">> OPTIMAL PSOR: Omega = ", best_omega_psor, min_iters_psor, " iters"
    write(csv_id, '(A10, I20, F20.6, F15.2)') "PSOR", min_iters_psor, best_c_time_psor, best_omega_psor

    ! Optimize LSOR
    min_iters_lsor = 999999
    best_omega_lsor = 1.0_wp
    best_c_time_lsor = 0.0_wp
    
    do w_int = 1, 19
        omega = real(w_int, wp) / 10.0_wp
        call reset_grid(T, imax, jmax)
        call solve_LSOR(T, imax, jmax, dx, dy, omega, iters, c_time)
        
        write(*, '(A10, I20, F20.6, F15.2)') "LSOR", iters, c_time, omega
        
        if (iters < min_iters_lsor) then
            min_iters_lsor = iters
            best_omega_lsor = omega
            best_c_time_lsor = c_time
        end if
    end do
    print *, ">> OPTIMAL LSOR: Omega = ", best_omega_lsor, min_iters_lsor, " iters"
    write(csv_id, '(A10, I20, F20.6, F15.2)') "LSOR", min_iters_lsor, best_c_time_lsor, best_omega_lsor

    close(csv_id)

! INTERNAL SUBROUTINES
contains

    ! Instantly resets the grid to 0.0 and re-applies Problem 1 Boundaries
    subroutine reset_grid(T_array, max_i, max_j)
        real(wp), intent(inout) :: T_array(:,:)
        integer, intent(in) :: max_i, max_j
        
        T_array = 0.0_wp 
        T_array(:, 1)     = 40.0_wp    ! Bottom
        T_array(1, :)     = 0.0_wp     ! Left
        T_array(:, max_j) = 10.0_wp    ! Top
        T_array(max_i, :) = 0.0_wp     ! Right
    end subroutine reset_grid

end program main