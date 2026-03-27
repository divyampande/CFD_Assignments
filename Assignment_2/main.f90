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
    integer :: iters, csv_id
    real(wp) :: c_time
    
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

    ! PSOR
    call reset_grid(T, imax, jmax)
    call solve_PSOR(T, imax, jmax, dx, dy, 1.5_wp, iters, c_time)
    write(csv_id, '(A10, I20, F20.6, F15.2)') "PSOR", iters, c_time, 1.5_wp

    ! LSOR
    call reset_grid(T, imax, jmax)
    call solve_LSOR(T, imax, jmax, dx, dy, 1.2_wp, iters, c_time)
    write(csv_id, '(A10, I20, F20.6, F15.2)') "LSOR", iters, c_time, 1.2_wp

    close(csv_id)

! ======================================================================
! INTERNAL SUBROUTINES
! ======================================================================
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