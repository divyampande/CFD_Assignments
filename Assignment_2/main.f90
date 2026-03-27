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
    call reset_grid(T, imax, jmax, BC)
    
    print *, "--- 2D Heat Conduction Solver Initialized ---"
    print *, "Grid: ", imax, "x", jmax
    print *, "dx = ", dx, "m, dy = ", dy, "m"
    print *, "---------------------------------------------"

    ! Problem 1: Solve using all 5 methods and record performance
    open(newunit=csv_id, file='results/performance.csv', status='replace')
    write(csv_id, '(A10, A20, A20, A15)') "Method", "Iterations", "Comp. Time (s)", "Omega"

    ! PGS
    call reset_grid(T, imax, jmax, BC)
    call solve_PGS(T, imax, jmax, dx, dy, iters, c_time)
    write(csv_id, '(A10, I20, F20.6, A15)') "PGS", iters, c_time, "N/A"

    ! LGS
    call reset_grid(T, imax, jmax, BC)
    call solve_LGS(T, imax, jmax, dx, dy, iters, c_time)
    write(csv_id, '(A10, I20, F20.6, A15)') "LGS", iters, c_time, "N/A"

    ! ADI
    call reset_grid(T, imax, jmax, BC)
    call solve_ADI(T, imax, jmax, dx, dy, iters, c_time)
    write(csv_id, '(A10, I20, F20.6, A15)') "ADI", iters, c_time, "N/A"

    ! RELAXATION FACTOR OPTIMIZERS (Sweeping omega 1.0 to 1.9)
    
    ! Optimize PSOR
    min_iters_psor = 999999
    best_omega_psor = 1.0_wp
    best_c_time_psor = 0.0_wp
    
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
            best_c_time_psor = c_time
        end if
    end do
    print *, ">> OPTIMAL PSOR: Omega = ", best_omega_psor, min_iters_psor, " iters"
    write(csv_id, '(A10, I20, F20.6, F15.3)') "PSOR", min_iters_psor, best_c_time_psor, best_omega_psor

    ! Optimize LSOR
    min_iters_lsor = 999999
    best_omega_lsor = 1.0_wp
    best_c_time_lsor = 0.0_wp
    
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
            best_c_time_lsor = c_time
        end if
    end do
    print *, ">> OPTIMAL LSOR: Omega = ", best_omega_lsor, min_iters_lsor, " iters"
    write(csv_id, '(A10, I20, F20.6, F15.3)') "LSOR", min_iters_lsor, best_c_time_lsor, best_omega_lsor

    close(csv_id)

    ! DATA EXPORT
    ! We will use LSOR to generate the final steady-state field since all 
    ! methods converge to the exact same physical answer.
    call reset_grid(T, imax, jmax, BC)
    call solve_LSOR(T, imax, jmax, dx, dy, best_omega_lsor, iters, c_time)
    call export_to_csv("results/prob1_results.csv", T, imax, jmax, dx, dy)

    ! Problem 2 Case a: Change top BC to 40.0 and solve again
    ! We will just use the best LSOR from Problem 1 to generate the final field for Problem 2a as well.
    call reset_grid(T, imax, jmax, BC_P2)
    call solve_LSOR(T, imax, jmax, dx, dy, best_omega_lsor, iters, c_time)
    open(newunit=csv_id, file='results/performance.csv', status='old', position='append', action='write')
    write(csv_id, '(A10, I20, F20.6, F15.3)') "LSOR_P2a", iters, c_time, best_omega_lsor
    close(csv_id)
    call export_to_csv("results/prob2_caseA_results.csv", T, imax, jmax, dx, dy)

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
end program main