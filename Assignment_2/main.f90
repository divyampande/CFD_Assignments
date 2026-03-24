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
    real(wp) :: T_prob1(imax, jmax)
    real(wp) :: T_prob2(imax, jmax)
    
    ! Tracking variables
    integer :: iters
    real(wp) :: c_time
    
    ! Pre-compute grid spacing
    dx = L / real(imax - 1, wp)
    dy = W / real(jmax - 1, wp)


    ! INITIALIZATION ROUTINES
    
    ! TODO: Write a quick subroutine here to initialize T arrays to 0.0
    ! and apply the specific boundary walls (T1, T2, T3, T4) for Prob 1 & 2.


    ! PROBLEM 1 EXECUTION

    print *, "--- Starting Problem 1 Computations ---"
    
    ! TODO: Loop through omegas for PSOR/LSOR to find the optimum
    
    ! TODO: Save the best result to a CSV file


    ! PROBLEM 2 EXECUTION

    print *, "--- Starting Problem 2 Computations ---"
    
    ! TODO: Set up the 0.15m x 0.2m grid for Part B
    ! TODO: Apply symmetry boundary conditions (Neumann BCs)
    ! TODO: Compare computational cost

end program main