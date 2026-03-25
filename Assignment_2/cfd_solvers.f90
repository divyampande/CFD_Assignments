! ======================================================================
! File: cfd_solvers.f90
! Author: Divyam Pandey (AE25M021)
! Course: AM5630 - Foundations of Computational Fluid Dynamics
! Description: Module containing numerical solvers for 2D steady-state 
!              heat conduction (Laplace equation).
! ======================================================================

module kinds_mod
    use, intrinsic :: iso_fortran_env, only: real64 
    implicit none

    ! Define our working precision for real numbers
    integer, parameter :: wp = real64 
end module kinds_mod


module cfd_solvers
    use kinds_mod
    implicit none
contains

    ! UTILITY FUNCTIONS

    ! Calculates the absolute error sum between iterations
    real(wp) function get_error(T_new, T_old, imax, jmax)
        real(wp), intent(in) :: T_new(:,:), T_old(:,:)
        integer, intent(in)  :: imax, jmax
        integer :: i, j
        
        get_error = 0.0_wp 
        
        do j = 2, jmax - 1
            do i = 2, imax - 1
                get_error = get_error + abs(T_new(i,j) - T_old(i,j))
            end do
        end do
    end function get_error

    ! Thomas Algorithm (TDMA) for solving tridiagonal systems
    subroutine tdma(a, b, c, d, x, n)
        integer, intent(in)  :: n
        real(wp), intent(in) :: a(n), b(n), c(n), d(n)
        real(wp), intent(out):: x(n)
        
        x = 0.0_wp ! Placeholder
        ! TODO: Implement Forward Elimination and Backward Substitution
    end subroutine tdma


    ! EXPLICIT SOLVERS (Point-by-Point)

    ! (a) Point Gauss-Seidel Method
    subroutine solve_PGS(T, imax, jmax, dx, dy, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time
        
        iter_count = 0
        comp_time  = 0.0_wp
        ! TODO: Implement the point-by-point Gauss-Seidel sweep
    end subroutine solve_PGS

    ! (c) Point Successive Over-Relaxation
    subroutine solve_PSOR(T, imax, jmax, dx, dy, omega, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy, omega
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time
        
        iter_count = 0
        comp_time  = 0.0_wp
        ! TODO: Implement PGS with the relaxation factor (omega)
    end subroutine solve_PSOR


    ! IMPLICIT SOLVERS (Line-by-Line)

    ! (b) Line Gauss-Seidel Method
    subroutine solve_LGS(T, imax, jmax, dx, dy, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time

        iter_count = 0
        comp_time  = 0.0_wp        
        ! TODO: Set up tridiagonal arrays for each line and call TDMA
    end subroutine solve_LGS

    ! (d) Line Successive Over-Relaxation
    subroutine solve_LSOR(T, imax, jmax, dx, dy, omega, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy, omega
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time

        iter_count = 0
        comp_time  = 0.0_wp
        ! TODO: Implement LGS with relaxation factor
    end subroutine solve_LSOR

    ! (e) Alternating Direction Implicit (ADI) Method
    subroutine solve_ADI(T, imax, jmax, dx, dy, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time

        iter_count = 0
        comp_time  = 0.0_wp
        ! TODO: Implement half-step x-sweep followed by half-step y-sweep
    end subroutine solve_ADI

end module cfd_solvers