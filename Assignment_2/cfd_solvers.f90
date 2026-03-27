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
    ! Equation form: a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i)
    subroutine tdma(a, b, c, d, x, n)
        integer, intent(in)  :: n
        real(wp), intent(in) :: a(n), b(n), c(n), d(n)
        real(wp), intent(out):: x(n)
        
        ! Arrays for the modified coefficients
        real(wp) :: cp(n), dp(n)
        real(wp) :: denom
        integer  :: i
        
        ! Modify the first row coefficients
        cp(1) = c(1) / b(1)
        dp(1) = d(1) / b(1)
        
        ! Loop through the rest of the rows
        do i = 2, n
            denom = b(i) - a(i) * cp(i-1)

            if (i < n) then
                cp(i) = c(i) / denom
            end if
            
            dp(i) = (d(i) - a(i) * dp(i-1)) / denom
        end do
        
        ! The last variable is fully solved now
        x(n) = dp(n)
        
        ! Substitute back up the line to solve the rest
        do i = n - 1, 1, -1
            x(i) = dp(i) - cp(i) * x(i+1)
        end do
        
    end subroutine tdma


    ! EXPLICIT SOLVERS (Point-by-Point)

    ! (a) Point Gauss-Seidel Method
    subroutine solve_PGS(T, imax, jmax, dx, dy, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time
        
        real(wp) :: T_old(imax, jmax)
        real(wp) :: error, beta2, denom
        integer  :: i, j
        integer  :: tick_start, tick_end, tick_rate

        ! Start the timer
        call system_clock(tick_start, tick_rate)

        beta2 = (dx / dy)**2
        denom = 2.0_wp * (1.0_wp + beta2)

        iter_count = 0
        error = 1.0_wp

        do while (error > 0.01_wp)
            iter_count = iter_count + 1
            
            ! Save the current state 
            T_old = T

            ! The Gauss-Seidel Sweep
            do j = 2, jmax - 1
                do i = 2, imax - 1
                    T(i,j) = ( T(i+1,j) + T(i-1,j) + beta2 * (T(i,j+1) + T(i,j-1)) ) / denom
                end do
            end do

            ! Check Convergence
            error = get_error(T, T_old, imax, jmax)

            ! Safety net against infinite loops
            if (iter_count > 50000) then
                print *, "WARNING: PGS did not converge!"
                exit
            end if
        end do

        ! Stop the timer
        call system_clock(tick_end)
        comp_time = real(tick_end - tick_start, wp) / real(tick_rate, wp)

    end subroutine solve_PGS

    ! (c) Point Successive Over-Relaxation (PSOR)
    subroutine solve_PSOR(T, imax, jmax, dx, dy, omega, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy, omega
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time
        
        real(wp) :: T_old(imax, jmax)
        real(wp) :: error, beta2, denom, T_GS
        integer  :: i, j
        integer  :: tick_start, tick_end, tick_rate

        ! Start timer
        call system_clock(tick_start, tick_rate)

        beta2 = (dx / dy)**2
        denom = 2.0_wp * (1.0_wp + beta2)

        iter_count = 0
        error = 1.0_wp 

        do while (error > 0.01_wp)
            iter_count = iter_count + 1
            T_old = T

            ! The PSOR Sweep
            do j = 2, jmax - 1
                do i = 2, imax - 1
                    ! Standard Gauss-Seidel guess
                    T_GS = ( T(i+1,j) + T(i-1,j) + beta2 * (T(i,j+1) + T(i,j-1)) ) / denom
                    
                    ! Relaxation factor
                    T(i,j) = (1.0_wp - omega) * T_old(i,j) + omega * T_GS
                end do
            end do

            ! Check Convergence
            error = get_error(T, T_old, imax, jmax)

            if (iter_count > 50000) then
                print *, "WARNING: PSOR did not converge with omega = ", omega
                exit
            end if
        end do

        ! Stop timer
        call system_clock(tick_end)
        comp_time = real(tick_end - tick_start, wp) / real(tick_rate, wp)

    end subroutine solve_PSOR


    ! IMPLICIT SOLVERS (Line-by-Line)

    ! (b) Line Gauss-Seidel Method (LGS)
    subroutine solve_LGS(T, imax, jmax, dx, dy, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time

        real(wp) :: T_old(imax, jmax)
        real(wp) :: error, beta2
        integer  :: i, j, n
        integer  :: tick_start, tick_end, tick_rate

        ! Arrays for the TDMA
        real(wp) :: a(imax-2), b(imax-2), c(imax-2), d(imax-2), x(imax-2)

        call system_clock(tick_start, tick_rate)

        beta2 = (dx / dy)**2
        n = imax - 2
        
        a = 1.0_wp
        b = -2.0_wp * (1.0_wp + beta2)
        c = 1.0_wp

        iter_count = 0
        error = 1.0_wp

        do while (error > 0.01_wp)
            iter_count = iter_count + 1
            T_old = T

            ! Sweep through the grid horizontally (line by line)
            do j = 2, jmax - 1
                
                ! Build the RHS (d) using the lines above and below
                do i = 2, imax - 1
                    d(i-1) = -beta2 * (T(i, j-1) + T(i, j+1))
                end do

                ! Apply Boundary Conditions
                d(1) = d(1) - T(1, j)       ! Left wall
                d(n) = d(n) - T(imax, j)    ! Right wall

                ! TDMA
                call tdma(a, b, c, d, x, n)

                do i = 2, imax - 1
                    T(i, j) = x(i-1)
                end do
            end do

            ! Check Convergence
            error = get_error(T, T_old, imax, jmax)

            if (iter_count > 50000) then
                print *, "WARNING: LGS did not converge!"
                exit
            end if
        end do

        call system_clock(tick_end)
        comp_time = real(tick_end - tick_start, wp) / real(tick_rate, wp)

    end subroutine solve_LGS

    ! (d) Line Successive Over-Relaxation (LSOR)
    subroutine solve_LSOR(T, imax, jmax, dx, dy, omega, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy, omega
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time

        real(wp) :: T_old(imax, jmax)
        real(wp) :: error, beta2
        integer  :: i, j, n
        integer  :: tick_start, tick_end, tick_rate

        real(wp) :: a(imax-2), b(imax-2), c(imax-2), d(imax-2), x(imax-2)

        call system_clock(tick_start, tick_rate)

        beta2 = (dx / dy)**2
        n = imax - 2 
        
        a = 1.0_wp
        b = -2.0_wp * (1.0_wp + beta2)
        c = 1.0_wp

        iter_count = 0
        error = 1.0_wp

        do while (error > 0.01_wp)
            iter_count = iter_count + 1
            T_old = T

            do j = 2, jmax - 1
                do i = 2, imax - 1
                    d(i-1) = -beta2 * (T(i, j-1) + T(i, j+1))
                end do

                d(1) = d(1) - T(1, j) 
                d(n) = d(n) - T(imax, j)

                call tdma(a, b, c, d, x, n)

                ! LSOR RELAXATION
                do i = 2, imax - 1
                    T(i, j) = (1.0_wp - omega) * T_old(i, j) + omega * x(i-1)
                end do
            end do

            error = get_error(T, T_old, imax, jmax)

            if (iter_count > 50000) then
                print *, "WARNING: LSOR did not converge with omega = ", omega
                exit
            end if
        end do

        call system_clock(tick_end)
        comp_time = real(tick_end - tick_start, wp) / real(tick_rate, wp)

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