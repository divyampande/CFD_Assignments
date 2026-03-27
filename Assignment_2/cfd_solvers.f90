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

        real(wp) :: T_old(imax, jmax)
        real(wp) :: T_half(imax, jmax) ! Array for the half-step
        real(wp) :: error, beta2, alpha2
        integer  :: i, j, nx, ny
        integer  :: tick_start, tick_end, tick_rate

        ! Arrays for the X-Sweep
        real(wp) :: ax(imax-2), bx(imax-2), cx(imax-2), rhs_x(imax-2), x_val(imax-2)
        ! Arrays for the Y-Sweep
        real(wp) :: ay(jmax-2), by(jmax-2), cy(jmax-2), rhs_y(jmax-2), y_val(jmax-2)

        call system_clock(tick_start, tick_rate)

        ! Geometric ratios
        beta2 = (dx / dy)**2
        alpha2 = (dy / dx)**2
        nx = imax - 2 
        ny = jmax - 2
        ax = 1.0_wp;  bx = -2.0_wp * (1.0_wp + beta2);  cx = 1.0_wp
        ay = 1.0_wp;  by = -2.0_wp * (1.0_wp + alpha2); cy = 1.0_wp

        iter_count = 0
        error = 1.0_wp
        T_half = T

        do while (error > 0.01_wp)
            iter_count = iter_count + 1
            T_old = T

            ! HALF-STEP 1: Implicit X-Sweep
            do j = 2, jmax - 1
                do i = 2, imax - 1
                    rhs_x(i-1) = -beta2 * (T(i, j-1) + T(i, j+1))
                end do

                ! Left and Right boundary conditions
                rhs_x(1)  = rhs_x(1)  - T_half(1, j) 
                rhs_x(nx) = rhs_x(nx) - T_half(imax, j)

                call tdma(ax, bx, cx, rhs_x, x_val, nx)

                do i = 2, imax - 1
                    T_half(i, j) = x_val(i-1)
                end do
            end do

            ! HALF-STEP 2: Implicit Y-Sweep
            do i = 2, imax - 1
                do j = 2, jmax - 1
                    rhs_y(j-1) = -alpha2 * (T_half(i-1, j) + T_half(i+1, j))
                end do

                ! Bottom and Top boundary conditions
                rhs_y(1)  = rhs_y(1)  - T(i, 1) 
                rhs_y(ny) = rhs_y(ny) - T(i, jmax)

                call tdma(ay, by, cy, rhs_y, y_val, ny)

                do j = 2, jmax - 1
                    T(i, j) = y_val(j-1)
                end do
            end do

            ! Check Convergence (Compare fully completed step against old step)
            error = get_error(T, T_old, imax, jmax)

            if (iter_count > 50000) then
                print *, "WARNING: ADI did not converge!"
                exit
            end if
        end do

        call system_clock(tick_end)
        comp_time = real(tick_end - tick_start, wp) / real(tick_rate, wp)

    end subroutine solve_ADI

    ! (f) LSOR for Quarter Domain (Symmetry Boundaries)
    subroutine solve_LSOR_sym(T, imax, jmax, dx, dy, omega, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy, omega
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time

        real(wp) :: T_old(imax, jmax), T_star(imax-1)
        real(wp) :: error, beta2
        integer  :: i, j, nx
        integer  :: tick_start, tick_end, tick_rate

        ! TDMA Arrays (Size is now imax - 1 because we solve for the right wall too!)
        real(wp) :: a(imax-1), b(imax-1), c(imax-1), rhs(imax-1)

        call system_clock(tick_start, tick_rate)

        beta2 = (dx / dy)**2
        nx = imax - 1 

        a = 1.0_wp
        b = -2.0_wp * (1.0_wp + beta2)
        c = 1.0_wp

        ! Apply Ghost Node modification for the Right Symmetry Wall (i = imax)
        a(nx) = 2.0_wp
        c(nx) = 0.0_wp

        iter_count = 0
        error = 1.0_wp

        do while (error > 0.01_wp)
            iter_count = iter_count + 1
            T_old = T

            ! Loop all the way to jmax (Top Symmetry Wall)
            do j = 2, jmax
                
                ! Build RHS
                do i = 2, imax
                    ! If we are at the top symmetry line, use the Y-direction ghost node
                    if (j == jmax) then
                        rhs(i-1) = -beta2 * (2.0_wp * T(i, j-1))
                    else
                        rhs(i-1) = -beta2 * (T(i, j-1) + T(i, j+1))
                    end if
                end do

                ! Apply Left Dirichlet Boundary (Fixed at 0 C)
                rhs(1) = rhs(1) - T(1, j) 

                call tdma(a, b, c, rhs, T_star, nx)

                ! Apply Over-Relaxation
                do i = 2, imax
                    T(i, j) = omega * T_star(i-1) + (1.0_wp - omega) * T_old(i, j)
                end do
            end do

            error = 0.0_wp
            do j = 2, jmax
                do i = 2, imax
                    error = error + abs(T(i, j) - T_old(i, j))
                end do
            end do

            if (iter_count > 50000) exit
        end do

        call system_clock(tick_end)
        comp_time = real(tick_end - tick_start, wp) / real(tick_rate, wp)

    end subroutine solve_LSOR_sym

end module cfd_solvers