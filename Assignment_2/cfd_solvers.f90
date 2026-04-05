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

    ! TDMA for constant-coefficient systems
    subroutine tdma_fast(a, cp, inv_denom, d, dp, x, n)
        integer, intent(in)  :: n
        real(wp), intent(in) :: a                 
        real(wp), intent(in) :: cp(n), inv_denom(n), d(n)
        real(wp), intent(out):: dp(n), x(n)       
        
        integer  :: i
        
        dp(1) = d(1) * inv_denom(1)
        
        do i = 2, n
            dp(i) = (d(i) - a * dp(i-1)) * inv_denom(i)
        end do
        
        x(n) = dp(n)
        do i = n - 1, 1, -1
            x(i) = dp(i) - cp(i) * x(i+1)
        end do
    end subroutine tdma_fast

    ! EXPLICIT SOLVERS (Point-by-Point)
    ! Point Successive Over-Relaxation (PSOR)
    subroutine solve_PSOR(T, imax, jmax, dx, dy, omega, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy, omega
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time
        
        real(wp) :: error, beta2, denom, T_GS, temp_T_old
        real(wp) :: omega_compl, omega_denom
        integer  :: i, j
        integer  :: tick_start, tick_end, tick_rate

        ! Start timer
        call system_clock(tick_start, tick_rate)

        ! Precompute domain constants
        beta2 = (dx / dy)**2
        denom = 2.0_wp * (1.0_wp + beta2)

        omega_compl = 1.0_wp - omega
        omega_denom = omega / denom

        iter_count = 0
        error = 1.0_wp ! Initialize > 0.01 to enter the loop

        do while (error > 0.01_wp)
            iter_count = iter_count + 1
            error = 0.0_wp ! Reset L1 error sum for the current iteration

            ! The Fused PSOR Sweep
            do j = 2, jmax - 1
                do i = 2, imax - 1
                    temp_T_old = T(i,j)
                    
                    ! Standard Gauss-Seidel sum
                    T_GS = T(i+1,j) + T(i-1,j) + beta2 * (T(i,j+1) + T(i,j-1))
                    
                    ! Apply Relaxation factor
                    T(i,j) = omega_compl * temp_T_old + omega_denom * T_GS
                    
                    ! Error calculation
                    error = error + abs(T(i,j) - temp_T_old)
                end do
            end do

            ! DIVERGENCE TRAP (Catches NaN and Infinity)
            if (error /= error .or. error > 1.0e10_wp) then
                iter_count = 999999
                exit
            end if

            ! Failsafe exit
            if (iter_count > 50000) then
                print *, "WARNING: PSOR did not converge with omega = ", omega, " Final error = ", error
                exit
            end if
        end do

        ! Stop timer
        call system_clock(tick_end)
        comp_time = real(tick_end - tick_start, wp) / real(tick_rate, wp)

    end subroutine solve_PSOR
    
    ! Line Successive Over-Relaxation (LSOR)
    subroutine solve_LSOR(T, imax, jmax, dx, dy, omega, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy, omega
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time

        real(wp) :: error, beta2, temp_T_old
        real(wp) :: a_val, b_val, c_val, denom_temp
        real(wp) :: b_base, relax_val
        integer  :: i, j, n
        integer  :: tick_start, tick_end, tick_rate

        real(wp) :: d(imax-2), x(imax-2), dp(imax-2)
        real(wp) :: cp(imax-2), inv_denom(imax-2) 

        call system_clock(tick_start, tick_rate)

        beta2 = (dx / dy)**2
        n = imax - 2 

        ! MATRIX SPLITTING GEOMETRY
        a_val = 1.0_wp
        c_val = 1.0_wp
        b_base = -2.0_wp * (1.0_wp + beta2)

        ! Divide the main diagonal by omega
        b_val = b_base / omega
        relax_val = b_base * (1.0_wp - omega) / omega

        inv_denom(1) = 1.0_wp / b_val
        cp(1) = c_val * inv_denom(1)
        
        do i = 2, n
            denom_temp = b_val - a_val * cp(i-1)
            inv_denom(i) = 1.0_wp / denom_temp  
            if (i < n) cp(i) = c_val * inv_denom(i)
        end do

        iter_count = 0
        error = 1.0_wp

        do while (error > 0.01_wp)
            iter_count = iter_count + 1
            error = 0.0_wp

            do j = 2, jmax - 1
                do i = 2, imax - 1
                    d(i-1) = -beta2 * (T(i, j-1) + T(i, j+1)) + relax_val * T(i, j)
                end do
                
                d(1) = d(1) - T(1, j) 
                d(n) = d(n) - T(imax, j)

                ! Solve using TDMA
                call tdma_fast(a_val, cp, inv_denom, d, dp, x, n)

                ! Error calculation
                do i = 2, imax - 1
                    temp_T_old = T(i, j)
                    T(i, j) = x(i-1)
                    error = error + abs(T(i, j) - temp_T_old)
                end do
            end do

            ! DIVERGENCE TRAP (Catches NaN and Infinity)
            if (error /= error .or. error > 1.0e10_wp) then
                iter_count = 999999
                exit
            end if

            if (iter_count > 50000) then
                print *, "WARNING: LSOR did not converge. Final error =", error
                exit
            end if
        end do

        call system_clock(tick_end)
        comp_time = real(tick_end - tick_start, wp) / real(tick_rate, wp)

    end subroutine solve_LSOR

    ! Alternating Direction Implicit with Relaxation (ADIR)
    subroutine solve_ADIR(T, imax, jmax, dx, dy, omega, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy, omega
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time

        real(wp) :: T_old(imax, jmax) 
        real(wp) :: error, beta2, alpha2, denom_temp
        integer  :: i, j, nx, ny
        integer  :: tick_start, tick_end, tick_rate

        ! Matrix Splitting Variables
        real(wp) :: bx_base, by_base, relax_x, relax_y

        ! X-Sweep TDMA Arrays
        real(wp) :: d_x(imax-2), x_val(imax-2), dp_x(imax-2)
        real(wp) :: cp_x(imax-2), inv_denom_x(imax-2)
        real(wp) :: ax_val, bx_val, cx_val

        ! Y-Sweep TDMA Arrays
        real(wp) :: d_y(jmax-2), y_val(jmax-2), dp_y(jmax-2)
        real(wp) :: cp_y(jmax-2), inv_denom_y(jmax-2)
        real(wp) :: ay_val, by_val, cy_val

        call system_clock(tick_start, tick_rate)

        beta2 = (dx / dy)**2
        alpha2 = (dy / dx)**2
        nx = imax - 2 
        ny = jmax - 2

        ! MATRIX SPLITTING
        ax_val = 1.0_wp; cx_val = 1.0_wp
        bx_base = -2.0_wp * (1.0_wp + beta2)
        
        ay_val = 1.0_wp; cy_val = 1.0_wp
        by_base = -2.0_wp * (1.0_wp + alpha2)

        ! Divide the main diagonals by omega
        bx_val = bx_base / omega
        by_val = by_base / omega
        relax_x = bx_base * (1.0_wp - omega) / omega
        relax_y = by_base * (1.0_wp - omega) / omega

        inv_denom_x(1) = 1.0_wp / bx_val
        cp_x(1) = cx_val * inv_denom_x(1)
        do i = 2, nx
            denom_temp = bx_val - ax_val * cp_x(i-1)
            inv_denom_x(i) = 1.0_wp / denom_temp
            if (i < nx) cp_x(i) = cx_val * inv_denom_x(i)
        end do

        inv_denom_y(1) = 1.0_wp / by_val
        cp_y(1) = cy_val * inv_denom_y(1)
        do j = 2, ny
            denom_temp = by_val - ay_val * cp_y(j-1)
            inv_denom_y(j) = 1.0_wp / denom_temp
            if (j < ny) cp_y(j) = cy_val * inv_denom_y(j)
        end do

        iter_count = 0
        error = 1.0_wp
        
        do while (error > 0.01_wp)
            iter_count = iter_count + 1
            T_old = T
            error = 0.0_wp

            ! HALF-STEP 1: Implicit X-Sweep
            do j = 2, jmax - 1
                do i = 2, imax - 1
                    d_x(i-1) = -beta2 * (T(i, j-1) + T(i, j+1)) + relax_x * T(i, j)
                end do

                d_x(1)  = d_x(1)  - T(1, j) 
                d_x(nx) = d_x(nx) - T(imax, j)

                call tdma_fast(ax_val, cp_x, inv_denom_x, d_x, dp_x, x_val, nx)

                do i = 2, imax - 1
                    T(i, j) = x_val(i-1)
                end do
            end do

            ! HALF-STEP 2: Implicit Y-Sweep
            do i = 2, imax - 1
                do j = 2, jmax - 1
                    d_y(j-1) = -alpha2 * (T(i-1, j) + T(i+1, j)) + relax_y * T(i, j)
                end do

                d_y(1)  = d_y(1)  - T(i, 1) 
                d_y(ny) = d_y(ny) - T(i, jmax)

                call tdma_fast(ay_val, cp_y, inv_denom_y, d_y, dp_y, y_val, ny)

                do j = 2, jmax - 1
                    T(i, j) = y_val(j-1) 
                    error = error + abs(T(i, j) - T_old(i, j))
                end do
            end do

            ! DIVERGENCE TRAP
            if (error /= error .or. error > 1.0e10_wp) then
                iter_count = 999999
                exit
            end if

            if (iter_count > 50000) then
                print *, "WARNING: ADIR did not converge. Final error =", error
                exit
            end if
        end do

        call system_clock(tick_end)
        comp_time = real(tick_end - tick_start, wp) / real(tick_rate, wp)

    end subroutine solve_ADIR

    ! LSOR for Quarter Domain (Symmetry Boundaries)
    subroutine solve_LSOR_sym(T, imax, jmax, dx, dy, omega, iter_count, comp_time)
        real(wp), intent(inout) :: T(:,:)
        integer, intent(in)     :: imax, jmax
        real(wp), intent(in)    :: dx, dy, omega
        integer, intent(out)    :: iter_count
        real(wp), intent(out)   :: comp_time

        real(wp) :: error, beta2, temp_T_old
        real(wp) :: a_val, b_val, c_val, denom_temp
        real(wp) :: b_base, relax_val
        integer  :: i, j, nx
        integer  :: tick_start, tick_end, tick_rate

        ! TDMA Precomputed Arrays and Work Arrays
        real(wp) :: d(imax-1), x_val(imax-1), dp(imax-1)
        real(wp) :: cp(imax-1), inv_denom(imax-1)

        call system_clock(tick_start, tick_rate)

        beta2 = (dx / dy)**2
        nx = imax - 1 

        a_val = 1.0_wp
        c_val = 1.0_wp
        b_base = -2.0_wp * (1.0_wp + beta2)

        ! Divide the main diagonal by omega
        b_val = b_base / omega

        ! Precalculate the RHS relaxation
        relax_val = b_base * (1.0_wp - omega) / omega

        ! PRECOMPUTE TDMA (Handling the Symmetry Wall!)
        inv_denom(1) = 1.0_wp / b_val
        cp(1) = c_val * inv_denom(1)
        
        ! Standard internal nodes
        do i = 2, nx - 1
            denom_temp = b_val - a_val * cp(i-1)
            inv_denom(i) = 1.0_wp / denom_temp
            cp(i) = c_val * inv_denom(i)
        end do
        
        ! Right Symmetry Wall Ghost Node (Index nx)
        ! Here, a = 2.0 and c = 0.0. The diagonal is still b_val.
        denom_temp = b_val - 2.0_wp * cp(nx-1)
        inv_denom(nx) = 1.0_wp / denom_temp
        cp(nx) = 0.0_wp

        iter_count = 0
        error = 1.0_wp

        do while (error > 0.01_wp)
            iter_count = iter_count + 1
            error = 0.0_wp

            ! STANDARD GRID LINES (j = 2 to jmax - 1)
            do j = 2, jmax - 1
                do i = 2, imax
                    d(i-1) = -beta2 * (T(i, j-1) + T(i, j+1)) + relax_val * T(i, j)
                end do
                
                ! Left Dirichlet Boundary (Fixed at 0 C)
                d(1) = d(1) - T(1, j) 

                ! INLINED FAST TDMA
                dp(1) = d(1) * inv_denom(1)
                do i = 2, nx - 1
                    dp(i) = (d(i) - a_val * dp(i-1)) * inv_denom(i)
                end do
                dp(nx) = (d(nx) - 2.0_wp * dp(nx-1)) * inv_denom(nx) ! Right wall uses a=2.0

                x_val(nx) = dp(nx)
                do i = nx - 1, 1, -1
                    x_val(i) = dp(i) - cp(i) * x_val(i+1)
                end do

                ! Error
                do i = 2, imax
                    temp_T_old = T(i, j)
                    T(i, j) = x_val(i-1)
                    error = error + abs(T(i, j) - temp_T_old)
                end do
            end do

            ! TOP SYMMETRY LINE ONLY (j = jmax)
            j = jmax
            
            ! Build RHS using Top Wall Ghost Node logic (2.0 * T(i, j-1))
            do i = 2, imax
                d(i-1) = -beta2 * (2.0_wp * T(i, j-1)) + relax_val * T(i, j)
            end do
            d(1) = d(1) - T(1, j) 

            ! INLINED FAST TDMA
            dp(1) = d(1) * inv_denom(1)
            do i = 2, nx - 1
                dp(i) = (d(i) - a_val * dp(i-1)) * inv_denom(i)
            end do
            dp(nx) = (d(nx) - 2.0_wp * dp(nx-1)) * inv_denom(nx)

            x_val(nx) = dp(nx)
            do i = nx - 1, 1, -1
                x_val(i) = dp(i) - cp(i) * x_val(i+1)
            end do

            ! Error
            do i = 2, imax
                temp_T_old = T(i, j)
                T(i, j) = x_val(i-1)
                error = error + abs(T(i, j) - temp_T_old)
            end do

            ! DIVERGENCE TRAP
            if (error /= error .or. error > 1.0e10_wp) then
                iter_count = 999999
                exit
            end if

            if (iter_count > 50000) then
                print *, "WARNING: LSOR_sym did not converge. Final error =", error
                exit
            end if
        end do

        call system_clock(tick_end)
        comp_time = real(tick_end - tick_start, wp) / real(tick_rate, wp)

    end subroutine solve_LSOR_sym

end module cfd_solvers