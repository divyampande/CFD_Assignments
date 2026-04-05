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
    integer :: imax, jmax
    real(wp) :: L, W
    real(wp) :: BC(4)        ! Bottom, Top, Left, Right
    real(wp) :: dx, dy

    ! Data Arrays
    real(wp), allocatable :: T(:,:)

    ! Data Arrays for Quarter Domain
    integer :: imax_q, jmax_q
    real(wp), allocatable :: T_q(:,:)

    ! Tracking & Timing variables
    integer :: iters, csv_id, w_int
    real(wp) :: c_time, omega
    integer :: min_iters_psor, min_iters_lsor, min_iters_adir
    real(wp) :: best_omega_psor, best_omega_lsor, best_omega_adir

    ! Constants
    integer, parameter :: N_REPEAT = 200

    ! Logical Variables
    Logical :: benchmark_mode, optimize_omega, quarter_domain

    ! File I/O variables
    integer :: file_unit
    character(len=256) :: filename
    character(len=256) :: base_name   ! Holds the name without extension
    character(len=256) :: out_file    ! Holds the final CSV path
    integer :: dot_idx, slash_idx     ! For string manipulation

    ! Read the input filename from command line arguments
    call get_command_argument(1, filename)

    if (trim(filename) == '') then
        print *, "ERROR: No input file provided!"
        print *, "Usage: ./cfd_solver <input_file.in>"
        stop
    end if

    ! Output file naming logic:
    ! Find the last slash to strip any folder paths (e.g., "inputs/prob1.in" to "prob1.in")
    slash_idx = index(trim(filename), '/', back=.true.)
    if (slash_idx > 0) then
        base_name = filename(slash_idx+1 : len_trim(filename))
    else
        base_name = trim(filename)
    end if
    
    ! Find the last dot to strip the extension (e.g., "prob1.in" to "prob1")
    dot_idx = index(trim(base_name), '.', back=.true.)
    if (dot_idx > 0) then
        base_name = base_name(1 : dot_idx-1)
    end if
    
    ! Glue it all together into the final relative path
    out_file = "results/" // trim(base_name) // "_results.csv"

    print *, "Reading parameters from: ", trim(filename)

    ! Read input file and extract parameters
    open(newunit=file_unit, file=trim(filename), status='old', action='read')

    read(file_unit, *) ! Skip header
    read(file_unit, *) ! Skip naming comment
    read(file_unit, *) ! Skip plotting comment
    read(file_unit, *) ! Skip separator
    read(file_unit, *) ! Skip section title
    read(file_unit, *) imax
    read(file_unit, *) jmax

    read(file_unit, *) ! Skip blank line
    read(file_unit, *) ! Skip section title
    read(file_unit, *) L
    read(file_unit, *) W

    read(file_unit, *) ! Skip blank line
    read(file_unit, *) ! Skip section title
    read(file_unit, *) BC(1) ! Bottom
    read(file_unit, *) BC(2) ! Top
    read(file_unit, *) BC(3) ! Left
    read(file_unit, *) BC(4) ! Right

    read(file_unit, *) ! Skip blank line
    read(file_unit, *) ! Skip section title
    read(file_unit, *) optimize_omega
    read(file_unit, *) benchmark_mode
    read(file_unit, *) quarter_domain

    read(file_unit, *) ! Skip blank line
    read(file_unit, *) ! Skip section title
    if (.NOT. (optimize_omega)) then
        read(file_unit, *) best_omega_psor
        read(file_unit, *) best_omega_lsor
        read(file_unit, *) best_omega_adir
    else
        read(file_unit, *) ! Skip PSOR Omega
        read(file_unit, *) ! Skip LSOR Omega
        read(file_unit, *) ! Skip ADIR Omega
    end if

    close(file_unit)

    ! Now that we know imax and jmax, we can create the arrays in memory
    allocate(T(imax, jmax))

    ! Setup for Quarter Domain (Problem 2b) dynamically based on full grid
    if (quarter_domain) then
        L = L / 2.0_wp
        W = W / 2.0_wp
    end if
    imax_q = (imax + 1) / 2
    jmax_q = (jmax + 1) / 2
    allocate(T_q(imax_q, jmax_q))

    ! Pre-compute grid spacing
    dx = L / real(imax - 1, wp)
    dy = W / real(jmax - 1, wp)

    ! INITIALIZATION & BOUNDARY CONDITIONS
    call reset_grid(T, imax, jmax, BC)
    
    print *, "--- 2D Heat Conduction Solver Initialized ---"
    print *, "Grid: ", imax, "x", jmax
    print *, "dx = ", dx, "m, dy = ", dy, "m"
    print *, "---------------------------------------------"

    ! RELAXATION FACTOR OPTIMIZERS (Sweeping omega 1.0 to 1.9)
    
    if (optimize_omega .and. .not. quarter_domain) then
         print *, "Optimizing relaxation factors for PSOR and LSOR..."
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
        ! between 1.25 and 1.3 to save time. If I had to do the entire range, I would just change the limits of this loop.
        do w_int = 1250, 1300
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

        ! Optimize LSOR
        min_iters_adir = 999999
        best_omega_adir = 1.0_wp

        ! Since I already ran between 1 to 1.9 with a step of 0.001,
        ! I already know the optimal omega is around 1.832, so I will just sweep,
        ! between 1.29 and 1.33 to save time. If I had to do the entire range, I would just change the limits of this loop.
        do w_int = 1290, 1330
            omega = real(w_int, wp) / 1000.0_wp
            call reset_grid(T, imax, jmax, BC)
            call solve_ADIR(T, imax, jmax, dx, dy, omega, iters, c_time)

            write(*, '(A10, I20, F20.6, F15.3)') "ADIR", iters, c_time, omega

            if (iters < min_iters_adir) then
                min_iters_adir = iters
                best_omega_adir = omega
            end if
        end do
        print *, ">> OPTIMAL ADIR: Omega = ", best_omega_adir, min_iters_adir, " iters"
        
    end if

    if (quarter_domain .and. optimize_omega) then
        print *, "Optimizing relaxation factor for LSOR with Quarter Domain Symmetry..."
        min_iters_lsor = 999999
        best_omega_lsor = 1.0_wp

        do w_int = 1250, 1300
            omega = real(w_int, wp) / 1000.0_wp
            T_q = 0.0_wp
            T_q(:, 1) = BC(1)     
            T_q(1, :) = BC(3)     
            T_q(1, 1) = (BC(1) + BC(3)) / 2.0_wp  
            call solve_LSOR_sym(T_q, imax_q, jmax_q, dx, dy, omega, iters, c_time)

            write(*, '(A10, I20, F20.6, F15.3)') "LSOR_Q", iters, c_time, omega

            if (iters < min_iters_lsor) then
                min_iters_lsor = iters
                best_omega_lsor = omega
            end if
        end do
        print *, ">> OPTIMAL LSOR with Quarter Symmetry: Omega = ", best_omega_lsor, min_iters_lsor, " iters"
    end if

    ! DATA EXPORT
    ! We will use LSOR to generate the final steady-state field since all 
    ! methods converge to the exact same physical answer.
    if (quarter_domain) then
        print *, "Running Quarter Domain Solver for: ", trim(base_name)
        ! Problem 2 Case B (Quarter Domain Symmetry)
        T_q = 0.0_wp
    
        ! Apply Dirichlet Boundaries
        T_q(:, 1) = BC(1)     ! Bottom Wall
        T_q(1, :) = BC(3)     ! Left Wall
        T_q(1, 1) = (BC(1) + BC(3)) / 2.0_wp  ! Corner Averaging

        call solve_LSOR_sym(T_q, imax_q, jmax_q, dx, dy, best_omega_lsor, iters, c_time)
        call export_to_csv(trim(out_file), T_q, imax_q, jmax_q, dx, dy)
    else
        print *, "Running Full Domain Solver for: ", trim(base_name)
        call reset_grid(T, imax, jmax, BC)
        call solve_LSOR(T, imax, jmax, dx, dy, best_omega_lsor, iters, c_time)
        call export_to_csv(trim(out_file), T, imax, jmax, dx, dy)
    end if
    
    ! TIMING ANALYSIS 
    ! For a more rigorous timing analysis, we could repeat the solver multiple times and average the time.
    ! This is especially useful for very fast methods where timing can be noisy.
    if (benchmark_mode) then
        open(newunit=csv_id, file='results/performance_' // trim(base_name) // ".csv", status='replace')
        write(csv_id, '(A20, A, A20, A, A20, A, A15)') "Method", ",", "Iterations", ",", "Comp. Time (ms)", ",", "Omega"
        close(csv_id)

        print *, "Running benchmarks for all solvers..."
        if (quarter_domain) then
            call benchmark_solver("LSOR_P2b", trim(base_name), best_omega_lsor)
        else
            call benchmark_solver("PGS", trim(base_name))
            call benchmark_solver("LGS", trim(base_name))
            call benchmark_solver("ADI", trim(base_name))
            call benchmark_solver("PSOR", trim(base_name), best_omega_psor)
            call benchmark_solver("LSOR", trim(base_name), best_omega_lsor)
            call benchmark_solver("ADIR", trim(base_name), best_omega_adir)
        end if
    end if    

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

        ! Corners (Average of adjacent walls)
        T_array(1, 1)         = (BC_list(1) +  BC_list(3)) / 2.0_wp  ! Bottom-Left
        T_array(max_i, 1)     = (BC_list(1) +  BC_list(4)) / 2.0_wp  ! Bottom-Right
        T_array(1, max_j)     = (BC_list(2) +  BC_list(3)) / 2.0_wp  ! Top-Left
        T_array(max_i, max_j) = (BC_list(2) +  BC_list(4)) / 2.0_wp  ! Top-Right
    end subroutine reset_grid

    ! Writes the 2D grid into a 1D X, Y, T format for Python
    subroutine export_to_csv(outfilename, T_array, max_i, max_j, d_x, d_y)
        character(len=*), intent(in) :: outfilename
        real(wp), intent(in)         :: T_array(:,:), d_x, d_y
        integer, intent(in)          :: max_i, max_j
        
        integer  :: i, j
        real(wp) :: x_pos, y_pos

        open(newunit=csv_id, file=outfilename, status='replace')
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
        print *, ">> Steady-state data exported to: ", outfilename
    end subroutine export_to_csv

    ! Benchmarker
    subroutine benchmark_solver(solver_name, base_name_val, omega_val)
        character(len=*), intent(in)   :: solver_name
        real(wp), intent(in), optional :: omega_val
        character(len=*), intent(in)   :: base_name_val

        real(wp) :: avg_time, dummy_c_time
        integer  :: rep
        integer(kind=8) :: t_start, t_end, t_rate 

        if (present(omega_val)) then
            omega = omega_val
        else
            omega = 1.0_wp
        end if
        
        base_name = base_name_val
        avg_time = 0.0_wp
        
        ! Start high-resolution timer
        call system_clock(t_start, t_rate)
        
        select case (solver_name)
            
            case ("PGS")
                do rep = 1, N_REPEAT
                    call reset_grid(T, imax, jmax, BC)
                    call solve_PSOR(T, imax, jmax, dx, dy, omega, iters, dummy_c_time)
                end do
                
            case ("LGS")
                do rep = 1, N_REPEAT
                    call reset_grid(T, imax, jmax, BC)
                    call solve_LSOR(T, imax, jmax, dx, dy, omega, iters, dummy_c_time)
                end do
                
            case ("ADI")
                do rep = 1, N_REPEAT
                    call reset_grid(T, imax, jmax, BC)
                    call solve_ADIR(T, imax, jmax, dx, dy, omega, iters, dummy_c_time)
                end do
                
            case ("PSOR")
                do rep = 1, N_REPEAT
                    call reset_grid(T, imax, jmax, BC)
                    call solve_PSOR(T, imax, jmax, dx, dy, omega, iters, dummy_c_time)
                end do
                
            case ("LSOR")
                do rep = 1, N_REPEAT
                    call reset_grid(T, imax, jmax, BC)
                    call solve_LSOR(T, imax, jmax, dx, dy, omega, iters, dummy_c_time)
                end do
                
            case ("ADIR")
                do rep = 1, N_REPEAT
                    call reset_grid(T, imax, jmax, BC)
                    call solve_ADIR(T, imax, jmax, dx, dy, omega, iters, dummy_c_time)
                end do

            case ("LSOR_P2b")
                do rep = 1, N_REPEAT
                    T_q = 0.0_wp
                    T_q(:, 1) = BC(1)     
                    T_q(1, :) = BC(3)     
                    call solve_LSOR_sym(T_q, imax_q, jmax_q, dx, dy, omega, iters, dummy_c_time)
                end do
                
            case default
                print *, "ERROR: Unknown solver passed to benchmark routine!"
                stop
        end select
        
        ! End high-resolution timer
        call system_clock(t_end)
        
        ! Calculate accurate average time in milliseconds
        avg_time = real(t_end - t_start, wp) / real(t_rate, wp) / real(N_REPEAT, wp) * 1000.0_wp 

        open(newunit=csv_id, file='results/performance_' // trim(base_name) // ".csv", status='old', position='append', action='write')
        if (present(omega_val)) then
            write(csv_id, '(A20, A, I20, A, F20.6, A, F15.3)') solver_name // "_" // trim(base_name), ",", iters, ",", avg_time, ",", omega
        else
            write(csv_id, '(A20, A, I20, A, F20.6, A, A15)') solver_name // "_" // trim(base_name), ",", iters, ",", avg_time, ",", "N/A"
        end if
        close(csv_id)
        
        print *, "Benchmarked: ", solver_name, " | Avg Time: ", avg_time, "ms | Iters: ", iters
    end subroutine benchmark_solver
end program main