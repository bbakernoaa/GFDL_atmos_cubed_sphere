!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module performance_measurement_mod
    use mpp_mod,           only: mpp_error, FATAL, mpp_broadcast, mpp_send, mpp_recv, mpp_sum, mpp_max
    use mpp_domains_mod,   only: domain2d
    use fv_arrays_mod,     only: fv_grid_type, fv_grid_bounds_type
    use fv_timing_mod,     only: timing_on, timing_off

    implicit none
    private

    public :: performance_timer, tracer_performance_test, memory_usage_monitor

    ! Performance timer structure
    type performance_timer_type
        real :: start_time
        real :: end_time
        real :: elapsed_time
        character(len=64) :: name
        integer :: call_count
    end type performance_timer_type

    type(performance_timer_type), allocatable, save :: performance_timers(:)
    integer, save :: num_timers = 0
    integer, parameter :: max_timers = 100
    real, parameter :: PERFORMANCE_EPSILON = 1.0e-12  ! Small value to avoid division by zero

contains

    ! Initialize performance timer array
    subroutine init_performance_timers()
        if (allocated(performance_timers)) then
            deallocate(performance_timers)
        endif
        allocate(performance_timers(max_timers))
        num_timers = 0
    end subroutine init_performance_timers

    ! Start a performance timer
    subroutine performance_timer_start(timer_name)
        character(len=*), intent(in) :: timer_name
        integer :: i, timer_idx

        ! Check if timer already exists
        timer_idx = -1
        do i = 1, num_timers
            if (trim(performance_timers(i)%name) == trim(timer_name)) then
                timer_idx = i
                exit
            endif
        end do

        ! If timer doesn't exist, create a new one
        if (timer_idx == -1) then
            if (num_timers >= max_timers) then
                call mpp_error(FATAL, 'performance_measurement_mod: Maximum number of timers exceeded')
            endif
            num_timers = num_timers + 1
            timer_idx = num_timers
            performance_timers(timer_idx)%name = trim(timer_name)
            performance_timers(timer_idx)%call_count = 0
        endif

        ! Start timing
        call cpu_time(performance_timers(timer_idx)%start_time)
        performance_timers(timer_idx)%call_count = performance_timers(timer_idx)%call_count + 1
    end subroutine performance_timer_start

    ! Stop a performance timer and record elapsed time
    subroutine performance_timer_stop(timer_name)
        character(len=*), intent(in) :: timer_name
        integer :: i, timer_idx

        ! Find the timer
        timer_idx = -1
        do i = 1, num_timers
            if (trim(performance_timers(i)%name) == trim(timer_name)) then
                timer_idx = i
                exit
            endif
        end do

        if (timer_idx == -1) then
            call mpp_error(FATAL, 'performance_measurement_mod: Timer not found: ' // trim(timer_name))
        endif

        ! Stop timing and record elapsed time
        call cpu_time(performance_timers(timer_idx)%end_time)
        performance_timers(timer_idx)%elapsed_time = &
            performance_timers(timer_idx)%end_time - performance_timers(timer_idx)%start_time
    end subroutine performance_timer_stop

    ! Get elapsed time for a timer
    function get_elapsed_time(timer_name) result(elapsed_time)
        character(len=*), intent(in) :: timer_name
        real :: elapsed_time
        integer :: i, timer_idx

        elapsed_time = -1.0  ! Initialize to -1 to indicate not found
        do i = 1, num_timers
            if (trim(performance_timers(i)%name) == trim(timer_name)) then
                elapsed_time = performance_timers(i)%elapsed_time
                exit
            endif
        end do
    end function get_elapsed_time

    ! Print performance summary
    subroutine print_performance_summary()
        integer :: i

        write(*,*) '=== Performance Measurement Summary ==='
        write(*,*) 'Timer Name                    | Elapsed Time (s) | Calls | Avg Time (s)'
        write(*,*) '------------------------------|------------------|-------|-------------'
        do i = 1, num_timers
            write(*,'(A30, F18.6, I8, F12.6)') &
                trim(performance_timers(i)%name), &
                performance_timers(i)%elapsed_time, &
                performance_timers(i)%call_count, &
                performance_timers(i)%elapsed_time / max(1, performance_timers(i)%call_count)
        end do
        write(*,*) '========================================'
    end subroutine print_performance_summary

    ! Memory usage monitoring
    subroutine memory_usage_monitor(current_usage, peak_usage, label)
        integer, intent(inout) :: current_usage, peak_usage
        character(len=*), intent(in) :: label
        integer :: temp_usage

        ! In a real implementation, this would interface with system calls
        ! to get actual memory usage. For now, we'll just simulate it.
        temp_usage = current_usage
        if (temp_usage > peak_usage) then
            peak_usage = temp_usage
        endif

        write(*,'(A, A, A, I10, A, I10)') &
            'Memory Usage - ', trim(label), &
            ': Current = ', current_usage, &
            ' bytes, Peak = ', peak_usage, ' bytes'
    end subroutine memory_usage_monitor

    ! Performance test for tracer transport
    subroutine tracer_performance_test(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                      npx, npy, npz, nq, hord, dt, lim_fac, test_iterations)
        type(fv_grid_bounds_type), intent(IN) :: bd
        integer, intent(IN) :: npx, npy, npz, nq, hord, test_iterations
        real, intent(IN) :: dt, lim_fac
        real, intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)
        real, intent(INOUT) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
        real, intent(INOUT) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,npz)
        real, intent(INOUT) :: mfy(bd%is:bd%ie,bd%js:bd%je+1,npz)
        real, intent(INOUT) :: cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz)
        real, intent(INOUT) :: cy(bd%isd:bd%ied,bd%js:bd%je+1,npz)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(domain2d), intent(INOUT) :: domain

        ! Local variables for timing
        integer :: i
        real :: total_time_old, total_time_new
        real :: avg_time_old, avg_time_new
        real :: speedup_factor

        write(*,*) 'Starting tracer performance test...'
        write(*,'(A, I0, A)') 'Running ', test_iterations, ' iterations'

        ! Test the original tracer implementation (simulated)
        call performance_timer_start('tracer_original')
        do i = 1, test_iterations
            ! Simulate original tracer computation
            call simulate_original_tracer_computation(q, dp1, mfx, mfy, cx, cy, &
                                                    gridstruct, bd, domain, &
                                                    npx, npy, npz, nq, hord, dt, lim_fac)
        end do
        call performance_timer_stop('tracer_original')
        total_time_old = get_elapsed_time('tracer_original')
        avg_time_old = total_time_old / test_iterations

        ! Test the optimized tracer implementation
        call performance_timer_start('tracer_optimized')
        do i = 1, test_iterations
            ! Simulate optimized tracer computation
            call simulate_optimized_tracer_computation(q, dp1, mfx, mfy, cx, cy, &
                                                    gridstruct, bd, domain, &
                                                    npx, npy, npz, nq, hord, dt, lim_fac)
        end do
        call performance_timer_stop('tracer_optimized')
        total_time_new = get_elapsed_time('tracer_optimized')
        avg_time_new = total_time_new / test_iterations

        ! Calculate performance improvement
        if (avg_time_new > 0.0) then
            speedup_factor = avg_time_old / avg_time_new
        else
            speedup_factor = 0.0
        endif

        ! Print results
        write(*,*) '=== Tracer Performance Test Results ==='
        write(*,'(A, F10.6, A)') 'Original average time per iteration: ', avg_time_old, ' seconds'
        write(*,'(A, F10.6, A)') 'Optimized average time per iteration: ', avg_time_new, ' seconds'
        write(*,'(A, F10.2, A)') 'Speedup factor: ', speedup_factor, 'x'
        write(*,'(A, F10.2, A)') 'Performance improvement: ', (1.0 - avg_time_new/avg_time_old)*100.0, '%'
        write(*,*) '======================================='

        call print_performance_summary()
    end subroutine tracer_performance_test

    ! Simulate original tracer computation (placeholder)
    subroutine simulate_original_tracer_computation(q, dp1, mfx, mfy, cx, cy, &
                                                   gridstruct, bd, domain, &
                                                   npx, npy, npz, nq, hord, dt, lim_fac)
        type(fv_grid_bounds_type), intent(IN) :: bd
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac
        real, intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)
        real, intent(INOUT) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
        real, intent(INOUT) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,npz)
        real, intent(INOUT) :: mfy(bd%is:bd%ie,bd%js:bd%je+1,npz)
        real, intent(INOUT) :: cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz)
        real, intent(INOUT) :: cy(bd%isd:bd%ied,bd%js:bd%je+1,npz)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(domain2d), intent(INOUT) :: domain

        ! Placeholder implementation - in real code this would call the original tracer routine
        ! For simulation purposes, we'll just perform some dummy operations
        integer :: i, j, k, iq

        do iq = 1, nq
            do k = 1, npz
                do j = bd%js, bd%je
                    do i = bd%is, bd%ie
                        ! Simulate some computation
                        q(i,j,k,iq) = q(i,j,k,iq) * 0.999 + dp1(i,j,k) * 0.001
                    end do
                end do
            end do
        end do
    end subroutine simulate_original_tracer_computation

    ! Simulate optimized tracer computation (placeholder)
    subroutine simulate_optimized_tracer_computation(q, dp1, mfx, mfy, cx, cy, &
                                                   gridstruct, bd, domain, &
                                                   npx, npy, npz, nq, hord, dt, lim_fac)
        type(fv_grid_bounds_type), intent(IN) :: bd
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac
        real, intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)
        real, intent(INOUT) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
        real, intent(INOUT) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,npz)
        real, intent(INOUT) :: mfy(bd%is:bd%ie,bd%js:bd%je+1,npz)
        real, intent(INOUT) :: cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz)
        real, intent(INOUT) :: cy(bd%isd:bd%ied,bd%js:bd%je+1,npz)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(domain2d), intent(INOUT) :: domain

        ! Placeholder implementation - in real code this would call the optimized tracer routine
        ! For simulation purposes, we'll just perform some dummy operations
        integer :: i, j, k, iq

        do iq = 1, nq
            do k = 1, npz
                do j = bd%js, bd%je
                    do i = bd%is, bd%ie
                        ! Simulate some computation with optimized approach
                        q(i,j,k,iq) = q(i,j,k,iq) * 0.9995 + dp1(i,j,k) * 0.0005
                    end do
                end do
            end do
        end do
    end subroutine simulate_optimized_tracer_computation

end module performance_measurement_mod