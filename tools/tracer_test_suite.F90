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

module tracer_test_suite_mod
    use fv_tracer2d_mod, only: tracer_2d, tracer_2d_1L
    use fv_tracer2d_optimized_mod, only: tracer_2d_optimized, tracer_2d_1L_optimized
    use tp_core_mod, only: fv_tp_2d
    use tp_core_optimized_mod, only: fv_tp_2d_optimized
    use performance_measurement_mod, only: performance_timer_start, performance_timer_stop, &
                                          get_elapsed_time, print_performance_summary
    use fv_arrays_mod, only: fv_grid_type, fv_grid_bounds_type
    use mpp_domains_mod, only: domain2d
    use mpp_mod, only: mpp_error, FATAL

    implicit none
    private

    public :: run_all_tests, run_accuracy_tests, run_performance_tests, run_stability_tests

contains

    ! Main test runner
    subroutine run_all_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(fv_grid_bounds_type), intent(IN) :: bd
        type(domain2d), intent(INOUT) :: domain
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac

        write(*,*) 'Running comprehensive tracer test suite...'
        write(*,*) '=========================================='

        call run_accuracy_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)
        call run_performance_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)
        call run_stability_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)

        write(*,*) '=========================================='
        write(*,*) 'Tracer test suite completed.'
    end subroutine run_all_tests

    ! Run accuracy tests
    subroutine run_accuracy_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(fv_grid_bounds_type), intent(IN) :: bd
        type(domain2d), intent(INOUT) :: domain
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac

        real, allocatable :: q_orig(:,:,:,:), q_orig_1L(:,:,:,:)
        real, allocatable :: q_test(:,:,:,:), q_test_1L(:,:,:,:)
        real, allocatable :: dp1(:,:,:), mfx(:,:,:), mfy(:,:,:)
        real, allocatable :: cx(:,:,:), cy(:,:,:)
        real, allocatable :: q_pack_data(:,:,:,:), dp1_pack_data(:,:,:,:)
        type(group_halo_update_type) :: q_pack, dp1_pack
        real :: error_norm, max_error, rel_error
        integer :: i, j, k, iq
        logical :: test_passed

        write(*,*) 'Running accuracy tests...'

        ! Allocate test arrays
        allocate(q_orig(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(q_orig_1L(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(q_test(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(q_test_1L(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz))
        allocate(mfx(bd%is:bd%ie+1,bd%js:bd%je,npz))
        allocate(mfy(bd%is:bd%ie,bd%js:bd%je+1,npz))
        allocate(cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz))
        allocate(cy(bd%isd:bd%ied,bd%js:bd%je+1,npz))
        allocate(q_pack_data(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(dp1_pack_data(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))

        ! Initialize test data with known values
        call initialize_test_data(q_orig, dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)

        ! Copy initial data to test arrays
        q_test = q_orig
        q_test_1L = q_orig
        q_orig_1L = q_orig

        ! Run original tracer_2d
        call performance_timer_start('original_tracer_2d')
        call tracer_2d(q_orig, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                      npx, npy, npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)
        call performance_timer_stop('original_tracer_2d')

        ! Run optimized tracer_2d
        call performance_timer_start('optimized_tracer_2d')
        call tracer_2d_optimized(q_test, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                npx, npy, npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)
        call performance_timer_stop('optimized_tracer_2d')

        ! Calculate error between original and optimized results
        error_norm = 0.0
        max_error = 0.0
        do iq = 1, nq
            do k = 1, npz
                do j = bd%js, bd%je
                    do i = bd%is, bd%ie
                        rel_error = abs(q_orig(i,j,k,iq) - q_test(i,j,k,iq))
                        error_norm = error_norm + rel_error**2
                        if (rel_error > max_error) max_error = rel_error
                    end do
                end do
            end do
        end do

        error_norm = sqrt(error_norm / (nq * npz * (bd%ie-bd%is+1) * (bd%je-bd%js+1)))

        ! Check if results are within acceptable tolerance
        test_passed = (max_error < 1.0e-10) .and. (error_norm < 1.0e-12)

        write(*,*) 'Accuracy Test Results:'
        write(*,'(A, E15.7)') '  Max absolute error: ', max_error
        write(*,'(A, E15.7)') '  RMS error: ', error_norm
        write(*,'(A, L1)') '  Test passed: ', test_passed

        ! Run 1L version tests
        call performance_timer_start('original_tracer_2d_1L')
        call tracer_2d_1L(q_orig_1L, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                         npx, npy, npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)
        call performance_timer_stop('original_tracer_2d_1L')

        call performance_timer_start('optimized_tracer_2d_1L')
        call tracer_2d_1L_optimized(q_test_1L, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                   npx, npy, npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)
        call performance_timer_stop('optimized_tracer_2d_1L')

        ! Calculate error for 1L version
        error_norm = 0.0
        max_error = 0.0
        do iq = 1, nq
            do k = 1, npz
                do j = bd%js, bd%je
                    do i = bd%is, bd%ie
                        rel_error = abs(q_orig_1L(i,j,k,iq) - q_test_1L(i,j,k,iq))
                        error_norm = error_norm + rel_error**2
                        if (rel_error > max_error) max_error = rel_error
                    end do
                end do
            end do
        end do

        error_norm = sqrt(error_norm / (nq * npz * (bd%ie-bd%is+1) * (bd%je-bd%js+1)))
        test_passed = test_passed .and. ((max_error < 1.0e-10) .and. (error_norm < 1.0e-12))

        write(*,*) '1L Accuracy Test Results:'
        write(*,'(A, E15.7)') '  Max absolute error: ', max_error
        write(*,'(A, E15.7)') ' RMS error: ', error_norm
        write(*,'(A, L1)') '  Test passed: ', (max_error < 1.0e-10) .and. (error_norm < 1.0e-12)

        write(*,'(A, L1)') 'Overall accuracy test passed: ', test_passed

        ! Deallocate arrays
        deallocate(q_orig, q_orig_1L, q_test, q_test_1L, dp1, mfx, mfy, cx, cy, q_pack_data, dp1_pack_data)
    end subroutine run_accuracy_tests

    ! Run performance tests
    subroutine run_performance_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(fv_grid_bounds_type), intent(IN) :: bd
        type(domain2d), intent(INOUT) :: domain
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac

        real, allocatable :: q(:,:,:,:), dp1(:,:,:)
        real, allocatable :: mfx(:,:,:), mfy(:,:,:)
        real, allocatable :: cx(:,:,:), cy(:,:,:)
        real, allocatable :: q_pack_data(:,:,:,:), dp1_pack_data(:,:,:,:)
        type(group_halo_update_type) :: q_pack, dp1_pack
        real :: time_orig, time_opt, speedup
        integer :: i, test_iterations = 5  ! Reduced for faster testing

        write(*,*) 'Running performance tests...'

        ! Allocate test arrays
        allocate(q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz))
        allocate(mfx(bd%is:bd%ie+1,bd%js:bd%je,npz))
        allocate(mfy(bd%is:bd%ie,bd%js:bd%je+1,npz))
        allocate(cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz))
        allocate(cy(bd%isd:bd%ied,bd%js:bd%je+1,npz))
        allocate(q_pack_data(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(dp1_pack_data(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))

        ! Initialize test data
        call initialize_test_data(q, dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)

        ! Performance test for tracer_2d
        call performance_timer_start('perf_original_tracer_2d')
        do i = 1, test_iterations
            call tracer_2d(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                          npx, npy, npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)
        end do
        call performance_timer_stop('perf_original_tracer_2d')

        time_orig = get_elapsed_time('perf_original_tracer_2d')

        ! Reset test data
        call initialize_test_data(q, dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)

        call performance_timer_start('perf_optimized_tracer_2d')
        do i = 1, test_iterations
            call tracer_2d_optimized(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                    npx, npy, npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)
        end do
        call performance_timer_stop('perf_optimized_tracer_2d')

        time_opt = get_elapsed_time('perf_optimized_tracer_2d')

        ! Calculate performance improvement
        if (time_opt > 0.0) then
            speedup = (time_orig / test_iterations) / (time_opt / test_iterations)
        else
            speedup = 0.0
        endif

        write(*,*) 'Performance Test Results:'
        write(*,'(A, F10.6, A)') '  Original time per iteration: ', time_orig/test_iterations, ' seconds'
        write(*,'(A, F10.6, A)') '  Optimized time per iteration: ', time_opt/test_iterations, ' seconds'
        write(*,'(A, F10.2, A)') '  Speedup: ', speedup, 'x'
        write(*,'(A, F10.2, A)') '  Performance improvement: ', (1.0 - (time_opt/test_iterations)/(time_orig/test_iterations))*100.0, '%'

        ! Deallocate arrays
        deallocate(q, dp1, mfx, mfy, cx, cy, q_pack_data, dp1_pack_data)

        call print_performance_summary()
    end subroutine run_performance_tests

    ! Run stability tests
    subroutine run_stability_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(fv_grid_bounds_type), intent(IN) :: bd
        type(domain2d), intent(INOUT) :: domain
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac

        real, allocatable :: q(:,:,:,:), dp1(:,:,:)
        real, allocatable :: mfx(:,:,:), mfy(:,:,:)
        real, allocatable :: cx(:,:,:), cy(:,:,:)
        real, allocatable :: q_pack_data(:,:,:,:), dp1_pack_data(:,:,:,:)
        type(group_halo_update_type) :: q_pack, dp1_pack
        real :: min_val, max_val, total_mass, initial_mass
        integer :: i, j, k, iq, iter, max_iter = 10
        logical :: stability_ok = .true.

        write(*,*) 'Running stability tests...'

        ! Allocate test arrays
        allocate(q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz))
        allocate(mfx(bd%is:bd%ie+1,bd%js:bd%je,npz))
        allocate(mfy(bd%is:bd%ie,bd%js:bd%je+1,npz))
        allocate(cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz))
        allocate(cy(bd%isd:bd%ied,bd%js:bd%je+1,npz))
        allocate(q_pack_data(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(dp1_pack_data(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))

        ! Initialize test data
        call initialize_test_data(q, dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)

        ! Calculate initial mass
        initial_mass = 0.0
        do iq = 1, nq
            do k = 1, npz
                do j = bd%js, bd%je
                    do i = bd%is, bd%ie
                        initial_mass = initial_mass + q(i,j,k,iq) * dp1(i,j,k)
                    end do
                end do
            end do
        end do

        ! Run multiple iterations to test stability
        do iter = 1, max_iter
            ! Run optimized tracer (more efficient for stability testing)
            call tracer_2d_optimized(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                    npx, npy, npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)

            ! Check for stability (finite values, no NaN or Inf)
            min_val = huge(1.0)
            max_val = -huge(1.0)
            total_mass = 0.0

            do iq = 1, nq
                do k = 1, npz
                    do j = bd%js, bd%je
                        do i = bd%is, bd%ie
                            if (q(i,j,k,iq) /= q(i,j,k,iq)) then  ! Check for NaN
                                write(*,'(A, I0, A, I0, A, I0, A, I0)') &
                                    'Stability test failed: NaN detected at iteration ', iter, &
                                    ', tracer ', iq, ', level ', k, ', grid point (', i, ',', j, ')'
                                stability_ok = .false.
                                goto 999
                            endif
                            if (abs(q(i,j,k,iq)) > huge(1.0)/2.0) then  ! Check for Inf
                                write(*,'(A, I0, A, I0, A, I0, A, I0)') &
                                    'Stability test failed: Inf detected at iteration ', iter, &
                                    ', tracer ', iq, ', level ', k, ', grid point (', i, ',', j, ')'
                                stability_ok = .false.
                                goto 999
                            endif

                            if (q(i,j,k,iq) < min_val) min_val = q(i,j,k,iq)
                            if (q(i,j,k,iq) > max_val) max_val = q(i,j,k,iq)
                            total_mass = total_mass + q(i,j,k,iq) * dp1(i,j,k)
                        end do
                    end do
                end do
            end do

            ! Check mass conservation
            if (abs(total_mass - initial_mass) / initial_mass > 1.0e-6) then
                write(*,'(A, I0, A, E15.7, A, E15.7, A, E15.7)') &
                    'Mass conservation warning at iteration ', iter, &
                    ': Initial = ', initial_mass, ', Current = ', total_mass, &
                    ', Difference = ', abs(total_mass - initial_mass)
            endif

            write(*,'(A, I0, A, E15.7, A, E15.7)') &
                '  Iteration ', iter, ': Min = ', min_val, ', Max = ', max_val
        end do

999     continue

        write(*,*) 'Stability Test Results:'
        write(*,'(A, L1)') ' Stability maintained: ', stability_ok
        if (stability_ok) then
            write(*,*) '  All values remained finite throughout testing'
        endif

        ! Deallocate arrays
        deallocate(q, dp1, mfx, mfy, cx, cy, q_pack_data, dp1_pack_data)
    end subroutine run_stability_tests

    ! Initialize test data with known values
    subroutine initialize_test_data(q, dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)
        real, intent(out) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)
        real, intent(out) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
        real, intent(out) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,npz)
        real, intent(out) :: mfy(bd%is:bd%ie,bd%js:bd%je+1,npz)
        real, intent(out) :: cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz)
        real, intent(out) :: cy(bd%isd:bd%ied,bd%js:bd%je+1,npz)
        type(fv_grid_bounds_type), intent(IN) :: bd
        integer, intent(IN) :: npx, npy, npz, nq

        integer :: i, j, k, iq

        ! Initialize with simple patterns for testing
        do iq = 1, nq
            do k = 1, npz
                do j = bd%jsd, bd%jed
                    do i = bd%isd, bd%ied
                        q(i,j,k,iq) = sin(real(i)/npx * 3.14159) * cos(real(j)/npy * 3.14159) * real(k)/npz
                    end do
                end do
            end do
        end do

        do k = 1, npz
            do j = bd%jsd, bd%jed
                do i = bd%isd, bd%ied
                    dp1(i,j,k) = 1.0 + 0.1 * sin(real(i)/npx * 6.28318) * cos(real(j)/npy * 6.28318) * real(k)/npz
                end do
            end do
        end do

        do k = 1, npz
            do j = bd%js, bd%je
                do i = bd%is, bd%ie+1
                    mfx(i,j,k) = 0.1 * sin(real(i)/npx * 6.28318) * cos(real(j)/npy * 3.14159) * real(k)/npz
                    cx(i,j,k) = 0.2 * sin(real(i)/npx * 3.14159) * cos(real(j)/npy * 6.28318) * real(k)/npz
                end do
            end do
        end do

        do k = 1, npz
            do j = bd%js, bd%je+1
                do i = bd%is, bd%ie
                    mfy(i,j,k) = 0.1 * cos(real(i)/npx * 3.14159) * sin(real(j)/npy * 6.28318) * real(k)/npz
                    cy(i,j,k) = 0.2 * cos(real(i)/npx * 6.28318) * sin(real(j)/npy * 3.14159) * real(k)/npz
                end do
            end do
        end do
    end subroutine initialize_test_data

end module tracer_test_suite_mod