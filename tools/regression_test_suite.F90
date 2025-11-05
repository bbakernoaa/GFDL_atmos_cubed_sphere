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

module regression_test_suite_mod
    use fv_tracer2d_mod, only: tracer_2d
    use fv_tracer2d_optimized_mod, only: tracer_2d_optimized
    use performance_measurement_mod, only: performance_timer_start, performance_timer_stop, &
                                          get_elapsed_time, print_performance_summary
    use benchmark_test_suite_mod, only: run_comprehensive_benchmark_tests
    use fv_arrays_mod, only: fv_grid_type, fv_grid_bounds_type
    use mpp_domains_mod, only: domain2d
    use mpp_mod, only: mpp_error, FATAL, mpp_sum, mpp_max

    implicit none
    private

    public :: run_regression_tests, run_performance_regression_tests

    ! Performance baseline values (these would be updated with actual measured values)
    real, parameter :: BASELINE_ORIG_TIME = 1.0      ! Placeholder baseline for original implementation
    real, parameter :: BASELINE_OPT_TIME = 0.7      ! Placeholder baseline for optimized implementation
    real, parameter :: BASELINE_SPEEDUP = 1.43     ! Expected speedup (orig/opt)
    real, parameter :: BASELINE_ACCURACY_ERROR = 1.0e-13  ! Maximum acceptable error
    real, parameter :: BASELINE_MASS_CONSERVATION_ERROR = 1.0e-12  ! Mass conservation tolerance

    ! Thresholds for regression detection
    real, parameter :: PERFORMANCE_DEGRADATION_THRESHOLD = 0.95  ! Performance should not drop below 95% of baseline
    real, parameter :: ACCURACY_DEGRADATION_THRESHOLD = 1.5e-13   ! Accuracy should not worsen beyond this
    
    ! Configuration parameters
    integer, parameter :: REGRESSION_TEST_ITERATIONS = 5

contains

    ! Main regression test runner
    subroutine run_regression_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(fv_grid_bounds_type), intent(IN) :: bd
        type(domain2d), intent(INOUT) :: domain
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac

        logical :: all_tests_passed = .true.

        write(*,*) 'Running regression tests...'
        write(*,*) '============================'

        ! Run performance regression tests
        if (.not. run_performance_regression_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)) then
            all_tests_passed = .false.
            write(*,*) 'Performance regression tests FAILED'
        else
            write(*,*) 'Performance regression tests PASSED'
        endif

        ! Run accuracy regression tests
        if (.not. run_accuracy_regression_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)) then
            all_tests_passed = .false.
            write(*,*) 'Accuracy regression tests FAILED'
        else
            write(*,*) 'Accuracy regression tests PASSED'
        endif

        ! Run stability regression tests
        if (.not. run_stability_regression_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)) then
            all_tests_passed = .false.
            write(*,*) 'Stability regression tests FAILED'
        else
            write(*,*) 'Stability regression tests PASSED'
        endif

        write(*,*) '============================'
        write(*,'(A, L1)') 'Overall regression tests passed: ', all_tests_passed

        if (.not. all_tests_passed) then
            call mpp_error(FATAL, 'Regression tests failed - performance or accuracy degradation detected')
        endif
    end subroutine run_regression_tests

    ! Run performance regression tests
    function run_performance_regression_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac) result(test_passed)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(fv_grid_bounds_type), intent(IN) :: bd
        type(domain2d), intent(INOUT) :: domain
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac
        logical :: test_passed

        real, allocatable :: q(:,:,:,:), dp1(:,:,:)
        real, allocatable :: mfx(:,:,:), mfy(:,:,:)
        real, allocatable :: cx(:,:,:), cy(:,:,:)
        type(group_halo_update_type) :: q_pack, dp1_pack
        real :: time_orig, time_opt, speedup
        integer :: i, test_iterations = 5
        real :: measured_speedup

        test_passed = .true.

        write(*,*) 'Running performance regression tests...'

        ! Allocate test arrays
        allocate(q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz))
        allocate(mfx(bd%is:bd%ie+1,bd%js:bd%je,npz))
        allocate(mfy(bd%is:bd%ie,bd%js:bd%je+1,npz))
        allocate(cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz))
        allocate(cy(bd%isd:bd%ied,bd%js:bd%je+1,npz))

        ! Initialize test data
        call initialize_test_data(q, dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)

        ! Test original implementation performance
        call performance_timer_start('regression_original')
        do i = 1, test_iterations
            call tracer_2d(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                          npx, npy, npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)
            ! Reset data between iterations to ensure consistent test conditions
            call initialize_test_data(q, dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)
        end do
        call performance_timer_stop('regression_original')
        time_orig = get_elapsed_time('regression_original')

        ! Test optimized implementation performance
        call performance_timer_start('regression_optimized')
        do i = 1, test_iterations
            call tracer_2d_optimized(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                    npx, npy, npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)
            ! Reset data between iterations to ensure consistent test conditions
            call initialize_test_data(q, dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)
        end do
        call performance_timer_stop('regression_optimized')
        time_opt = get_elapsed_time('regression_optimized')

        ! Calculate performance metrics
        if (time_opt > 0.0) then
            speedup = (time_orig / test_iterations) / (time_opt / test_iterations)
        else
            speedup = 0.0
        endif

        measured_speedup = speedup

        write(*,'(A, F10.6, A, F10.6, A, F10.2)') &
            ' Performance - Original: ', time_orig/test_iterations, &
            's, Optimized: ', time_opt/test_iterations, &
            's, Speedup: ', measured_speedup

        ! Check if performance meets baseline expectations
        if (measured_speedup < BASELINE_SPEEDUP * PERFORMANCE_DEGRADATION_THRESHOLD) then
            test_passed = .false.
            write(*,'(A, F10.2, A, F10.2)') &
                '  PERFORMANCE REGRESSION: Expected speedup >= ', BASELINE_SPEEDUP * PERFORMANCE_DEGRADATION_THRESHOLD, &
                ', but got: ', measured_speedup
        else
            write(*,'(A, F10.2, A, F10.2)') &
                ' Performance meets threshold: Expected >= ', BASELINE_SPEEDUP * PERFORMANCE_DEGRADATION_THRESHOLD, &
                ', Got: ', measured_speedup
        endif

        ! Additional scaling tests
        if (.not. run_scaling_regression_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)) then
            test_passed = .false.
        endif

        deallocate(q, dp1, mfx, mfy, cx, cy)
    end function run_performance_regression_tests

    ! Run scaling regression tests
    function run_scaling_regression_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac) result(test_passed)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(fv_grid_bounds_type), intent(IN) :: bd
        type(domain2d), intent(INOUT) :: domain
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac
        logical :: test_passed

        integer, parameter :: num_configs = 3
        integer :: configs_npx(num_configs) = [10, 20, 30]
        integer :: configs_npy(num_configs) = [10, 20, 30]
        integer :: configs_npz(num_configs) = [5, 10, 15]
        real :: speedups(num_configs)
        integer :: i
        real, allocatable :: q(:,:,:,:), dp1(:,:,:)
        real, allocatable :: mfx(:,:,:), mfy(:,:,:)
        real, allocatable :: cx(:,:,:), cy(:,:,:)
        type(group_halo_update_type) :: q_pack, dp1_pack
        real :: time_orig, time_opt
        integer :: current_npx, current_npy, current_npz

        test_passed = .true.
        write(*,*) '  Running scaling regression tests...'

        do i = 1, num_configs
            current_npx = configs_npx(i)
            current_npy = configs_npy(i)
            current_npz = configs_npz(i)

            ! Adjust bounds based on current config
            if (current_npx > npx .or. current_npy > npy .or. current_npz > npz) then
                cycle  ! Skip if current config exceeds available grid size
            endif

            allocate(q(bd%isd:bd%ied,bd%jsd:bd%jed,current_npz,nq))
            allocate(dp1(bd%isd:bd%ied,bd%jsd:bd%jed,current_npz))
            allocate(mfx(bd%is:bd%ie+1,bd%js:bd%je,current_npz))
            allocate(mfy(bd%is:bd%ie,bd%js:bd%je+1,current_npz))
            allocate(cx(bd%is:bd%ie+1,bd%jsd:bd%jed,current_npz))
            allocate(cy(bd%isd:bd%ied,bd%js:bd%je+1,current_npz))

            ! Initialize test data
            call initialize_test_data(q, dp1, mfx, mfy, cx, cy, bd, current_npx, current_npy, current_npz, nq)

            ! Test original implementation
            call performance_timer_start('scaling_orig_' // trim(str(i)))
            call tracer_2d(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                          current_npx, current_npy, current_npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)
            call performance_timer_stop('scaling_orig_' // trim(str(i)))
            time_orig = get_elapsed_time('scaling_orig_' // trim(str(i)))

            ! Reset data
            call initialize_test_data(q, dp1, mfx, mfy, cx, cy, bd, current_npx, current_npy, current_npz, nq)

            ! Test optimized implementation
            call performance_timer_start('scaling_opt_' // trim(str(i)))
            call tracer_2d_optimized(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                    current_npx, current_npy, current_npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)
            call performance_timer_stop('scaling_opt_' // trim(str(i)))
            time_opt = get_elapsed_time('scaling_opt_' // trim(str(i)))

            ! Calculate speedup
            if (time_opt > 0.0) then
                speedups(i) = time_orig / time_opt
            else
                speedups(i) = 0.0
            endif

            write(*,'(A, I0, A, I0, A, I0, A, F10.2, A)') &
                '    Config (', current_npx, 'x', current_npy, 'x', current_npz, &
                '): Speedup = ', speedups(i), 'x'

            ! Check if scaling performance is acceptable
            if (speedups(i) < 1.0) then  ! Should at least maintain performance
                test_passed = .false.
                write(*,'(A, I0, A, F10.2)') &
                    '    SCALING REGRESSION at config ', i, ': Speedup = ', speedups(i)
            endif

            deallocate(q, dp1, mfx, mfy, cx, cy)
        end do

        if (test_passed) then
            write(*,*) '    Scaling tests PASSED'
        else
            write(*,*) '    Scaling tests FAILED'
        endif
    end function run_scaling_regression_tests

    ! Run accuracy regression tests
    function run_accuracy_regression_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac) result(test_passed)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(fv_grid_bounds_type), intent(IN) :: bd
        type(domain2d), intent(INOUT) :: domain
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac
        logical :: test_passed

        real, allocatable :: q_orig(:,:,:,:), q_opt(:,:,:,:)
        real, allocatable :: dp1(:,:,:)
        real, allocatable :: mfx(:,:,:), mfy(:,:,:)
        real, allocatable :: cx(:,:,:), cy(:,:,:)
        type(group_halo_update_type) :: q_pack, dp1_pack
        real :: error_norm, max_error, rel_error
        integer :: i, j, k, iq
        real :: measured_accuracy_error

        test_passed = .true.

        write(*,*) 'Running accuracy regression tests...'

        ! Allocate test arrays
        allocate(q_orig(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(q_opt(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz))
        allocate(mfx(bd%is:bd%ie+1,bd%js:bd%je,npz))
        allocate(mfy(bd%is:bd%ie,bd%js:bd%je+1,npz))
        allocate(cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz))
        allocate(cy(bd%isd:bd%ied,bd%js:bd%je+1,npz))

        ! Initialize test data with known values
        call initialize_test_data(q_orig, dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)
        q_opt = q_orig

        ! Run original and optimized implementations
        call tracer_2d(q_orig, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                      npx, npy, npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)

        ! Reset dp1, mfx, mfy, cx, cy since they may have been modified
        call initialize_test_data_params(dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)
        call tracer_2d_optimized(q_opt, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                npx, npy, npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)

        ! Calculate difference between original and optimized results
        error_norm = 0.0
        max_error = 0.0
        do iq = 1, nq
            do k = 1, npz
                do j = bd%js, bd%je
                    do i = bd%is, bd%ie
                        rel_error = abs(q_orig(i,j,k,iq) - q_opt(i,j,k,iq))
                        error_norm = error_norm + rel_error**2
                        if (rel_error > max_error) max_error = rel_error
                    end do
                end do
            end do
        end do

        error_norm = sqrt(error_norm / (nq * npz * (bd%ie-bd%is+1) * (bd%je-bd%js+1)))
        measured_accuracy_error = max(max_error, error_norm)

        write(*,'(A, E15.7, A, E15.7)') &
            ' Accuracy - Max error: ', max_error, &
            ', RMS error: ', error_norm

        ! Check if accuracy meets baseline expectations
        if (measured_accuracy_error > ACCURACY_DEGRADATION_THRESHOLD) then
            test_passed = .false.
            write(*,'(A, E15.7, A, E15.7)') &
                '  ACCURACY REGRESSION: Expected error <= ', ACCURACY_DEGRADATION_THRESHOLD, &
                ', but got: ', measured_accuracy_error
        else
            write(*,'(A, E15.7, A, E15.7)') &
                '  Accuracy meets threshold: Expected <= ', ACCURACY_DEGRADATION_THRESHOLD, &
                ', Got: ', measured_accuracy_error
        endif

        ! Run mass conservation regression test
        if (.not. run_mass_conservation_regression_test(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)) then
            test_passed = .false.
        endif

        deallocate(q_orig, q_opt, dp1, mfx, mfy, cx, cy)
    end function run_accuracy_regression_tests

    ! Run mass conservation regression test
    function run_mass_conservation_regression_test(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac) result(test_passed)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(fv_grid_bounds_type), intent(IN) :: bd
        type(domain2d), intent(INOUT) :: domain
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac
        logical :: test_passed

        real, allocatable :: q(:,:,:,:), dp1(:,:,:)
        real, allocatable :: mfx(:,:,:), mfy(:,:,:)
        real, allocatable :: cx(:,:,:), cy(:,:,:)
        type(group_halo_update_type) :: q_pack, dp1_pack
        real :: initial_mass, final_mass, conservation_error
        integer :: i, j, k, iq

        test_passed = .true.

        write(*,*) '  Running mass conservation regression test...'

        allocate(q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz))
        allocate(mfx(bd%is:bd%ie+1,bd%js:bd%je,npz))
        allocate(mfy(bd%is:bd%ie,bd%js:bd%je+1,npz))
        allocate(cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz))
        allocate(cy(bd%isd:bd%ied,bd%js:bd%je+1,npz))

        ! Initialize test data
        call initialize_test_data(q, dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)

        ! Calculate initial mass for first tracer
        initial_mass = 0.0
        do k = 1, npz
            do j = bd%js, bd%je
                do i = bd%is, bd%ie
                    initial_mass = initial_mass + q(i,j,k,1) * dp1(i,j,k)
                end do
            end do
        end do

        ! Run optimized tracer transport
        call tracer_2d_optimized(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                npx, npy, npz, 1, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)

        ! Calculate final mass
        final_mass = 0.0
        do k = 1, npz
            do j = bd%js, bd%je
                do i = bd%is, bd%ie
                    final_mass = final_mass + q(i,j,k,1) * dp1(i,j,k)
                end do
            end do
        end do

        conservation_error = abs(final_mass - initial_mass) / abs(initial_mass)

        write(*,'(A, E15.7)') '    Mass conservation error: ', conservation_error

        if (conservation_error > BASELINE_MASS_CONSERVATION_ERROR) then
            test_passed = .false.
            write(*,'(A, E15.7, A, E15.7)') &
                '    MASS CONSERVATION REGRESSION: Expected error <= ', BASELINE_MASS_CONSERVATION_ERROR, &
                ', but got: ', conservation_error
        else
            write(*,'(A, E15.7, A, E15.7)') &
                '    Mass conservation meets threshold: Expected <= ', BASELINE_MASS_CONSERVATION_ERROR, &
                ', Got: ', conservation_error
        endif

        deallocate(q, dp1, mfx, mfy, cx, cy)
    end function run_mass_conservation_regression_test

    ! Run stability regression tests
    function run_stability_regression_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac) result(test_passed)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(fv_grid_bounds_type), intent(IN) :: bd
        type(domain2d), intent(INOUT) :: domain
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac
        logical :: test_passed

        real, allocatable :: q(:,:,:,:), dp1(:,:,:)
        real, allocatable :: mfx(:,:,:), mfy(:,:,:)
        real, allocatable :: cx(:,:,:), cy(:,:,:)
        type(group_halo_update_type) :: q_pack, dp1_pack
        real :: min_val, max_val
        integer :: i, j, k, iq, iter, max_iter = 10
        real :: cfl_max

        test_passed = .true.

        write(*,*) 'Running stability regression tests...'

        ! Allocate test arrays
        allocate(q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq))
        allocate(dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz))
        allocate(mfx(bd%is:bd%ie+1,bd%js:bd%je,npz))
        allocate(mfy(bd%is:bd%ie,bd%js:bd%je+1,npz))
        allocate(cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz))
        allocate(cy(bd%isd:bd%ied,bd%js:bd%je+1,npz))

        ! Initialize test data
        call initialize_test_data(q, dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)

        ! Calculate maximum CFL number
        cfl_max = 0.0
        do k = 1, npz
            do j = bd%js, bd%je
                do i = bd%is, bd%ie+1
                    if (abs(cx(i,j,k)) > cfl_max) cfl_max = abs(cx(i,j,k))
                end do
            do j = bd%js, bd%je+1
                do i = bd%is, bd%ie
                    if (abs(cy(i,j,k)) > cfl_max) cfl_max = abs(cy(i,j,k))
                end do
            end do
        end do

        write(*,'(A, F10.6)') '  Maximum CFL in test: ', cfl_max

        ! Run multiple iterations to test stability
        do iter = 1, max_iter
            ! Run optimized tracer
            call tracer_2d_optimized(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                    npx, npy, npz, nq, hord, 0, dt, -1, q_pack, dp1_pack, 0, 0.0, lim_fac)

            ! Check for stability (finite values, no NaN or Inf)
            min_val = huge(1.0)
            max_val = -huge(1.0)

            do iq = 1, nq
                do k = 1, npz
                    do j = bd%js, bd%je
                        do i = bd%is, bd%ie
                            if (q(i,j,k,iq) /= q(i,j,k,iq)) then  ! Check for NaN
                                test_passed = .false.
                                write(*,'(A, I0, A, I0, A, I0, A, I0)') &
                                    '  STABILITY REGRESSION: NaN detected at iteration ', iter, &
                                    ', tracer ', iq, ', level ', k, ', grid point (', i, ',', j, ')'
                                goto 999
                            endif
                            if (abs(q(i,j,k,iq)) > huge(1.0)/2.0) then  ! Check for Inf
                                test_passed = .false.
                                write(*,'(A, I0, A, I0, A, I0, A, I0)') &
                                    '  STABILITY REGRESSION: Inf detected at iteration ', iter, &
                                    ', tracer ', iq, ', level ', k, ', grid point (', i, ',', j, ')'
                                goto 999
                            endif

                            if (q(i,j,k,iq) < min_val) min_val = q(i,j,k,iq)
                            if (q(i,j,k,iq) > max_val) max_val = q(i,j,k,iq)
                        end do
                    end do
                end do
            end do

            write(*,'(A, I0, A, E15.7, A, E15.7)') &
                '    Iteration ', iter, ': Min = ', min_val, ', Max = ', max_val
        end do

 999     continue

        if (test_passed) then
            write(*,*) '  Stability regression tests PASSED'
        else
            write(*,*) '  Stability regression tests FAILED'
        endif

        deallocate(q, dp1, mfx, mfy, cx, cy)
    end function run_stability_regression_tests

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

        call initialize_test_data_params(dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)
    end subroutine initialize_test_data

    ! Initialize parameter arrays
    subroutine initialize_test_data_params(dp1, mfx, mfy, cx, cy, bd, npx, npy, npz, nq)
        real, intent(out) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
        real, intent(out) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,npz)
        real, intent(out) :: mfy(bd%is:bd%ie,bd%js:bd%je+1,npz)
        real, intent(out) :: cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz)
        real, intent(out) :: cy(bd%isd:bd%ied,bd%js:bd%je+1,npz)
        type(fv_grid_bounds_type), intent(IN) :: bd
        integer, intent(IN) :: npx, npy, npz, nq

        integer :: i, j, k

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
    end subroutine initialize_test_data_params

    ! Helper function to convert integer to string
    function str(i) result(s)
        integer, intent(in) :: i
        character(len=20) :: s
        write(s, '(I0)') i
        s = adjustl(s)
    end function str

end module regression_test_suite_mod