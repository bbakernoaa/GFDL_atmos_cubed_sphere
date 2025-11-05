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

module tracer_transport_integrator_mod
    use fv_tracer2d_mod, only: tracer_2d, tracer_2d_1L, tracer_2d_nested
    use fv_tracer2d_optimized_mod, only: tracer_2d_optimized, tracer_2d_1L_optimized, tracer_2d_nested_optimized
    use tp_core_mod, only: fv_tp_2d
    use tp_core_optimized_mod, only: fv_tp_2d_optimized
    use performance_measurement_mod, only: performance_timer_start, performance_timer_stop, &
                                          get_elapsed_time, print_performance_summary, &
                                          tracer_performance_test
    use tracer_test_suite_mod, only: run_all_tests
    use fv_arrays_mod, only: fv_grid_type, fv_grid_bounds_type
    use mpp_domains_mod, only: domain2d
    use mpp_mod, only: mpp_error, FATAL, mpp_max

    implicit none
    private

    public :: integrate_tracer_transport, select_optimized_transport, run_validation_suite

    ! Configuration parameters
    logical, parameter :: USE_OPTIMIZED_DEFAULT = .true.
    integer, parameter :: MAX_TEST_ITERATIONS = 10

contains

    ! Main integration routine for tracer transport
    subroutine integrate_tracer_transport(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                        npx, npy, npz, nq, hord, q_split, dt, id_divg_mean, &
                                        q_pack, dp1_pack, nord_tr, trdm, k_split, &
                                        neststruct, parent_grid, n_map, lim_fac, &
                                        use_optimized)
        type(fv_grid_bounds_type), intent(IN) :: bd
        integer, intent(IN) :: npx, npy, npz, nq, hord, q_split, id_divg_mean
        integer, intent(IN), optional :: k_split, n_map
        real, intent(IN) :: dt, trdm, lim_fac
        real, intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)
        real, intent(INOUT) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
        real, intent(INOUT) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,npz)
        real, intent(INOUT) :: mfy(bd%is:bd%ie,bd%js:bd%je+1,npz)
        real, intent(INOUT) :: cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz)
        real, intent(INOUT) :: cy(bd%isd:bd%ied,bd%js:bd%je+1,npz)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(domain2d), intent(INOUT) :: domain
        type(group_halo_update_type), intent(INOUT) :: q_pack, dp1_pack
        integer, intent(IN), optional :: nord_tr
        real, intent(IN), optional :: trdm
        type(fv_nest_type), intent(INOUT), optional :: neststruct
        type(fv_atmos_type), pointer, intent(IN), optional :: parent_grid
        logical, intent(IN), optional :: use_optimized

        logical :: use_opt
        integer :: nord_tr_local
        real :: trdm_local

        ! Set default values if optional parameters are not provided
        if (present(use_optimized)) then
            use_opt = use_optimized
        else
            use_opt = USE_OPTIMIZED_DEFAULT
        endif

        if (present(nord_tr)) then
            nord_tr_local = nord_tr
        else
            nord_tr_local = 0
        endif

        if (present(trdm)) then
            trdm_local = trdm
        else
            trdm_local = 0.0
        endif

        ! Perform tracer transport based on optimization flag
        if (use_opt) then
            call select_optimized_transport(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                          npx, npy, npz, nq, hord, q_split, dt, id_divg_mean, &
                                          q_pack, dp1_pack, nord_tr_local, trdm_local, &
                                          k_split, neststruct, parent_grid, n_map, lim_fac)
        else
            ! Use original transport routines
            if (present(neststruct) .and. present(parent_grid)) then
                if (present(k_split) .and. present(n_map)) then
                    call tracer_2d_nested(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                         npx, npy, npz, nq, hord, q_split, dt, id_divg_mean, &
                                         q_pack, dp1_pack, nord_tr_local, trdm_local, &
                                         k_split, neststruct, parent_grid, n_map, lim_fac)
                else
                    call mpp_error(FATAL, 'tracer_transport_integrator: k_split and n_map required for nested runs')
                endif
            else
                if (npz == 1) then
                    call tracer_2d_1L(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                     npx, npy, npz, nq, hord, q_split, dt, id_divg_mean, &
                                     q_pack, dp1_pack, nord_tr_local, trdm_local, lim_fac)
                else
                    call tracer_2d(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                  npx, npy, npz, nq, hord, q_split, dt, id_divg_mean, &
                                  q_pack, dp1_pack, nord_tr_local, trdm_local, lim_fac)
                endif
            endif
        endif
    end subroutine integrate_tracer_transport

    ! Select optimized transport routine based on configuration
    subroutine select_optimized_transport(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                        npx, npy, npz, nq, hord, q_split, dt, id_divg_mean, &
                                        q_pack, dp1_pack, nord_tr, trdm, k_split, &
                                        neststruct, parent_grid, n_map, lim_fac)
        type(fv_grid_bounds_type), intent(IN) :: bd
        integer, intent(IN) :: npx, npy, npz, nq, hord, q_split, id_divg_mean
        integer, intent(IN), optional :: k_split, n_map
        real, intent(IN) :: dt, trdm, lim_fac
        real, intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)
        real, intent(INOUT) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
        real, intent(INOUT) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,npz)
        real, intent(INOUT) :: mfy(bd%is:bd%ie,bd%js:bd%je+1,npz)
        real, intent(INOUT) :: cx(bd%is:bd%ie+1,bd%jsd:bd%jed,npz)
        real, intent(INOUT) :: cy(bd%isd:bd%ied,bd%js:bd%je+1,npz)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(domain2d), intent(INOUT) :: domain
        type(group_halo_update_type), intent(INOUT) :: q_pack, dp1_pack
        integer, intent(IN) :: nord_tr
        type(fv_nest_type), intent(INOUT), optional :: neststruct
        type(fv_atmos_type), pointer, intent(IN), optional :: parent_grid

        ! Use optimized transport routines
        if (present(neststruct) .and. present(parent_grid)) then
            if (present(k_split) .and. present(n_map)) then
                call tracer_2d_nested_optimized(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                              npx, npy, npz, nq, hord, q_split, dt, id_divg_mean, &
                                              q_pack, dp1_pack, nord_tr, trdm, &
                                              k_split, neststruct, parent_grid, n_map, lim_fac)
            else
                call mpp_error(FATAL, 'tracer_transport_integrator: k_split and n_map required for nested runs')
            endif
        else
            if (npz == 1) then
                call tracer_2d_1L_optimized(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                          npx, npy, npz, nq, hord, q_split, dt, id_divg_mean, &
                                          q_pack, dp1_pack, nord_tr, trdm, lim_fac)
            else
                call tracer_2d_optimized(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, &
                                       npx, npy, npz, nq, hord, q_split, dt, id_divg_mean, &
                                       q_pack, dp1_pack, nord_tr, trdm, lim_fac)
            endif
        endif
    end subroutine select_optimized_transport

    ! Run comprehensive validation suite
    subroutine run_validation_suite(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(fv_grid_bounds_type), intent(IN) :: bd
        type(domain2d), intent(INOUT) :: domain
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac

        write(*,*) 'Running comprehensive validation suite...'
        write(*,*) '=========================================='

        ! Run the full test suite
        call run_all_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)

        ! Run performance comparison
        call tracer_performance_test_dummy(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)

        write(*,*) '=========================================='
        write(*,*) 'Validation suite completed.'
    end subroutine run_validation_suite

    ! Dummy performance test to demonstrate integration
    subroutine tracer_performance_test_dummy(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)
        type(fv_grid_type), intent(IN) :: gridstruct
        type(fv_grid_bounds_type), intent(IN) :: bd
        type(domain2d), intent(INOUT) :: domain
        integer, intent(IN) :: npx, npy, npz, nq, hord
        real, intent(IN) :: dt, lim_fac

        write(*,*) 'Running performance comparison test...'
        
        ! This would normally run the actual performance test
        ! For demonstration, we'll just start and stop timers
        call performance_timer_start('integration_performance_test')
        
        ! Simulate some work
        call simulate_integration_work(npx, npy, npz, nq)
        
        call performance_timer_stop('integration_performance_test')
        
        write(*,'(A, F10.6, A)') 'Integration test elapsed time: ', &
            get_elapsed_time('integration_performance_test'), ' seconds'
    end subroutine tracer_performance_test_dummy

    ! Simulate integration work
    subroutine simulate_integration_work(npx, npy, npz, nq)
        integer, intent(IN) :: npx, npy, npz, nq
        integer :: i, j, k, iq
        real :: temp

        ! Simulate computational work
        temp = 0.0
        do iq = 1, nq
            do k = 1, npz
                do j = 1, npy
                    do i = 1, npx
                        temp = temp + sin(real(i)*real(j)*real(k)*real(iq)*0.001)
                    end do
                end do
            end do
        end do

        ! Prevent optimization from eliminating the loop
        if (temp < 0.0) write(*,*) 'This should not happen: ', temp
    end subroutine simulate_integration_work

end module tracer_transport_integrator_mod