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

program run_comprehensive_tests
    use benchmark_test_suite_mod
    use performance_measurement_mod
    use fv_arrays_mod, only: fv_grid_type, fv_grid_bounds_type
    use mpp_domains_mod, only: domain2d
    use mpp_mod, only: mpp_error, FATAL, mpp_pe, mpp_npes, mpp_init, mpp_sync

    implicit none

    ! Grid and domain variables
    type(fv_grid_type) :: gridstruct
    type(fv_grid_bounds_type) :: bd
    type(domain2d) :: domain
    integer :: npx, npy, npz, nq, hord
    real :: dt, lim_fac

    ! Initialize MPI and model components
    call mpp_init()
    write(*,*) 'FV3 Tracer Transport Optimization - Comprehensive Test Suite'
    write(*,*) '==========================================================='
    write(*,'(A, I0, A, I0)') 'Running on ', mpp_pe(), ' of ', mpp_npes()-1, ' processors'

    ! Set up test parameters
    npx = 20    ! Number of grid points in x direction
    npy = 20    ! Number of grid points in y direction  
    npz = 10    ! Number of vertical levels
    nq = 5      ! Number of tracers
    hord = 5    ! Horizontal order of accuracy
    dt = 900.0 ! Time step (in seconds)
    lim_fac = 1.0 ! Limiter factor

    ! Initialize bounds structure
    bd%is = 2
    bd%ie = npx - 1
    bd%js = 2  
    bd%je = npy - 1
    bd%isd = 1
    bd%ied = npx
    bd%jsd = 1
    bd%jed = npy

    ! Initialize grid structure with basic values
    call initialize_grid_structure(gridstruct, bd, npx, npy)

    ! Initialize domain (simplified)
    call initialize_domain(domain, bd)

    ! Initialize performance timers
    call init_performance_timers()

    ! Run comprehensive benchmark tests
    call run_comprehensive_benchmark_tests(gridstruct, bd, domain, npx, npy, npz, nq, hord, dt, lim_fac)

    ! Synchronize and finalize
    call mpp_sync()
    write(*,*) 'All tests completed successfully.'

contains

    ! Initialize grid structure with basic values
    subroutine initialize_grid_structure(gridstruct, bd, npx, npy)
        type(fv_grid_type), intent(inout) :: gridstruct
        type(fv_grid_bounds_type), intent(in) :: bd
        integer, intent(in) :: npx, npy

        integer :: i, j

        ! Allocate and initialize grid arrays
        allocate(gridstruct%area(bd%isd:bd%ied, bd%jsd:bd%jed))
        allocate(gridstruct%rarea(bd%isd:bd%ied, bd%jsd:bd%jed))
        allocate(gridstruct%sin_sg(bd%isd:bd%ied, bd%jsd:bd%jed, 5))
        allocate(gridstruct%dxa(bd%is:bd%ie+1, bd%js:bd%je))
        allocate(gridstruct%dya(bd%is:bd%ie, bd%js:bd%je+1))
        allocate(gridstruct%dx(bd%is:bd%ie+1, bd%js:bd%je+1))
        allocate(gridstruct%dy(bd%is:bd%ie+1, bd%js:bd%je+1))

        ! Initialize with simple values
        do j = bd%jsd, bd%jed
            do i = bd%isd, bd%ied
                gridstruct%area(i,j) = 1.0
                gridstruct%rarea(i,j) = 1.0
                gridstruct%sin_sg(i,j,1) = 0.1
                gridstruct%sin_sg(i,j,2) = 0.1
                gridstruct%sin_sg(i,j,3) = 0.1
                gridstruct%sin_sg(i,j,4) = 0.1
                gridstruct%sin_sg(i,j,5) = 0.9
            end do
        end do

        do j = bd%js, bd%je
            do i = bd%is, bd%ie+1
                gridstruct%dxa(i,j) = 1.0
            end do
        end do

        do j = bd%js, bd%je+1
            do i = bd%is, bd%ie
                gridstruct%dya(i,j) = 1.0
            end do
        end do

        do j = bd%js, bd%je+1
            do i = bd%is, bd%ie+1
                gridstruct%dx(i,j) = 1.0
                gridstruct%dy(i,j) = 1.0
            end do
        end do

        ! Set other grid properties
        gridstruct%nested = .false.
        gridstruct%regional = .false.
    end subroutine initialize_grid_structure

    ! Initialize domain structure (simplified)
    subroutine initialize_domain(domain, bd)
        type(domain2d), intent(inout) :: domain
        type(fv_grid_bounds_type), intent(in) :: bd
        
        ! Initialize domain with dummy values
        ! In a real application, this would be properly initialized
    end subroutine initialize_domain

end program run_comprehensive_tests