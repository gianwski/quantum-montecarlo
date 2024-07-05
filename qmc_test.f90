! Author: Gianluca Regni
! Copyright (c) 2024 Gianluca Regni
! License: GPL-3.0 (https://www.gnu.org/licenses/)
! Credits: Please cite the author if you use or reference this program.

! Description:
! This program is built on the principles of Test-Driven Development (TDD).
! It includes a suite of automated tests designed to validate the functionality 
! of quantum Monte Carlo (QMC) routines.
! The program tests various components, including psi functions, potential 
! energy calculations, kinetic energy calculations, and atomic number 
! identification

! Usage (compilation and run):
! gfortran -o qmc_test qmc_routines.f90 qmc_test.f90
! ./qmc_test

subroutine test_psi()

    use qmc_routines, only: psi

    implicit none

    double precision              :: a, expected_psi
    double precision, allocatable :: r(:,:), geom(:,:)
    integer                       :: natoms, nelec
     
    ! TEST 1
    a = 1.2d0
    natoms = 2
    nelec = 2

    allocate(r(nelec,3), geom(natoms,3))

    r = 1.d0
    geom = 0.d0

    expected_psi = (dexp(-1.2d0*dsqrt(3.d0)) + dexp(-1.2d0*dsqrt(3.d0))) * &
                   (dexp(-1.2d0*dsqrt(3.d0)) + dexp(-1.2d0*dsqrt(3.d0)))

    if (psi(a,r,natoms,geom,nelec) /= expected_psi) stop "PSI TEST FAILED"

    deallocate(r, geom)

    ! TEST 2
    a = 2.d0
    natoms = 2
    nelec = 2

    allocate(r(nelec,3), geom(natoms,3))

    r = 1.d0
    geom = 1.5d0

    expected_psi = (dexp(-2.d0*dsqrt(0.75d0)) + dexp(-2.d0*dsqrt(0.75d0))) * &
                   (dexp(-2.d0*dsqrt(0.75d0)) + dexp(-2.d0*dsqrt(0.75d0)))

    if (psi(a,r,natoms,geom,nelec) /= expected_psi) stop "PSI TEST FAILED"
 
    deallocate(r, geom)

    ! TEST 3
    a = 1.5d0
    natoms = 1
    nelec = 2

    allocate(r(nelec,3), geom(natoms,3))

    r = 0.d0
    geom = 2.d0

    expected_psi = dexp(-1.5d0*dsqrt(12.d0)) * dexp(-1.5d0*dsqrt(12.d0))

    if (psi(a,r,natoms,geom,nelec) /= expected_psi) stop "PSI TEST FAILED"

    print *, "1) Psi test passed!"


end subroutine test_psi

!*******************************************************************************

subroutine test_potential()

    use qmc_routines, only: potential

    implicit none

    double precision              :: expected_potential
    double precision, allocatable :: r(:,:), geom(:,:)
    integer         , allocatable :: z(:)
    integer                       :: natoms, nelec
 
    ! TEST 1
    natoms = 1
    nelec = 1
    
    allocate(r(nelec,3), geom(natoms,3), z(natoms))

    expected_potential = -1.d0/15.d0

    r(1,:) = (/ 2.d0, 5.d0, 14.d0 /)
    geom = 0.d0
    z = 1

    if (potential(r,natoms,geom,nelec,z) /= expected_potential) stop "POTENTIAL TEST FAILED"

    deallocate(r, geom, z)

    ! TEST 2
    natoms = 2
    nelec = 2

    allocate(r(nelec,3), geom(natoms,3), z(natoms))

    expected_potential = -1.0D+308

    r = 0.d0
    geom = 0.d0
    z = 1

    if (potential(r,natoms,geom,nelec,z) >= expected_potential) stop "POTENTIAL TEST FAILED"


    print *, "2) Potential test passed!"


end subroutine test_potential

!*******************************************************************************

subroutine test_psi1e1N()

    use qmc_routines, only: psi_1e1N

    implicit none

    double precision              :: expected_psi1e1N, a
    double precision, allocatable :: r(:,:), geom(:,:)
    integer                       :: natoms, nelec, i, j
 
    ! TEST 1
    a = 1.2d0
    natoms = 2
    nelec = 2
    i = 1
    j = 1

    allocate(r(nelec,3), geom(natoms,3))

    expected_psi1e1N = dexp(-1.2d0 * dsqrt(3.d0))

    r = 1.d0
    geom = 0.d0

    if (psi_1e1N(a,r,geom,nelec,natoms,i,j) /= expected_psi1e1N) stop "PSI_1E1N TEST FAILED"
    
    ! TEST 2

    r = 2.d0
    geom = 1.d0

    if (psi_1e1N(a,r,geom,nelec,natoms,i,j) /= expected_psi1e1N) stop "PSI_1E1N TEST FAILED"

    ! TEST 3

    r = 3.d0
    geom = 2.d0
    i = 2
    j = 2

    if (psi_1e1N(a,r,geom,nelec,natoms,i,j) /= expected_psi1e1N) stop "PSI_1E1N TEST FAILED"


    print *, "3) Psi_1e1N test passed!"

end subroutine test_psi1e1N

!*******************************************************************************

subroutine test_kinetic()

use qmc_routines, only: kinetic

    implicit none

    double precision              :: expected_kinetic, a
    double precision, allocatable :: r(:,:), geom(:,:)
    integer                       :: natoms, nelec

    ! TEST 1
    a = 1.0d0
    natoms = 1
    nelec = 1

    allocate(r(nelec,3),geom(natoms,3))

    expected_kinetic = 0.d0

    r = 0.d0
    geom = 0.d0

    r(1,1) = 2.d0

    if (kinetic(a,r,natoms,geom,nelec) /= expected_kinetic) stop "KINETIC TEST FAILED"

    deallocate(r,geom)

    ! TEST 2
    natoms = 2
    nelec = 2

    allocate(r(nelec,3),geom(natoms,3))

    expected_kinetic = 0.d0

    r = 0.d0
    geom = 0.d0

    r(1,1) = 2.d0
    r(2,1) = 2.d0

    if (kinetic(a,r,natoms,geom,nelec) /= expected_kinetic) stop "KINETIC TEST FAILED"

    print *, "4) Kinetic test passed!"

end subroutine test_kinetic

!*******************************************************************************

subroutine test_find_atomic_number()

    use qmc_routines, only: find_atomic_number

    implicit none

    character(len=2), allocatable  :: input_element(:)
    integer                        :: nelem, i
    integer, allocatable           :: z(:), expected_z(:)

    ! TEST 1

    nelem = 2

    allocate(input_element(nelem),z(nelem),expected_z(nelem))

    input_element = (/ "H ", "He" /)

    call find_atomic_number(input_element,nelem,z)

    expected_z = (/ 1, 2 /)

    do i = 1, nelem
        if (z(i) /= expected_z(i)) stop "FIND ATOMIC NUMBER TEST FAILED"
    enddo

    deallocate(input_element,z,expected_z)

    ! TEST 2

    nelem = 2

    allocate(input_element(nelem),z(nelem),expected_z(nelem))

    input_element = (/ "H ", " H" /)

    call find_atomic_number(input_element,nelem,z)

    expected_z = (/ 1, 1 /)

    do i = 1, nelem
        if (z(i) /= expected_z(i)) stop "FIND ATOMIC NUMBER TEST FAILED"
    enddo

    deallocate(input_element,z,expected_z)

    ! TEST 3

    nelem = 1

    allocate(input_element(nelem),z(nelem),expected_z(nelem))

    input_element = (/ "A " /)

    call find_atomic_number(input_element,nelem,z)

    expected_z = (/ 0 /)

    do i = 1, nelem
        if (z(i) /= expected_z(i)) stop "FIND ATOMIC NUMBER TEST FAILED"
    enddo

    print *, "5) Find atomic number test passed!"


end subroutine test_find_atomic_number

!*******************************************************************************

program qmc_test

    ! This program runs all the test subroutines for psi, potential,
    ! psi_1e1N, kinetic, and find_atomic_number.

    implicit none

    call test_psi()
    call test_potential()
    call test_psi1e1N()
    call test_kinetic()
    call test_find_atomic_number()
    
end program qmc_test
