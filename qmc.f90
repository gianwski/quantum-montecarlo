! Author: Gianluca Regni
! Copyright (c) 2024 Gianluca Regni
! License: GPL-3.0 (https://www.gnu.org/licenses/)
! Credits: Please cite the author if you use or reference this program.

! Description: 
! The program employs Quantum Monte Carlo (QMC) techniques to estimate 
! the average energy and error of quantum systems, specifically
! developed to accommodate five key cases: H, He, H2, H2+, and H3+.
! It supports Variational Monte Carlo (VMC) and Pure Diffusion Monte
! Carlo (PDMC) methods for simulations.
! The program is a compination of two codes: a main program (qmc.f90)
! and a module with specific routines for QMC calculations (qmc_routines.f90).

! Usage (compilation and run):
! gfortran -c qmc_routines.f90
! gfortran qmc_routines.o qmc.f90 -o qmc
! ./qmc


program qmc

  use qmc_routines

  implicit none

  ! Parameter for unit conversions to atomic unit

  double precision, parameter  :: ang_to_bohr = 1.8897259886
  double precision, parameter  :: bohr_to_ang = 0.529177249

  ! Definition of variables

  character(len=4)     :: method
  double precision     :: a, dt, E_ref, tau
  double precision     :: ave_e, ave_a, err_e, err_a
  integer*8            :: nmax
  integer              :: nruns, irun, i, natoms, charge, nelec

  double precision, allocatable :: X(:), accep(:), geom(:,:)
  character(len=2), allocatable :: atom(:)
  integer         , allocatable :: z(:)

  open (2,file="qmc.out")

  ! Read input variables

  call read_input(method,a,dt,E_ref,tau,nruns,nmax,natoms,atom,geom,charge)

  ! Check cases (H, H2, H2+, H3+, He)

  if (nelec > 2) stop "Error: number of electron must be < 3"
  if (natoms > 3) stop "Error: number of atoms must be < 4"

  ! Convert nuclei position in atomic unit

  geom = geom * ang_to_bohr

  ! Assign atomic number to elements

  allocate(z(natoms))
  
  call find_atomic_number(atom,natoms,z)

  ! Calculate the number of electrons

  nelec = sum(z) - charge

  ! Calculation using the Monte Carlo method specified in the input

  allocate(X(nruns), accep(nruns))

  if (method == "vmc") then

      do irun = 1, nruns
        call vmc(a,dt,nmax,X(irun),accep(irun),natoms,geom,nelec,z)
      enddo

  elseif (method == "pdmc") then

      do irun = 1, nruns
        call pdmc(a,dt,nmax,X(irun),accep(irun),tau,E_ref,natoms,geom,nelec,z)
      enddo

  else

      print *, "Error: method is not set correctly"
      stop

  endif

  ! Compute average energy and acceptance ratio with respective errors

  call ave_error(X,nruns,ave_e,err_e)
  call ave_error(accep,nruns,ave_a,err_a)

  ! Write outputs

  geom = bohr_to_ang * geom

  call output(method,a,dt,E_ref,tau,nruns,nmax,natoms,atom,geom,z,charge, &
              ave_e,ave_a,err_e,err_a)

end program qmc