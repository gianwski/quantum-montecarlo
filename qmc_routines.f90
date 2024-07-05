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

module qmc_routines

   public :: psi
   public :: potential
   public :: kinetic
   public :: e_loc
   public :: ave_error
   public :: random_gauss
   public :: drift
   public :: vmc
   public :: pdmc
   public :: read_input
   public :: output
   public :: find_atomic_number
   public :: psi_1e1N

contains

!****************************************************************************

function psi(a,r,natoms,geom,nelec)

   ! The function calculates and returns the wave function value of a specific system
   ! configuration


   ! a      : wavefunction parameter
   ! r      : electron positions
   ! natoms : number of atoms (or nuclei)
   ! geom   : nuclei positions
   ! nelec  : number of electrons
   ! psi    : wavefunction


   implicit none

   double precision, intent(in) :: a, r(nelec,3), geom(natoms,3)
   integer         , intent(in) :: natoms, nelec
   integer                      :: i, j
   double precision             :: psi, e, distance

   psi = 1.d0

   do i = 1, nelec

      e = 0.d0

      ! This loop iterates over each atom (nucleus) for a given electron 'i'.
      ! For each pair of electron 'i' and nucleus 'j', it calculates the exponential
      ! decay of their interaction based on the distance between them. This is done
      ! by first computing the Euclidean distance between electron 'i' and nucleus 'j',
      ! then applying an exponential decay function modulated by the parameter 'a'.
      ! The result reflects the contribution of nucleus 'j' to the electron's state.
      ! These contributions are summed over all nuclei to get a total interaction effect
      ! for electron 'i', stored in 'e'.

      do j = 1, natoms
         e = e + dexp(-a * dsqrt((r(i,1)-geom(j,1))**2.d0 + &
                                 (r(i,2)-geom(j,2))**2.d0 + &
                                 (r(i,3)-geom(j,3))**2.d0))
      enddo


      ! After summing the contributions from all nuclei for electron 'i',
      ! this total is then multiplied with the running product stored in 'psi'.
      ! This step effectively accumulates the interaction effect across all electrons,
      ! considering their respective interactions with all nuclei. The final value
      ! of 'psi' after completing the loop over all electrons represents the
      ! wave function of the system.
      ! The psi is approximated as the product of one-electron wave functions.

      psi = psi * e

   enddo

   end function psi

!****************************************************************************

function potential(r,natoms,geom,nelec,z)

   ! This function calculates the local total potential energy of a quantum system consisting
   ! of electrons and nuclei. It includes contributions from the attractive interactions
   ! between electrons and nuclei, as well as repulsive interactions among electrons and
   ! among nuclei. These interactions are modeled according to Coulomb's law, where the
   ! attractive and repulsive energies are inversely proportional to the distance between
   ! charges.

   ! r              : electron positions
   ! natoms         : number of atoms (or nuclei)
   ! geom           : nuclei positions
   ! nelec          : number of electrons
   ! z              : atomic number
   ! distance       : distance between electron and nucleus
   ! attractive     : attractive part of coulomb interaction
   ! repulsive_elec : electron repulsive part of coulomb interaction
   ! repulsive_nuc  : nuclear repulsive part of coulomb interaction
   ! potential      : local total potential energy

   implicit none

   double precision, intent(in) :: r(nelec,3), geom(natoms,3)
   integer         , intent(in) :: natoms, nelec, z(natoms)
   double precision             :: potential, distance, attractive, repulsive_elec, repulsive_nuc, eps
   integer                      :: i, j

   ! Calculation of attractive part of interaction (nucleus-electron)

   eps = 1.d-308       ! using huge(1.d0) command can lead to overflow
   attractive = 0.d0

   do i = 1, natoms
      do j = 1, nelec

         distance = dsqrt((geom(i,1)-r(j,1))**2.d0 + &
                           (geom(i,2)-r(j,2))**2.d0 + &
                           (geom(i,3)-r(j,3))**2.d0)

         if (distance > 0.d0) then
            attractive = attractive - z(i) / distance
         else
            print *, 'Warning: zero distance encountered'
            attractive = attractive - z(i) / eps
         end if

      enddo
   enddo

   ! Calculation of repulsive part of interaction (electron-electron)

   repulsive_elec = 0.d0

   if (nelec > 1) then
      do i = 1, nelec-1
         do j = i+1, nelec  ! avoid double counting and interaction with itself
            distance = dsqrt((r(i,1)-r(j,1))**2.d0 + &
                              (r(i,2)-r(j,2))**2.d0 + &
                              (r(i,3)-r(j,3))**2.d0)

            if (distance > 0.d0) then
               repulsive_elec = repulsive_elec + 1.d0 / distance
            else
               print *, 'Warning: zero distance encountered'
               repulsive_elec = repulsive_elec + 1.d0 / eps
            end if

         enddo
      enddo
   endif

   ! Calculation of repulsive part of interaction (nucleus-nucleus)

   repulsive_nuc = 0.d0

   if (natoms > 1) then
      do i = 1, natoms-1
         do j = i+1, natoms  ! avoid double counting and interaction with itself
            distance = dsqrt((geom(i,1)-geom(j,1))**2.d0 + &
                              (geom(i,2)-geom(j,2))**2.d0 + &
                              (geom(i,3)-geom(j,3))**2.d0)

            if (distance > 0.d0) then
               repulsive_nuc = repulsive_nuc + z(i)*z(j) / distance
            else
               print *, 'Warning: zero distance encountered'
               repulsive_nuc = repulsive_nuc + z(i)*z(j) / eps
            end if

         enddo
      enddo
   endif

   !  Calculation of the potential

   potential = attractive + repulsive_elec + repulsive_nuc

end function potential

 !****************************************************************************

function psi_1e1N(a,r,geom,nelec,natoms,i,j)

   ! This function calculates the contribution to the wavefunction of one electron and one nucleus.
   ! It is used to calculate the derivatives (or Laplacian) of psi for the calculation of local 
   ! kinetic energy and drift vector.

   ! a        : wavefunction parameter
   ! r        : electron positions
   ! geom     : nuclei positions
   ! nelec    : number of electrons
   ! natoms   : number of atoms (or nuclei)
   ! i        : electron index
   ! j        : nucleus index
   ! psi_1e1N : exp(-a * |r_e-R_N|)

   implicit none

   double precision, intent(in) :: a, r(nelec,3), geom(natoms,3)
   integer         , intent(in) :: i, j, nelec, natoms
   double precision             :: psi_1e1N


   psi_1e1N = dexp(-a*dsqrt((r(i,1)-geom(j,1))**2.d0 + &
                            (r(i,2)-geom(j,2))**2.d0 + &
                            (r(i,3)-geom(j,3))**2.d0))

end function psi_1e1N


!****************************************************************************

function kinetic(a,r,natoms,geom,nelec)

   ! This function computes local kinetic energy for the 5 key cases: H, H2+, He, H2, H3+.
   ! It calculates the local kinetic energy of the wave function on a case-by-case basis.

   ! Input parameters
   ! a        : wavefunction parameter
   ! r        : electron positions
   ! geom     : nuclei positions
   ! nelec    : number of electrons
   ! natoms   : number of atoms (or nuclei)

   ! Output parameters
   ! kinetic : local kinetic energy

   implicit none

   double precision, intent(in) :: a, r(nelec,3), geom(natoms,3)
   integer         , intent(in) :: natoms, nelec
   double precision             :: distance(nelec,natoms), prod(nelec,natoms), kinetic
   double precision             :: term, sum_prod
   integer                      :: i, j, k

   ! Compute distances and products
   do i = 1, nelec
      do j = 1, natoms
         distance(i,j) = dsqrt(sum((r(i,:) - geom(j,:))**2.d0))
         prod(i,j) = psi_1e1N(a,r,geom,nelec,natoms,1,i) * psi_1e1N(a,r,geom,nelec,natoms,2,j)
      enddo
   enddo

   kinetic = 0.d0

   ! Calculate kinetic energy based on the specific case determined by the number of electrons and atoms.
   select case (nelec * 10 + natoms) 

   ! H, H2+ cases
   case(1,12)
   do j = 1, natoms
      kinetic = kinetic + (a*a - 2*a/distance(1,j)) * psi_1e1N(a,r,geom,nelec,natoms,1,j)
   enddo

   ! He, H2, H3+ cases
   case(21,22,23)

   do i = 1, nelec
      do j = 1, natoms
         term = (a*a - 2.d0*a/distance(i,j))
         sum_prod = 0.d0
         do k = 1, natoms
            if (i == 1) then
               sum_prod = sum_prod + prod(j,k)
            elseif (i == 2) then
               sum_prod = sum_prod + prod(k,j) 
            endif
         enddo
         kinetic = kinetic + term * sum_prod
      enddo
   enddo

   end select

   kinetic = -0.5d0 * kinetic / psi(a,r,natoms,geom,nelec)

end function kinetic

!****************************************************************************

function e_loc(a,r,natoms,geom,nelec,z)

   ! The function calculates the local energy of a specific system
   ! configuration

   ! kinetic   : local kinetic energy
   ! potential : local total potential energy
   ! e_loc     : local energy

   implicit none

   double precision, intent(in) :: a, r(nelec,3), geom(natoms,3)
   integer         , intent(in) :: natoms, nelec, z(natoms)
   double precision             :: e_loc

   e_loc = kinetic(a,r,natoms,geom,nelec) + potential(r,natoms,geom,nelec,z)

end function e_loc

!****************************************************************************

subroutine ave_error(x,n,ave,err)

   ! This function returns the average and the statistical error of an input
   ! array x

   ! x        : the data set for which the average and error are to be computed
   ! n        : the number of elements in the array x
   ! ave      : the calculated average of the input array
   ! err      : the calculated statistical error of the average
   ! variance : the variance of the input data set

   implicit none

   integer         , intent(in)  :: n 
   double precision, intent(in)  :: x(n) 
   double precision, intent(out) :: ave, err
   double precision              :: variance

   ! Terminates the program if the array size is less than 1
   if (n < 1) then             
      stop 'n<1 in ave_error'

   ! Sets the average to the only element if array size is 1 and error to 0
   else if (n == 1) then       
      ave = x(1)
      err = 0.d0

   ! Calculates the average, variance, and error for arrays with more than one element
   else                        
      ave      = sum(x(:)) / dble(n)
      variance = sum((x(:) - ave)**2) / dble(n-1)
      err      = dsqrt(variance/dble(n))
      
   endif

end subroutine ave_error

!****************************************************************************

subroutine random_gauss(z,n,nset)

   ! This subroutine generates a matrix 'z' of Gaussian-distributed random numbers 
   ! for 'nset' sets, each with 'n' numbers, using the Box-Muller transform.

   ! n     : number of random numbers per set
   ! nset  : number of sets
   ! z     : output matrix (nset by n) of Gaussian-distributed random numbers

   implicit none

   double precision, parameter   :: two_pi = 2.d0*dacos(-1.d0)

   integer         , intent(in)  :: n, nset
   double precision, intent(out) :: z(nset,n)
   double precision              :: u(n+1)
   integer                       :: i, j

   do j = 1, nset

      call random_number(u)

      if (iand(n,1) == 0) then
         ! n is even
         do i=1,n,2
            z(j,i)   = dsqrt(-2.d0*dlog(u(i))) 
            z(j,i+1) = z(j,i) * dsin( two_pi*u(i+1) )
            z(j,i)   = z(j,i) * dcos( two_pi*u(i+1) )
         end do
      else
         ! n is odd
         do i=1,n-1,2
            z(j,i)   = dsqrt(-2.d0*dlog(u(i))) 
            z(j,i+1) = z(j,i) * dsin( two_pi*u(i+1) )
            z(j,i)   = z(j,i) * dcos( two_pi*u(i+1) )
         end do
         z(j,n)   = dsqrt(-2.d0*dlog(u(n))) 
         z(j,n)   = z(j,n) * dcos( two_pi*u(n+1) )
      end if

   enddo

end subroutine random_gauss

!****************************************************************************

subroutine drift(a,r,D,nelec,natoms,geom)

   ! This function computes the drift vector

   ! Input parameters
   ! a        : wavefunction parameter
   ! r        : electron positions
   ! geom     : nuclei positions
   ! nelec    : number of electrons
   ! natoms   : number of atoms (or nuclei)

   ! Output parameters
   ! D : drift vector

   implicit none

   integer         , intent(in)  :: nelec, natoms
   double precision, intent(in)  :: a, r(nelec,3), geom(natoms,3)
   double precision, intent(out) :: D(nelec,3)
   double precision              :: distance(nelec,natoms), prod(nelec,natoms)
   double precision              :: sum_prod, term(3)
   integer                       :: i, j, k


   ! Compute distances and products
   do i = 1, nelec
      do j = 1, natoms
         distance(i,j) = dsqrt(sum((r(i,:) - geom(j,:))**2.d0))
         prod(i,j) = psi_1e1N(a,r,geom,nelec,natoms,1,i) * psi_1e1N(a,r,geom,nelec,natoms,2,j)
      enddo
   enddo

   D = 0.d0

   ! Calculate kinetic energy based on the specific case determined by the number of electrons and atoms
   select case (nelec * 10 + natoms) 

   ! H, H2+ cases
   case(11,12)

      do j = 1, natoms
            D(1,:) =  D(1,:) - a * dabs(r(1,:)-geom(j,:)) / distance(1,j) * psi_1e1N(a,r,geom,nelec,natoms,1,j)
      enddo

   ! He, H2, H3+ cases
   case(21,22,23)

   do i = 1, nelec
      do j = 1, natoms
         term(:) = -a * dabs(r(i,:)-geom(j,:)) / distance(i,j)
         sum_prod = 0.d0
         do k = 1, natoms
            if (i == 1) then
               sum_prod = sum_prod + prod(j,k)
            elseif (i == 2) then
               sum_prod = sum_prod + prod(k,j)
            endif
         enddo
         D(i,:) = D(i,:) + term(:) * sum_prod
      enddo
   enddo

   end select
   
   D(:,:) = D(:,:) / psi(a,r,natoms,geom,nelec)

end subroutine drift

!****************************************************************************

subroutine vmc(a,dt,nmax,energy,accep,natoms,geom,nelec,z)

   ! This subroutine implements the Variational Monte Carlo (VMC) method to
   ! estimate the average energy

   implicit none

   ! Input parameters:
   ! a      : parameter of wave function
   ! dt     : step size
   ! geom   : geometry of system (xyz coordinates)
   ! nmax   : Monte Carlo steps
   ! natoms : number of atoms in the system
   ! nelec  : number of electrons in the system
   ! z      : atomic numbers of the atoms in the system

   ! Output parameters:
   ! energy : estimated average energy of the system
   ! accep  : acceptance ratio of the Monte Carlo step simualtion


   double precision, intent(in)  :: a, dt, geom(natoms,3)
   integer*8       , intent(in)  :: nmax 
   integer         , intent(in)  :: natoms, nelec, z(natoms)
   double precision, intent(out) :: energy, accep

   integer*8                            :: istep, n_accep
   integer                              :: i
   double precision                     :: sq_dt, u, psi_old, psi_new, argexpo, q
   double precision, dimension(nelec,3) :: r_old, r_new, d_old, d_new, chi
   double precision, dimension(nelec)   :: prod, d2_old, d2_new

   sq_dt = dsqrt(dt)

   ! Initialization
   energy  = 0.d0
   n_accep = 0_8

   ! Generate inital gaussian-distributed 3-dimensional vector for electron positions
   call random_gauss(r_old,3,nelec)

   ! Calculates initial drift vector for initial positions
   call drift(a,r_old,d_old,nelec,natoms,geom)

   ! Calculate square of drift vector d_old
   d2_old(:) = d_old(:,1)*d_old(:,1) + &
               d_old(:,2)*d_old(:,2) + &
               d_old(:,3)*d_old(:,3)

   ! Evaluate initial wave function value for the system
   psi_old = psi(a,r_old,natoms,geom,nelec)

   ! Main Monte Carlo loop
   do istep = 1, nmax

      ! Evaluate local energy and accumulate it
      energy = energy + e_loc(a,r_old,natoms,geom,nelec,z)

      ! Propose a new configuration by moving electrons
      call random_gauss(chi,3,nelec)
      do i = 1, nelec
         r_new(i,:) = r_old(i,:) + dt*d_old(i,:) + chi(i,:)*sq_dt
      enddo

      ! Calculates drift vector d_new for new position
      call drift(a,r_new,d_new,nelec,natoms,geom)

      ! Calculate square of the new drift vector (d_new)
      d2_new(:) = d_new(:,1)*d_new(:,1) + &
                  d_new(:,2)*d_new(:,2) + &
                  d_new(:,3)*d_new(:,3)

      ! Evaluate wave function for the new configuration
      psi_new = psi(a,r_new,natoms,geom,nelec)

      ! Metropolis acceptance criterion calculation
      prod(:) = (d_new(:,1) + d_old(:,1))*(r_new(:,1) - r_old(:,1)) + &
                (d_new(:,2) + d_old(:,2))*(r_new(:,2) - r_old(:,2)) + &
                (d_new(:,3) + d_old(:,3))*(r_new(:,3) - r_old(:,3))

      if (nelec == 1) then
         argexpo = 0.5d0 * (d2_new(1) - d2_old(1))*dt + prod(1)
      else
         argexpo = 0.5d0 * sum((d2_new(:) - d2_old(:))*dt + prod(:))
      endif
      
      q = psi_new / psi_old
      q = dexp(-argexpo) * q*q

      call random_number(u)

      if (u <= q) then
         n_accep = n_accep + 1_8
         r_old   = r_new
         d_old   = d_new
         d2_old  = d2_new
         psi_old = psi_new
      end if

   end do

   ! Compute energy and acceptance ratio 
   energy = energy / dble(nmax)
   accep  = dble(n_accep) / dble(nmax)

end subroutine vmc

!****************************************************************************

subroutine pdmc(a,dt,nmax,energy,accep,tau,E_ref,natoms,geom,nelec,z)

   ! This subroutine implements the Pure Diffusion Monte Carlo (PDMC) method to
   ! estimate the average energy

   implicit none

   ! Input parameters:
   ! a      : parameter of wave function
   ! dt     : step size
   ! tau    : imaginary time variable
   ! geom   : geometry of system (xyz coordinates)
   ! nmax   : Monte Carlo steps
   ! natoms : number of atoms in the system
   ! nelec  : number of electrons in the system
   ! z      : atomic numbers of the atoms in the system
   ! E_ref  : energy shift

   ! Output parameters:
   ! energy : estimated average energy of the system
   ! accep  : acceptance ratio of the Monte Carlo step simualtion

   double precision, intent(in)  :: a, dt, tau, E_ref, geom(natoms,3)
   integer*8       , intent(in)  :: nmax 
   integer         , intent(in)  :: natoms, nelec, z(natoms)
   double precision, intent(out) :: energy, accep

   integer*8                            :: istep, n_accep
   integer                              :: i
   double precision                     :: sq_dt, u, psi_old, psi_new, argexpo, q
   double precision                     :: e, w, normalization, tau_current
   double precision, dimension(nelec,3) :: r_old, r_new, d_old, d_new, chi
   double precision, dimension(nelec)   :: prod, d2_old, d2_new

   sq_dt = dsqrt(dt)

   ! Initialization
   energy  = 0.d0
   n_accep = 0_8
   normalization = 0.d0

   w           = 1.d0
   tau_current = 0.d0

   ! Initialize electron positions using a Gaussian distribution
   call random_gauss(r_old,3,nelec)

   ! Compute the initial drift vector based on the current positions
   call drift(a,r_old,d_old,nelec,natoms,geom)  

   ! Compute the sum of the square of drift vectors
   d2_old(:)  =  d_old(:,1)*d_old(:,1) + &
                  d_old(:,2)*d_old(:,2) + &
                  d_old(:,3)*d_old(:,3)

   ! Calculate the initial wave function value
   psi_old = psi(a,r_old,natoms,geom,nelec)

   ! Main Monte Carlo loop
   do istep = 1, nmax

      ! Evaluate local energy
      e = e_loc(a,r_old,natoms,geom,nelec,z)

      ! Update the weight factor using the difference from the reference energy
      w = w * dexp(-dt*(e - E_ref))

      ! Accumulate weighted energy and weight for normalization
      normalization = normalization + w
      energy = energy + w*e

      ! Increment current tau value
      tau_current = tau_current + dt

      ! Reset when tau is reached
      if (tau_current > tau) then
         w           = 1.d0
         tau_current = 0.d0
      endif

      ! Propose a new configuration by moving electrons
      call random_gauss(chi,3,nelec)
      do i = 1, nelec
         r_new(i,:) = r_old(i,:) + dt*d_old(i,:) + chi(i,:)*sq_dt
      enddo

      ! Recalculate drift vector for the new positions
      call drift(a,r_new,d_new,nelec,natoms,geom) 

      ! Compute the square of new drift vector (d_new)
      d2_new(:) = d_new(:,1)*d_new(:,1) + &
                  d_new(:,2)*d_new(:,2) + &
                  d_new(:,3)*d_new(:,3)

      ! Evaluate wave function at the new positions
      psi_new = psi(a,r_new,natoms,geom,nelec)

      ! Metropolis acceptance criterion calculation
      prod(:) = (d_new(:,1) + d_old(:,1))*(r_new(:,1) - r_old(:,1)) + &
                (d_new(:,2) + d_old(:,2))*(r_new(:,2) - r_old(:,2)) + &
                (d_new(:,3) + d_old(:,3))*(r_new(:,3) - r_old(:,3))

      if (nelec == 1) then
         argexpo = 0.5d0 * (d2_new(1) - d2_old(1))*dt + prod(1)
      else
         argexpo = 0.5d0 * sum((d2_new(:) - d2_old(:))*dt + prod(:))
      endif

      q = psi_new / psi_old
      q = dexp(-argexpo) * q*q

      call random_number(u)

      if (u <= q) then
         n_accep = n_accep + 1_8
         r_old   = r_new
         d_old   = d_new
         d2_old  = d2_new
         psi_old = psi_new
      end if

   end do

   ! Compute energy and acceptance ratio 
   energy = energy / normalization
   accep  = dble(n_accep) / dble(nmax)

end subroutine pdmc

!****************************************************************************

subroutine read_input(method,a,dt,E_ref,tau,nruns,nmax,natoms,atom,geom,charge)

   ! This subroutine read input parameters from the input file.

   implicit none

   character(len=4)              :: method
   double precision, intent(out) :: a, dt, E_ref, tau
   integer*8       , intent(out) :: nmax
   integer         , intent(out) :: nruns, natoms, charge
   integer                       :: i

   character(len=2), allocatable, intent(out) :: atom(:)
   double precision, allocatable, intent(out) :: geom(:,:)

   ! Open input file 

   open(1, file="qmc.in", status="old")

   ! Read input parameters

   read(1,*) 
   read(1,*) natoms

   allocate(atom(natoms), geom(natoms,3))

   read(1,*) charge

   do i = 1, natoms
      read(1,*) atom(i), geom(i,1), geom(i,2), geom(i,3)
   enddo

   read(1,*) 
   read(1,*)
   read(1,*) method
   read(1,*)
   read(1,*)
   read(1,*) a
   read(1,*) dt
   read(1,*) nmax
   read(1,*) nruns

   if (method == "pdmc") then
      read(1,*) E_ref
      read(1,*) tau
   endif

   ! Close input file

   close(1)

end subroutine read_input

!****************************************************************************

subroutine output(method,a,dt,E_ref,tau,nruns,nmax,natoms,atom,geom,z,charge, &
                  ave_e,ave_a,err_e,err_a)

   ! This subroutine is designed to print and write to a file the input variables, 
   ! and the results of the Quantum Monte Carlo simulation.

   implicit none

   character(len=4), intent(in) :: method
   character(len=2), intent(in) :: atom(natoms)
   double precision, intent(in) :: a, dt, E_ref, tau, geom(natoms,3), &
                                   ave_e, ave_a, err_e, err_a
   integer*8                    :: nmax
   integer                      :: nruns, i, natoms, charge, z(natoms)

   ! Print input variables

   print *, " "
   print *, "--------- System -----------"
   print *, " "
   print *, "- Molecular Charge: ", charge
   print *, " "
   print *, "Atom, Atomic Number (Z),  Cartesian Coordinates (Å)"
   do i = 1, natoms
      print *, atom(i), z(i), geom(i,:)
   enddo
   print *, " "
   print *, "-- Simulation Parameters -- "
   print *, " "
   print *, "- Method: ", method
   print *, "- Wavefunction Parameter (a): ", a
   print *, "- Step Size (dt): ", dt
   if (method == "pdmc") then
      print *, "- Energy Shift (E_ref): ", E_ref
      print *, "- Imaginary Time Variable (tau): ", tau
   endif
   print *, "- Number of MC Steps: ", nmax
   print *, "- Number of MC Runs: ", nruns
   print *, " "

   ! Write input variables to the summary output

   write(2,"(A)") "==================== Simulation Summary ===================="
   write(2,"(A)") " "
   write(2,"(A)") "----------------------------"
   write(2,"(A)") "          System            "
   write(2,"(A)") "----------------------------"
   write(2,"(A)") " "
   write(2,"(A,I2)") "- Molecular Charge: ", charge
   write(2,"(A)") " "
   write(2,"(A)") "Atom, Atomic Number (Z),  Cartesian Coordinates (Å)"
   do i = 1, natoms
      write(2,"(A2,1X,I2,3(1X,F9.6))") atom(i), z(i), geom(i,:)
   enddo
   write(2,"(A)") " "
   write(2,"(A)") " "
   write(2,"(A)") "----------------------------"
   write(2,"(A)") "         Simulation         "
   write(2,"(A)") "----------------------------"
   write(2,"(A)") " "
   write(2,"(A,A)") "- Method: ", method
   write(2,"(A,F7.3)") "- Wavefunction Parameter (a): ", a
   write(2,"(A,F6.3)") "- Step Size (dt): ", dt
   if (method == "pdmc") then
      write(2,"(A,F9.6)") "- Energy Shift (E_ref): ", E_ref
      write(2,"(A,F12.6)") "- Imaginary Time Variable (tau): ", tau
   endif
   write(2,"(A,I10)") "- Number of MC Steps: ", nmax
   write(2,"(A,I10)") "- Number of MC Runs: ", nruns
   write(2,"(A)") " "

   ! Print results

   print *, "--------- Results ---------"
   print *, " "
   print *, 'E = ', ave_e, '+/-', err_e
   print *, 'A = ', ave_a, '+/-', err_a
   print *, " "

   ! Write results to the summary output

   write(2,"(A)") " "
   write(2,"(A)") "----------------------------"
   write(2,"(A)") "          Results           "
   write(2,"(A)") "----------------------------"
   write(2,"(A)") " "
   write(2,"(A,F9.6,A,E7.1)") "Energy (Eh): ", ave_e, " +/- ", err_e
   write(2,"(A,F9.6,A,E7.1)") "Acceptance Ratio: ", ave_a, " +/- ", err_a

end subroutine output

!****************************************************************************

subroutine find_atomic_number(input_element,nelem,z)

   ! This subroutine finds the atomic numbers for a list of chemical elements.
   ! It searches each element in the input list within a predefined array of all
   ! 118 chemical elements ordered by atomic number and assigns the corresponding
   ! atomic number to the output array 'z'.

   ! input_element : array of chemical element symbols to find atomic numbers for
   ! nelem         : number of elements in the 'input_element' array
   ! z             : output array of atomic numbers corresponding to 'input_element'
   ! elements      : ordered array of all 118 element symbols by atomic number

   implicit none

   character(len=2), intent(in)     :: input_element(nelem)
   integer         , intent(in)     :: nelem
   integer         , intent(out)    :: z(nelem)
   character(len=2), dimension(118) :: elements
   integer                          :: i, j

   ! Initializing the array of elements with the symbols of the 118 chemical elements
   elements = (/ &
         "H ", "He", &
         "Li", "Be", &
         "B ", "C ", "N ", "O ", "F ", "Ne", &
         "Na", "Mg", &
         "Al", "Si", "P ", "S ", "Cl", "Ar", &
         "K ", "Ca", &
         "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
         "Ga", "Ge", "As", "Se", "Br", "Kr", &
         "Rb", "Sr", &
         "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", &
         "In", "Sn", "Sb", "Te", "I ", "Xe", &
         "Cs", "Ba", &
         "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
         "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
         "Tl", "Pb", "Bi", "Po", "At", "Rn", &
         "Fr", "Ra", &
         "Ac", "Th", "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", &
         "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", &
         "Nh", "Fl", "Mc", "Lv", "Ts", "Og"/)

   do j = 1, nelem

      z(j) = 0  ! Initialize the atomic number to 0 to indicate that the element has not been found

      ! Search through the 'elements' array to find each input element and assign its atomic number
      do i = 1, 118
         if (trim(adjustl(elements(i))) == trim(adjustl(input_element(j)))) then
            z(j) = i
            exit
         endif
      enddo 

      if (z(j) <= 0) then
         print *, "Error: element not found."
      elseif (z(j) > 2) then
         print *, "Warning: QMC program is specific for H, He, H2, H2+, and H3+ systems"
      endif
      
   enddo

end subroutine find_atomic_number

!****************************************************************************

end module qmc_routines
