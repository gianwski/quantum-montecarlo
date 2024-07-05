![Header image](./images/chemical-banner_short.jpg)

# Quantum Monte Carlo (QMC) Program
Quantum Monte Carlo (QMC) program, written in Fortran 90, is designed to estimate the average energy and its associated error of system through Monte Carlo simulations. It is specifically developed to accommodate five key cases: $H$, $He$, $H_2$, $H_2^+$, and $H_3^+$. It supports both Variational Monte Carlo (VMC) and Pure Diffusion Monte Carlo (PDMC) methods, allowing for flexible analysis of quantum systems. This program reads input parameters and system configurations, performs Monte Carlo simulations, and outputs the calculated average energy of the system and the acceptance ratio with the respective error.
The program is actually a combination of two distinct codes. The first code _qmc.f90_ serves as the main program, while the second _qmc_routines.f90_ is a module that comprises a set of functions and subroutines required for the calculation.
All quantities must be given in _atomic units_ (except for geometry in ångström). The results will be in _atomic units_.

## Prerequisites
Before compiling the QMC program, ensure you have a Fortran compiler installed on your system. The program has been tested with _gfortran_, which is part of the GNU Compiler Collection (GCC). Other Fortran compilers may also work.

For detailed instructions on installing _gfortran_, please refer to this [guide](https://fortran-lang.org/learn/os_setup/install_gfortran/).  
Further information about _gfortran_ can be found on its [official web page](https://gcc.gnu.org/fortran/).

## Components

### _Main Program (qmc)_
The main program, _qmc_, performs quantum Monte Carlo simulations. It reads input parameters, including simulation method, system parameters, and atomic geometry. After converting units and assigning atomic numbers, it executes the simulation using either the VMC or PDMC method, based on a input keyword. It calculates the average energy and acceptance ratio with their associated errors and writes the results to the outputs.

### _QMC Routines Module (qmc_routines)_
This module contains routines required for the QMC simulations, including:

- ***psi*** : Calculates the value of the wave function for a specific system configuration.
- ***potential*** : Computes the local potential energy of the system considering nucleus-electron, electron-electron and nucleus-nucleus interactions.
- ***psi_1e_1N*** : Calculates the contribution to the wavefunction of one electron and one nucleus ($ e^{-a \cdot | r_e - R_N|} $).
- ***kinetic*** : Calculates the local kinetic energy of the system.
- ***e_loc*** : Computes the local energy for a given configuration of the system.
- ***ave_error*** : Calculates the average and statistical error of an input array.
- ***random_gauss*** : Generates normally distributed random numbers using the Box-Muller transform.
- ***drift*** : Computes the drift vector.
- ***vmc*** : Perform the VMC method.
- ***pdmc*** : Perform the PDMC method.
- ***read_input*** : Reads simulation parameters from the input file.
- ***output*** : Print and write to a summary output file the input parameters and the results of the simulation.
- ***find_atomic_number*** : Assign atomic numbers to a given list of chemical elements.

## Instructions
1. Create the input file named _qmc.in_. Refer to the documentation for input configuration.
2. Compile the program using the following commands:

   ```bash
   gfortran -c qmc_routines.f90
   gfortran qmc_routines.o qmc.f90 -o qmc
   ```
3. Run the executable:

   ```bash
   ./qmc
   ```

4. A summary output file named _qmc.out_ will be generated, containing the input variables and results of the QMC simulation.

## Test Program (_qmc_test_)
The _qmc_test_ program is developed following the principles of _Test-Driven Development_ (TDD). It includes a comprehensive suite of automated tests aimed at validating the functionality of various QMC routines. This program tests several components essential for QMC simulations, such as:

- psi functions
- potential energy calculations
- kinetic energy calculations
- atomic number identification

This program is an essential tool for maintaining the integrity and reliability of the QMC simulation software by systematically verifying each component's correct functionality.

### Usage
First, compile the _qmc_test.f90_ program:
```bash
gfortran -o qmc_test qmc_routines.f90 qmc_test.f90
 ```
To run the compiled program and execute the tests, use:
```bash
./qmc_test
 ```

## Credits
If you use the QMC program in your research or any publication, please cite the author(s).

## License
Distributed under the GPL-3.0 License. See _LICENSE_ for more information.

## Contact
Gianluca Regni - [_gianluca.regni@studenti.unipg.it_](mailto:gianluca.regni@studenti.unipg.it)
