# Skyrmions

Skyrmions are topological excitations, that were originally proposed by Skyrme to describe particles (hadrons). In the paper Scientific Reports 5, 7692 (2015), we proposed to create skyrmion excitations in a quantum emulator of ultracold atoms that are in the superfluid state and have an order parameter of non-trivial topology. The system is described in detail in the text.

How to use the code:
The code consists of three parts: 
(1) Initialization using a MATLAB routine, where the system parameters can be set up (temperature, number of lattice sites etc.)
This routine generates three files: 
- [...]_Filenames: contains the name of some files that the C++ simulation (2)
- [...]_SystemParameters: contains all parameters of the system, that are read by the 
- [...]_RunFile: runfile for running the simulation on the cluster
The beginning of the files '[...]' is set up by this routine, and contains 

(2) C++ solver that calculates the ground state of the skyrmion in the topologically non-trivial state. This can be compiled using the makefile that is provided with the code, producing the executable 'imagTimeGPESolver'. 

The simulation can be run as './imagTimeGPESolver [...]_SystemParameters'
