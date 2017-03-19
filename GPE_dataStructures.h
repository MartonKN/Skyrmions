// This header file contains the definitions of the data structures that are necessary for the imaginary time
// GPE solver. It has to be included in each .cpp file in this program.


// *********************************************************************************************************


// NECESSARY INCLUDE FILES
// FFTW libraries
  #include<fftw3.h>
  // #include<fftw_threads.h>

// C libraries
  #include<string.h>
  #include<stdio.h>

// C++ libraries
  #include<iostream>
  #include<string>
  #include<cmath>
  #include<fstream>
  #include<cstdlib>
  using namespace std;
  
// OpenMP
#ifndef __USE_GNU
	#define __USE_GNU
#endif
#include<omp.h>
#include<sched.h>

// Whether to save densities
#define SAVE_DENSITIES 1

// *********************************************************************************************************


#ifndef GPE_DATA_STRUCTURES
#define GPE_DATA_STRUCTURES


// *********************************************************************************************************

#include "dcomplex.h"
#include "matrix.h"

// SYSTEM PARAMETERS
struct System_Parameters {
	int Mx,My,Mz;                     // number of points in each direction of the grid
	int M;                            // M=Mx*My*Mz
	double ax,ay,az;                  // length of one spatial step
	double tmax,dt;                   // maximal imaginary time and imaginary time step of the simulation
	
	// All quantities of energy dimension are measured in units 
	// of the average interaction strength V0=(V11+V22)/2:=1.
	double T;                         // dimensionless temperature
	double V0,V2;                     // dimensionless interaction terms
	double J;                         // dimensionless parameters of the hopping
	double mu_0;                      // dimensionless chemical potential in the middle of the trap
	                                  // -mu[i][j][k] = -mu_0 + trapx*(i-(Mx-1)/2)^2 + trapy*(j-(My-1)/2)^2 + trapz*(k-(Mz-1)/2)^2,
	double trapx,trapy,trapz;         // trapx=m*omega_x^2*ax^2/2
	
	// Parameters for the splines
	int Nmax;                         // the maximal number of atoms per site
	double AbsPsi_Max;                // the edges of the cube on which the spline is defined
	double SqrtMinusDmu_Max;           
	int AbsPsi_Steps,AbsF_Steps;      // number of sample points for the spline in each direction
	int SqrtMinusDmu_Steps;
	double D_AbsPsi,D_AbsF;           // stepsizes of the grid
	double D_SqrtMinusDmu;            // D_AbsPsi1=AbsPsi1_Max/(AbsPsi1_Steps-1.0), etc.
	
	// Representation of SO(3)
	dcomplex **Fx,**Fy,**Fz;
	
	// Interaction part of the Hamiltonian and matrices of the particle numbers
	dcomplex **HInt;
	dcomplex ****N;                  // N[j1][j2][][] is the operator phi^+[j1] phi[j2]
	dcomplex ***phi;                 // phi[j][][] is the annihilation operator phi[j]
	int sizeOfHamiltonian;           // the size of the Hamiltonian matrix: 0, ..., (sizeHamiltonian-1)
	
	// Add random noise
	double *addNoiseTimes;
	int length_addNoiseTimes_array;
	int *noOfFreqs;
	double **maxFreqRate;
	double **randomnessRate;
	long idum; // GPE_ran seed. It is set in function load_SystemParameters().
	// Algorithm:
	// Psis.Psi_k += amplitude*exp(i*(freqx*jx+freqy*jy+freqz*jz));  (k=0,1,2)
	// noOfFreqs: the number of times this random part will be added at randomly chosen quasimomenta (freqx,freqy,freqz).
	// maxFreqRate[3]: the random frequencies will come from an interval width endpoints at 
	//                 +-maxFreqRate[k]*(the maximal possible quasimomentum)
	// randomnessRate: amplitude=randomnessRate[k]*|Psis.Psi_k|_2
	
	// Parameters for saving
	double *saving_times;             
	int length_saving_times_array;    
	int dispform_saving_times[2];     // saving times will appear in the filenames in 
	                                  // %dispform_saving_times[0].dispform_saving_times[1]f format
	int nthreads;                     // number of threads used during parallelization.
};


// *********************************************************************************************************


// GROSS-PITAEVSKII EQUATION (GPE) PARAMETERS
struct GPE_Parameters {
	// The arrays of size sysparams.M will be allocated in the program.
	// The arrays are stored in "row-major format": A[jx][jy][jz] --> A[jz + Mz*(jy + My*jx)].
	double *EXP;
	// This array contains the terms appearing in the exponentialized diffusion exponent:
	// exp(-dt LAPLACE/(2 MASS)) = exp((1-cos(kx ax)) dt/(MASS ax^2)) exp((1-cos(ky ay)) dt/(MASS ay^2)) exp((1-cos(kz az)) dt/(MASS az^2)) 
	//                           ~ exp(dt kx^2/(2 MASS)) exp(dt ky^2/(2 MASS)) exp(dt kz^2/(2 MASS))
	// Plus they also contain the necesary 1/(sysparams.Mx*sysparams.My*sysparams.Mz) 
	// normalization of the Fourier transform.
	double t; // dimensionless imaginary time.
	
	                                          // Arrays of size SqrtMinusDmuSteps*Psi1Steps*Psi2Steps.
	                                          // The arrays are stored in "row-major format": 
	                                          // A[jSqrtMinusDmu_2][jPsi1_2][jPsi2_2] --> A[jPsi2_2 + Psi2_2Steps*(jPsi1_2 + Psi1_2Steps*jSqrtMinusDmu_2)].
	double ***SeffLDA;                        // SeffLDA: the effective action without the kinetic term
	double ***dSeff_dPsi2,***dSeff_dF2;       // dSeffLDA/d(|Psi|^2) and dSeffLDA/d((Psi* F Psi)^2).
	dcomplex *****n;                          // Average local densities. (<phi^+[i1] phi[i2]> = n[j_SqrtMinusDmu][j_alpha_p][j_alpha_m][i1][i2].)
	dcomplex ****phi_avg;                     // Average of the annihilation operators (<phi[i]> = phi_avg[j_SqrtMinusDmu][j_alpha_p][j_alpha_m][i].)
};


// *********************************************************************************************************


// ORDER PARAMETERS
struct GPE_Psis {
	fftw_complex *Psi1, *PsiI1, *PsiK1, *Psi2, *PsiI2, *PsiK2, *Psi3, *PsiI3, *PsiK3;
	// The arrays of size sysparams.M will be allocated in the program.
	// The arrays are stored in "row-major format": A[jx][jy][jz] --> A[jz + Mz*(jy + My*jx)].
	// Psi is the order parameter in real space, PsiI is that in interaction picture and PsiK is
	// a temporary variable.
};


// *********************************************************************************************************


// VARIABLES OF THE FFTW FAST FOURIER TRANSFORM
struct FFTW_Struct{
	fftw_plan plan;
	fftw_plan inverseplan;
};


// *********************************************************************************************************


// COMPLEX NUMBER TYPE:
//   We use the built-in fftw_complex type of FFTW, a double[2] array.
//   fftw_complex z; z[0]=real part of z, z[1]=imaginary part.


// *********************************************************************************************************


// OTHERS
#define sqr_norm(x) ((x)[0]*(x)[0] + (x)[1]*(x)[1])
// square of an fftw_complex number.
#define __min(x,y) ((x)<(y) ? (x) : (y))
#define __max(x,y) ((x)>(y) ? (x) : (y))

#define USE_MAX_NUM_PROCS -1
// Do not modify these two constants!

#define PI 3.1415926535897932384626433

// *********************************************************************************************************


#endif
