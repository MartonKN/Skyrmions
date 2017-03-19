void monitoring(const char filenameMonitoring[],GPE_Psis Psis,GPE_Parameters gpe,
		System_Parameters sysparams,FFTW_Struct FFTW_Variables);
// Displays energy, norms and angular momenta and also saves them to file.

void kineticEnergy(double *kinEnergy1,double *kinEnergy2,double *kinEnergy3,GPE_Psis Psis,
	           System_Parameters sysparams,FFTW_Struct FFTW_Variables);
// kinEnergy1=\Sum_{ij} a^2*\beta/(6^2*J) Psi1_i^* \Laplacian_{ij} Psi1_j.
// kinEnergy2=\Sum_{ij} a^2*\beta/(6^2*J) Psi2_i^* \Laplacian_{ij} Psi2_j.

void potentialEnergy(double *potEnergy,GPE_Psis Psis,System_Parameters sysparams,GPE_Parameters gpe);
// potEnergy=\Sum_i Seff(\mu1_i,\mu2_i,Psi1_i,Psi2_i)

void angularMomenta(double *Lx,double *Ly,double *Lz,GPE_Psis Psis,System_Parameters sysparams);
// L = \Sum Psi1_i^* (r x grad)_i Psi1_i

void norms(double *normPsi1,double *normPsi2,double *normPsi3,GPE_Psis Psis,System_Parameters sysparams);
// normPsi1=sqrt(\Sum_i |Psi1_i|^2)
// normPsi2=sqrt(\Sum_i |Psi2_i|^2)
// normPsi3=sqrt(\Sum_i |Psi3_i|^2)

dcomplex PontrjaginIndex(GPE_Psis Psis,System_Parameters sysparams);
// Calculates the Pontrjagin index of the configuration.

dcomplex mult_dcomplex(dcomplex a,dcomplex b);
// Multiplies two dcomplex numbers.

dcomplex PontrAuxiliaryFunction(double epsilon[3][3][3],dcomplex derivatives[3][3]);
// Auxiliary routine used only by function PontrjaginIndex().
// It was needed basically because the compiler was slow in compiling the 9 times embedded 'for' loops.