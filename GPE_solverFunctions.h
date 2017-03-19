// These functions are the heart of the GPE solver. 
// 'iteration()' performs one iteration step of imaginary time sysparams.dt of the GPE.
// It uses only the functions listed below.

void iteration(GPE_Psis Psis,GPE_Parameters *gpe,System_Parameters sysparams,FFTW_Struct FFTW_Variables);
// Performs a single dt imaginary timestep iteration of the GPE using the RK4IP method.

void copying(fftw_complex* Psi1_final,fftw_complex* Psi2_final,fftw_complex* Psi3_final,
             fftw_complex* Psi1_initial,fftw_complex* Psi2_initial,fftw_complex* Psi3_initial,
             System_Parameters sysparams);
// Psi1_final <-- Psi1_initial, Psi2_final <-- Psi2_initial.
// NOTE: Parallelized function.

void nondiffusion(fftw_complex* Psi1,fftw_complex* Psi2,fftw_complex* Psi3,GPE_Parameters gpe,System_Parameters sysparams);
// (Psi) <-- N(Psi), where Psi=(Psi1,Psi2)
// NOTE: Parallelized function.

void weighted_sum(fftw_complex *PsiA1,fftw_complex *PsiA2,fftw_complex *PsiA3,double A,
                  fftw_complex *PsiB1,fftw_complex *PsiB2,fftw_complex *PsiB3,double B,
                  System_Parameters sysparams);
// PsiA <-- A*PsiA + B*PsiB, where PsiA=(PsiA1,PsiA2), PsiB=(PsiB1,PsiB2).
// NOTE: Parallelized function.

void diffusion_exponential(fftw_complex* Psi,GPE_Parameters gpe,System_Parameters sysparams,FFTW_Struct FFTW_Variables);
// Psi <-- exp(-dt LAPLACE/(2 MASS)) Psi.
// NOTE: Parallelized function.
