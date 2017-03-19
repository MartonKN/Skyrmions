void allocateSO3Matrices(System_Parameters& sysparams);
// Allocates and fills up Fx, Fy and Fz.
void freeSO3Matrices(System_Parameters sysparams);
// Frees Fx, Fy and Fz.
void allocate_GPE_Parameters(GPE_Parameters *gpe_ptr,System_Parameters sysparams);
// Allocates the array gpe.EXP of size sizeof(double)*sysparams.M 
// and the arrays gpe.SeffLDA, gpe.dSeff_dPsi2 and dSeff_dF2 for the spline.
void allocate_GPE_Psis(GPE_Psis *Psis_ptr,System_Parameters sysparams);
// Allocates the arrays of fftw_complex numbers of size 2*sizeof(double)*sysparams.M
// for the variables Psis1 and Psis2.
void create_FFTW_plans(FFTW_Struct *FFTW_Variables_ptr,GPE_Psis Psis,System_Parameters sysparams);
// Before using any of the FFTW functions, so called plans have to be created. 
// These tell the functions how to optimize their performance on the given computer.
void free_SystemParameters(System_Parameters sysparams);
// This short function solely frees the array 'sysparam.saving_times' allocated by
// the function 'load_SystemParameters()' in GPE_IO.cpp.
void free_GPE_Parameters(GPE_Parameters gpe,System_Parameters sysparams);
// Frees arrays allocated by function 'allocate_GPE_Parameters()'.
void free_GPE_Psis(GPE_Psis Psis);
// Frees arrays allocated by function 'allocate_GPE_Psis()'.
void destroy_FFTW_plans(FFTW_Struct FFTW_Variables);
// After using them, the FFTW plans have to be destroyed using functions supported
// by the FFTW library.
