dcomplex evaluateSpline(int Option,double SqrtMinusDmu,dcomplex PsiVec[3],GPE_Parameters gpe,System_Parameters sysparams);
// Using the data stored in 'gpe' it calculates the following splines at points (Dmu_2,PsiVec):
// Option==0: SeffLDA, Option==1: dSeff_d(Psi1*), Option==2: dSeff_d(Psi2*), Option==3: dSeff_d(Psi3*)

cmatrix evaluateSplineDensities(double SqrtMinusDmu,dcomplex PsiVec[3],GPE_Parameters gpe,System_Parameters sysparams);
// Calculates densities using spline parameters gpe.n.
// First it transforms Psi onto the subspace spanned by alpha_p*v_p+alpha_m*v_m, where v_p is the +1, v_m is the -1 eigenvector
// of the matrix Fz. alpha_p, alpha_m>0. Then evaluates the spline value for the transformed matrix of densities and transforms it back.

cmatrix evaluateSplinePhi(double SqrtMinusDmu,dcomplex PsiVec[3],GPE_Parameters gpe,System_Parameters sysparams);
// Calculates the expectation values of the annihilation operators using spline parameters gpe.phi_avg.
// First it transforms Psi onto the subspace spanned by alpha_p*v_p+alpha_m*v_m, where v_p is the +1, v_m is the -1 eigenvector
// of the matrix Fz. alpha_p, alpha_m>0. Then evaluates the spline value for the transformed phi_average and transforms it back.

void initializeSpline(GPE_Parameters gpe,System_Parameters sysparams);
// Fills up the spline data arrays in 'gpe':
// gpe.SeffLDA:     the non-kinetic part of effective action on each site
// gpe.dSeff_dPsi2: derivative of SeffLDA with respect to |Psi|^2.
// gpe.dSeff_dF2:   derivative of SeffLDA with respect to <F>^2.
// Note: <F> in [0,1], ie. the average is taken by a Psi of norm 1.
// NOTE: Parallelized function.

void calcAbsPsiAndAbsF(dcomplex PsiVec[3],double& AbsPsi,double& AbsF,double FVec[3],System_Parameters sysparams);
// Determines AbsPsi and AbsF for a given PsiVec[] configuration.

void setPsi(double AbsPsi,double AbsF,dcomplex PsiVec[3]);
// Small routine used only by functions initializeSpline() and resetSplineParameters().
// Sets a Psi configuration with PsiVec'*PsiVec=AbsPsi^2 
// and ((PsiVec'*Fx*PsiVec)^2+(PsiVec'*Fy*PsiVec)^2+(PsiVec'*Fz*PsiVec)^2)/(PsiVec'*PsiVec)^2=AbsF^2

void transformPsi(dcomplex PsiVec[3],double& alpha_p,double& alpha_m,dmatrix& U);
// Calculates the appropriate SO(3) transformation that brings PsiVec into 
// the subspace spanned by v_p and v_m, the +1 and -1 eigenvectors of Fz.
// U*PsiVec=(alpha_p*v_p+alpha_m*v_m)*exp(i*tau), where alpha_p and alpha_m
// are positive numbers.

double calculateSeff(double SqrtMinusDmu,dcomplex PsiVec[3],System_Parameters sysparams);
// Calculates the non-kinetic part of Seff (SeffLDA).
// The inclusion of the 3J terms in the redefinition of the chemical potentials is a matter of taste here.
// Either you choose (-mu+3J)(N1+N2+N3), or just -mu(N1+N2+N3). These can be set in the program.

cmatrix calculateDensitiesAndPhi(int Option,double SqrtMinusDmu,dcomplex PsiVec[3],System_Parameters sysparams);
// Option==0: densities
// Option==1: phi_average
// Calculates the matrix of densities/the vector containing the average values of the annihilation operator.
// The inclusion of the 3J terms in the redefinition of the chemical potentials is a matter of taste here.
// Either you choose (-mu+3J)(N1+N2+N3), or just -mu(N1+N2+N3). These can be set in the program.

void phiAvereages(fftw_complex* phiAverage1,fftw_complex* phiAverage2,fftw_complex* phiAverage3,
                  GPE_Psis Psis,GPE_Parameters gpe,System_Parameters sysparams);
// Calculates the expectation values of the annihilation operators phi1, phi2 and phi3,
// and puts the results in phiAverage1, phiAverage2 etc.
// These arrays shall be of the size sysparams.M.

void densityCorrelations(fftw_complex* n11,fftw_complex* n12,fftw_complex* n13,fftw_complex* n22,fftw_complex* n23,fftw_complex* n33,
                         GPE_Psis Psis,GPE_Parameters gpe,System_Parameters sysparams,FFTW_Struct FFTW_Variables);
// This function calculates the Fourier transform of the density correlations and puts them into n11, n12 etc.
// The arrays n11, n22 and n33 shall be of the size sysparams.M.
// It is a good idea to put them into Psis.PsiI1 etc. or Psis.PsiK1 etc., since these are only auxiliary arrays.

void resetSplineParameters(System_Parameters& sysparams);
// Sets the boundaries of the spline to appropriate values.
// Finds the minimum of Seff in terms of AbsPsi and AbsF, and 
// sets sysparams.AbsPsi_Max to be twice the value of this.
// It also sets sysparams.SqrtMinusDmu_Max to be the maximal 
// possible value at the edge of the cloud.
// NOTE: Parallelized function.

void createHIntAndNs(System_Parameters& sysparams);
// Generates the interaction part of the Hamiltonian (sysparams.HInt) containing the terms with V0 and V2.
// (Chemical potentials and Psis are not included.)
// sysparams.N1, sysparams.N2, sysparams.N3 are also filled up.
// Note: sysparams.HInt, N1, N2 and N3 are allocated here, therefore they have to be freed 
// using freeHIntAndNs().

void multMatrices(dcomplex** A,dcomplex** B,dcomplex** C,int sizeMx);
void MxAdjoint(dcomplex** A,dcomplex** A_dagger,int sizeMx);

void freeHIntAndNs(System_Parameters sysparams);
// Frees array sysparams.HInt, sysparams.N1, sysparams.N2 and sysparams.N3.

void annihilate(int flavor,int n[3],double& coefficient,System_Parameters sysparams);
// Annihilation operator.

void create(int flavor,int n[3],double& coefficient,System_Parameters sysparams);
// Creation operator.

int index(int n1,int n2,int n3,int& n,System_Parameters sysparams);
// (n1,n2,n3) --> n
// Converts the filling numbers (n1,n2,n3) to indices of the Hamiltonian.

int inv_index(int& n1,int& n2,int& n3,int n,System_Parameters sysparams);
// n --> (n1,n2,n3)
// The inverse of function index(). Tells what filling numbers correspond
// to a given row of the Hamiltonian,

int diagonalize(dcomplex **HermitianMx,double* eigenvalues,int sizeMx);
// Diagonalizes HermitianMx. The array HermitianMx will be overwritten by the eigenvectors.
// HamiltonianMx=U*eigenvalues*U'
// HamiltonianMx <-- U