#include<GPE_dataStructures.h>
#include<GPE_IO.h>
#include<GPE_allocations.h>
#include<GPE_solverFunctions.h>
#include<GPE_spline.h>
#include<GPE_monitoring.h>

#define SAVING 1
// Whether to save data during run


int main(int argc, char *argv[]) {
	System_Parameters sysparams;
	GPE_Parameters gpeparams;
	GPE_Psis Psis;
	FFTW_Struct FFTW_Variables;
	
	if(argc!=2) {
		cerr<<"The name of the file, that contains the other necessary filenames, has to be given. "<<endl;
		cerr<<"(Incorrect number of parameters in function main().)"<<endl;
		return 0;
	}
	char filenameSystemParameters[200],filenameBaseOrderParameters[200];
	if(!create_filenames(argv[1],filenameSystemParameters,filenameBaseOrderParameters))
		return 0;
	
	// Initialization of FFTW threading
	if(!fftw_init_threads()){
		cerr<<"Problem occured calling function fftw_init_threads()."<<endl;
		return 0;
	}
	// Load system parameters
	if(!load_SystemParameters(filenameSystemParameters,sysparams)) {
		cerr<<"Problem occured calling function load_SystemParameters(). "<<endl;
		return 0;
	}
	// OpenMP
	if(sysparams.nthreads==USE_MAX_NUM_PROCS || sysparams.nthreads>omp_get_num_procs()) {
		sysparams.nthreads=omp_get_num_procs();
	}
	omp_set_num_threads(sysparams.nthreads);
	// Allocate memory
	allocate_GPE_Parameters(&gpeparams,sysparams);
	allocate_GPE_Psis(&Psis,sysparams);
	allocateSO3Matrices(sysparams);
	createHIntAndNs(sysparams);
	cout<<"Allocations done."<<endl;
	
	
	// Initialize splines
	resetSplineParameters(sysparams);
	initializeSpline(gpeparams,sysparams);
	cout<<"Spline initializations done."<<endl;
	// Initialize order parameters and gpeparams
	initializePsisAndGPE(gpeparams,sysparams,Psis);
	cout<<"Initialization of Psis done."<<endl;
	// Set the FFTW plans and set t=0.
	create_FFTW_plans(&FFTW_Variables,Psis,sysparams);
	gpeparams.t=0.0;
	// Reset time step size
	reset_dt(gpeparams,sysparams,Psis,FFTW_Variables);
	saveSystemParameters(filenameSystemParameters,sysparams);
	cout<<"System parameters have been reset."<<endl;
	
	// Iterations
	cout<<"Iterations start."<<endl;
	do{
		addRandomNoise(Psis,gpeparams,sysparams);
		#if SAVING
			if(save_Psis(filenameBaseOrderParameters,Psis,gpeparams,sysparams,FFTW_Variables)==0) return 0;
		#endif
		iteration(Psis,&gpeparams,sysparams,FFTW_Variables);
	}while(gpeparams.t<=sysparams.tmax);
	cout<<endl;
	cout<<"Iterations are over."<<endl;
	
	// Free allocated memory
	free_GPE_Parameters(gpeparams,sysparams);
	free_GPE_Psis(Psis);
	destroy_FFTW_plans(FFTW_Variables);
	freeHIntAndNs(sysparams);
	freeSO3Matrices(sysparams);
	return 1;
}

