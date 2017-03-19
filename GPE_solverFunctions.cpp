#include<GPE_dataStructures.h>
#include<GPE_solverFunctions.h>
#include<GPE_spline.h>

void iteration(GPE_Psis Psis,GPE_Parameters *gpe,System_Parameters sysparams,FFTW_Struct FFTW_Variables) {
	// Performs a single dt imaginary timestep iteration of the GPE using the RK4IP method.
	copying(Psis.PsiK1,Psis.PsiK2,Psis.PsiK3, Psis.Psi1,Psis.Psi2,Psis.Psi3,sysparams);
	diffusion_exponential(Psis.Psi1,*gpe,sysparams,FFTW_Variables);
	diffusion_exponential(Psis.Psi2,*gpe,sysparams,FFTW_Variables);
	diffusion_exponential(Psis.Psi3,*gpe,sysparams,FFTW_Variables);
	copying(Psis.PsiI1,Psis.PsiI2,Psis.PsiI3,Psis.Psi1,Psis.Psi2,Psis.Psi3,sysparams);
	nondiffusion(Psis.PsiK1,Psis.PsiK2,Psis.PsiK3,*gpe,sysparams);
	diffusion_exponential(Psis.PsiK1,*gpe,sysparams,FFTW_Variables);
	diffusion_exponential(Psis.PsiK2,*gpe,sysparams,FFTW_Variables);
	diffusion_exponential(Psis.PsiK3,*gpe,sysparams,FFTW_Variables);
	weighted_sum(Psis.Psi1,Psis.Psi2,Psis.Psi3,1.0,Psis.PsiK1,Psis.PsiK2,Psis.PsiK3,1.0/6.0,sysparams);
	weighted_sum(Psis.PsiK1,Psis.PsiK2,Psis.PsiK3,0.5,Psis.PsiI1,Psis.PsiI2,Psis.PsiI3,1.0,sysparams);
	(*gpe).t+=0.5*sysparams.dt;
	nondiffusion(Psis.PsiK1,Psis.PsiK2,Psis.PsiK3,*gpe,sysparams);
	weighted_sum(Psis.Psi1,Psis.Psi2,Psis.Psi3,1.0,Psis.PsiK1,Psis.PsiK2,Psis.PsiK3,1.0/3.0,sysparams);
	weighted_sum(Psis.PsiK1,Psis.PsiK2,Psis.PsiK3,0.5,Psis.PsiI1,Psis.PsiI2,Psis.PsiI3,1.0,sysparams);
	nondiffusion(Psis.PsiK1,Psis.PsiK2,Psis.PsiK3,*gpe,sysparams);
	weighted_sum(Psis.Psi1,Psis.Psi2,Psis.Psi3,1.0,Psis.PsiK1,Psis.PsiK2,Psis.PsiK3,1.0/3.0,sysparams);
	weighted_sum(Psis.PsiK1,Psis.PsiK2,Psis.PsiK3,1.0,Psis.PsiI1,Psis.PsiI2,Psis.PsiI3,1.0,sysparams);
	diffusion_exponential(Psis.PsiK1,*gpe,sysparams,FFTW_Variables);
	diffusion_exponential(Psis.PsiK2,*gpe,sysparams,FFTW_Variables);
	diffusion_exponential(Psis.PsiK3,*gpe,sysparams,FFTW_Variables);
	diffusion_exponential(Psis.Psi1,*gpe,sysparams,FFTW_Variables);
	diffusion_exponential(Psis.Psi2,*gpe,sysparams,FFTW_Variables);
	diffusion_exponential(Psis.Psi3,*gpe,sysparams,FFTW_Variables);
	(*gpe).t+=0.5*sysparams.dt;
	nondiffusion(Psis.PsiK1,Psis.PsiK2,Psis.PsiK3,*gpe,sysparams);
	weighted_sum(Psis.Psi1,Psis.Psi2,Psis.Psi3,1.0,Psis.PsiK1,Psis.PsiK2,Psis.PsiK3,1.0/6.0,sysparams);
}

void copying(fftw_complex* Psi1_final,  fftw_complex* Psi2_final,  fftw_complex* Psi3_final,
	     fftw_complex* Psi1_initial,fftw_complex* Psi2_initial,fftw_complex* Psi3_initial,
	     System_Parameters sysparams) {
	// Psi1_final <-- Psi1_initial, Psi2_final <-- Psi2_initial, Psi3_final <-- Psi3_initial.
	int threadID;
	register int j;
	#pragma omp parallel private(threadID) 
	{
		threadID=omp_get_thread_num();
		cpu_set_t new_mask;
		cpu_set_t was_mask;
		CPU_ZERO(&new_mask);		// Makes the new_mask CPU set empty
		CPU_SET(threadID,&new_mask);	// Orders new_mask to threadID
		if (sched_getaffinity(0,sizeof(was_mask),&was_mask) == -1) {
			cerr<<"Error: sched_getaffinity("<<threadID<<", sizeof(was_mask), &was_mask)"<<endl;
		}
		if (sched_setaffinity(0, sizeof(new_mask), &new_mask) == -1) {
			cerr<<"Error: sched_setaffinity("<<threadID<<", sizeof(was_mask), &was_mask)"<<endl;
		}
		#pragma omp for schedule(dynamic) nowait private(j)
		for(j=0;j<sysparams.M;j++) {
			Psi1_final[j][0]=Psi1_initial[j][0];
			Psi1_final[j][1]=Psi1_initial[j][1];
		};
		#pragma omp for schedule(dynamic) nowait private(j)
		for(j=0;j<sysparams.M;j++) {
			Psi2_final[j][0]=Psi2_initial[j][0];
			Psi2_final[j][1]=Psi2_initial[j][1];
		};
		#pragma omp for schedule(dynamic) nowait private(j)
		for(j=0;j<sysparams.M;j++) {
			Psi3_final[j][0]=Psi3_initial[j][0];
			Psi3_final[j][1]=Psi3_initial[j][1];
		};
	}
}

void nondiffusion(fftw_complex* Psi1,fftw_complex* Psi2,fftw_complex* Psi3,GPE_Parameters gpe,System_Parameters sysparams) {
	// (Psi) <-- N(Psi), where Psi=(Psi1,Psi2,Psi3) and N(Psi)_i=-sysparams.dt*d(SeffNonKinetic)_d(Psi_i^*)
	int threadID;
	register double SqrtMinusDmu;
	register dcomplex PsiVec[3];
	register int ix,iy,iz,j;
	register double diffix,diffiy,diffiz;
	register dcomplex tmp;
	#pragma omp parallel private(threadID) 
	{
		threadID=omp_get_thread_num();
		cpu_set_t new_mask;
		cpu_set_t was_mask;
		CPU_ZERO(&new_mask);		// Makes the new_mask CPU set empty
		CPU_SET(threadID,&new_mask);	// Orders new_mask to threadID
		if (sched_getaffinity(0,sizeof(was_mask),&was_mask) == -1) {
			cerr<<"Error: sched_getaffinity("<<threadID<<", sizeof(was_mask), &was_mask)"<<endl;
		}
		if (sched_setaffinity(0, sizeof(new_mask), &new_mask) == -1) {
			cerr<<"Error: sched_setaffinity("<<threadID<<", sizeof(was_mask), &was_mask)"<<endl;
		}
		#pragma omp for schedule(dynamic) nowait private(ix,iy,iz,j,diffix,diffiy,diffiz,PsiVec,SqrtMinusDmu,tmp)
		for(ix=0;ix<sysparams.Mx;ix++) {
			diffix=ix-0.5*(sysparams.Mx-1.0);
			for(iy=0;iy<sysparams.My;iy++) {
				diffiy=iy-0.5*(sysparams.My-1.0);
				for(iz=0;iz<sysparams.Mz;iz++) {
					diffiz=iz-0.5*(sysparams.Mz-1.0);
					SqrtMinusDmu=sqrt(sysparams.trapx*diffix*diffix
					                 +sysparams.trapy*diffiy*diffiy
					                 +sysparams.trapz*diffiz*diffiz);
					j=iz + sysparams.Mz*(iy + sysparams.My*ix);
					PsiVec[0]=Psi1[j];
					PsiVec[1]=Psi2[j];
					PsiVec[2]=Psi3[j];
					tmp=evaluateSpline(1,SqrtMinusDmu,PsiVec,gpe,sysparams);
					Psi1[j][0]=-sysparams.dt*tmp.re;
					Psi1[j][1]=-sysparams.dt*tmp.im;
					tmp=evaluateSpline(2,SqrtMinusDmu,PsiVec,gpe,sysparams);
					Psi2[j][0]=-sysparams.dt*tmp.re;
					Psi2[j][1]=-sysparams.dt*tmp.im;
					tmp=evaluateSpline(3,SqrtMinusDmu,PsiVec,gpe,sysparams);
					Psi3[j][0]=-sysparams.dt*tmp.re;
					Psi3[j][1]=-sysparams.dt*tmp.im;
				};
			};
		};
	}
}

void weighted_sum(fftw_complex *PsiA1,fftw_complex *PsiA2,fftw_complex *PsiA3,double A,
		  fftw_complex *PsiB1,fftw_complex *PsiB2,fftw_complex *PsiB3,double B,
		  System_Parameters sysparams) {
	// PsiA <-- A*PsiA + B*PsiB, where PsiA=(PsiA1,PsiA2,PsiA3), PsiB=(PsiB1,PsiB2,PsiB3).
	int threadID;
	register int j;
	#pragma omp parallel private(threadID) 
	{
		threadID=omp_get_thread_num();
		cpu_set_t new_mask;
		cpu_set_t was_mask;
		CPU_ZERO(&new_mask);		// Makes the new_mask CPU set empty
		CPU_SET(threadID,&new_mask);	// Orders new_mask to threadID
		if (sched_getaffinity(0,sizeof(was_mask),&was_mask) == -1) {
			cerr<<"Error: sched_getaffinity("<<threadID<<", sizeof(was_mask), &was_mask)"<<endl;
		}
		if (sched_setaffinity(0, sizeof(new_mask), &new_mask) == -1) {
			cerr<<"Error: sched_setaffinity("<<threadID<<", sizeof(was_mask), &was_mask)"<<endl;
		}
		#pragma omp for schedule(dynamic) nowait private(j)
		for(j=0;j<sysparams.M;j++) {
			PsiA1[j][0]=A*PsiA1[j][0] + B*PsiB1[j][0];
			PsiA1[j][1]=A*PsiA1[j][1] + B*PsiB1[j][1];
		};
		#pragma omp for schedule(dynamic) nowait private(j)
		for(j=0;j<sysparams.M;j++) {
			PsiA2[j][0]=A*PsiA2[j][0] + B*PsiB2[j][0];
			PsiA2[j][1]=A*PsiA2[j][1] + B*PsiB2[j][1];
		};
		#pragma omp for schedule(dynamic) nowait private(j)
		for(j=0;j<sysparams.M;j++) {
			PsiA3[j][0]=A*PsiA3[j][0] + B*PsiB3[j][0];
			PsiA3[j][1]=A*PsiA3[j][1] + B*PsiB3[j][1];
		};
	}
}

void diffusion_exponential(fftw_complex* Psi,GPE_Parameters gpe,System_Parameters sysparams,FFTW_Struct FFTW_Variables) {
	// Calculates exp(-dt LAPLACE/(2 MASS)) Psi.
	int threadID;
	register int j;
	fftw_execute_dft(FFTW_Variables.plan,Psi,Psi);
	#pragma omp parallel private(threadID) 
	{
		threadID=omp_get_thread_num();
		cpu_set_t new_mask;
		cpu_set_t was_mask;
		CPU_ZERO(&new_mask);		// Makes the new_mask CPU set empty
		CPU_SET(threadID,&new_mask);	// Orders new_mask to threadID
		if (sched_getaffinity(0,sizeof(was_mask),&was_mask) == -1) {
			cerr<<"Error: sched_getaffinity("<<threadID<<", sizeof(was_mask), &was_mask)"<<endl;
		}
		if (sched_setaffinity(0, sizeof(new_mask), &new_mask) == -1) {
			cerr<<"Error: sched_setaffinity("<<threadID<<", sizeof(was_mask), &was_mask)"<<endl;
		}
		#pragma omp for schedule(dynamic) nowait private(j)
		for(j=0;j<sysparams.M;j++) {
			Psi[j][0]*=gpe.EXP[j];
			Psi[j][1]*=gpe.EXP[j];
		};
	}
	fftw_execute_dft(FFTW_Variables.inverseplan,Psi,Psi);
}