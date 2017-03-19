#include<GPE_dataStructures.h>

void allocateSO3Matrices(System_Parameters& sysparams) {
	int i1,i2;
	sysparams.Fx=new dcomplex* [3];
	sysparams.Fy=new dcomplex* [3];
	sysparams.Fz=new dcomplex* [3];
	for(i1=0;i1<3;i1++){
		(sysparams.Fx)[i1]=new dcomplex[3];
		(sysparams.Fy)[i1]=new dcomplex[3];
		(sysparams.Fz)[i1]=new dcomplex[3];
		for(i2=0;i2<3;i2++) {
			sysparams.Fx[i1][i2]=0.0;
			sysparams.Fy[i1][i2]=0.0;
			sysparams.Fz[i1][i2]=0.0;
		};
	};
	sysparams.Fx[1][2].im=-1.0;
	sysparams.Fx[2][1].im=1.0;
	sysparams.Fy[2][0].im=-1.0;
	sysparams.Fy[0][2].im=1.0;
	sysparams.Fz[0][1].im=-1.0;
	sysparams.Fz[1][0].im=1.0;
}

void freeSO3Matrices(System_Parameters sysparams) {
	int i;
	for(i=0;i<3;i++) {
		delete[] sysparams.Fx[i];
		delete[] sysparams.Fy[i];
		delete[] sysparams.Fz[i];
	}
	delete[] sysparams.Fx;
	delete[] sysparams.Fy;
	delete[] sysparams.Fz;
}

void allocate_GPE_Parameters(GPE_Parameters *gpe_ptr,System_Parameters sysparams) {
	register int i1,i2,i3,i4;
	(*gpe_ptr).EXP=new (nothrow) double[sysparams.M];
	(*gpe_ptr).SeffLDA    =new (nothrow) double ** [(sysparams.SqrtMinusDmu_Steps)];
	(*gpe_ptr).dSeff_dPsi2=new (nothrow) double ** [(sysparams.SqrtMinusDmu_Steps)];
	(*gpe_ptr).dSeff_dF2  =new (nothrow) double ** [(sysparams.SqrtMinusDmu_Steps)];
	(*gpe_ptr).n=new (nothrow) dcomplex **** [(sysparams.SqrtMinusDmu_Steps)];
	(*gpe_ptr).phi_avg=new (nothrow) dcomplex *** [(sysparams.SqrtMinusDmu_Steps)];
	for(i1=0;i1<sysparams.SqrtMinusDmu_Steps;i1++) {
		((*gpe_ptr).SeffLDA    )[i1]=new (nothrow) double * [(sysparams.AbsPsi_Steps)];
		((*gpe_ptr).dSeff_dPsi2)[i1]=new (nothrow) double * [(sysparams.AbsPsi_Steps)];
		((*gpe_ptr).dSeff_dF2  )[i1]=new (nothrow) double * [(sysparams.AbsPsi_Steps)];
		((*gpe_ptr).n)[i1]=new (nothrow) dcomplex *** [(sysparams.AbsPsi_Steps)];
		((*gpe_ptr).phi_avg)[i1]=new (nothrow) dcomplex ** [(sysparams.AbsPsi_Steps)];
		for(i2=0;i2<sysparams.AbsPsi_Steps;i2++) {
			((*gpe_ptr).SeffLDA    )[i1][i2]=new (nothrow) double [(sysparams.AbsF_Steps)];
			((*gpe_ptr).dSeff_dPsi2)[i1][i2]=new (nothrow) double [(sysparams.AbsF_Steps)];
			((*gpe_ptr).dSeff_dF2  )[i1][i2]=new (nothrow) double [(sysparams.AbsF_Steps)];
			((*gpe_ptr).n)[i1][i2]=new (nothrow) dcomplex** [(sysparams.AbsPsi_Steps)];
			((*gpe_ptr).phi_avg)[i1][i2]=new (nothrow) dcomplex* [(sysparams.AbsPsi_Steps)];
			for(i3=0;i3<sysparams.AbsPsi_Steps;i3++) {
				((*gpe_ptr).n)[i1][i2][i3]=new (nothrow) dcomplex* [3];
				((*gpe_ptr).phi_avg)[i1][i2][i3]=new (nothrow) dcomplex [3];
				for(i4=0;i4<3;i4++)
					((*gpe_ptr).n)[i1][i2][i3][i4]=new (nothrow) dcomplex [3];
			};
		};
	};
}

void allocate_GPE_Psis(GPE_Psis *Psis_ptr,System_Parameters sysparams) {
	(*Psis_ptr).Psi1 =(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sysparams.M);
	(*Psis_ptr).PsiI1=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sysparams.M);
	(*Psis_ptr).PsiK1=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sysparams.M);
	(*Psis_ptr).Psi2 =(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sysparams.M);
	(*Psis_ptr).PsiI2=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sysparams.M);
	(*Psis_ptr).PsiK2=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sysparams.M);
	(*Psis_ptr).Psi3 =(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sysparams.M);
	(*Psis_ptr).PsiI3=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sysparams.M);
	(*Psis_ptr).PsiK3=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*sysparams.M);
}

void create_FFTW_plans(FFTW_Struct *FFTW_Variables_ptr,GPE_Psis Psis,System_Parameters sysparams) {
	register int j;
	for(j=0;j<sysparams.M;j++) {
		Psis.PsiK1[j][0]=Psis.Psi1[j][0];
		Psis.PsiK1[j][1]=Psis.Psi1[j][1];
	}; // PsiK1 will serve as the temporary array with the use of which one can define the plans.
	fftw_plan_with_nthreads(sysparams.nthreads);
	(*FFTW_Variables_ptr).plan       =fftw_plan_dft_3d(sysparams.Mx,sysparams.My,sysparams.Mz,
							   Psis.PsiK1,Psis.PsiK1,FFTW_FORWARD,FFTW_MEASURE);
	(*FFTW_Variables_ptr).inverseplan=fftw_plan_dft_3d(sysparams.Mx,sysparams.My,sysparams.Mz,
							   Psis.PsiK1,Psis.PsiK1,FFTW_BACKWARD,FFTW_MEASURE);
	for(j=0;j<sysparams.M;j++) {
		Psis.PsiK1[j][0]=0.0;
		Psis.PsiK1[j][1]=0.0;
	};
};

void free_SystemParameters(System_Parameters sysparams) {
	int j;
	delete[] sysparams.saving_times;
	delete[] sysparams.addNoiseTimes;
	delete[] sysparams.noOfFreqs;
	for(j=0;j<sysparams.length_addNoiseTimes_array;j++) {
		delete[] sysparams.randomnessRate[j];
		delete[] sysparams.maxFreqRate[j];
	};
	delete[] sysparams.randomnessRate;
	delete[] sysparams.maxFreqRate;
}

void free_GPE_Parameters(GPE_Parameters gpe,System_Parameters sysparams) {
	register int i1,i2,i3,i4;
	delete[] gpe.EXP;
	for(i1=0;i1<sysparams.SqrtMinusDmu_Steps;i1++) {
		for(i2=0;i2<sysparams.AbsPsi_Steps;i2++) {
			for(i3=0;i3<sysparams.AbsPsi_Steps;i3++) {
				for(i4=0;i4<3;i4++)
					delete[] gpe.n[i1][i2][i3][i4];
				delete[] gpe.n[i1][i2][i3];
				delete[] gpe.phi_avg[i1][i2][i3];
			};
			delete[] gpe.SeffLDA    [i1][i2];
			delete[] gpe.dSeff_dPsi2[i1][i2];
			delete[] gpe.dSeff_dF2  [i1][i2];
			delete[] gpe.n[i1][i2];
			delete[] gpe.phi_avg[i1][i2];
		};
		delete[] gpe.SeffLDA    [i1];
		delete[] gpe.dSeff_dPsi2[i1];
		delete[] gpe.dSeff_dF2  [i1];
		delete[] gpe.n[i1];
		delete[] gpe.phi_avg[i1];
	};
	delete[] gpe.SeffLDA;
	delete[] gpe.dSeff_dPsi2;
	delete[] gpe.dSeff_dF2;
	delete[] gpe.n;
	delete[] gpe.phi_avg;
}

void free_GPE_Psis(GPE_Psis Psis) {
	fftw_free(Psis.Psi1 );
	fftw_free(Psis.PsiI1);
	fftw_free(Psis.PsiK1);
	fftw_free(Psis.Psi2 );
	fftw_free(Psis.PsiI2);
	fftw_free(Psis.PsiK2);
	fftw_free(Psis.Psi3 );
	fftw_free(Psis.PsiI3);
	fftw_free(Psis.PsiK3);
}

void destroy_FFTW_plans(FFTW_Struct FFTW_Variables) {
	// Destroy FFTW plans
	fftw_destroy_plan(FFTW_Variables.plan);
	fftw_destroy_plan(FFTW_Variables.inverseplan);
};

