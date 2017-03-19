#include<GPE_dataStructures.h>
#include<GPE_spline.h>

#define SPLINE_ERR 0
dcomplex evaluateSpline(int Option,double SqrtMinusDmu,dcomplex PsiVec[3],GPE_Parameters gpe,System_Parameters sysparams) {
	// Option==0: SeffLDA, Option==1: dSeff_d(Psi1*), Option==2: dSeff_d(Psi2*), Option==3: dSeff_d(Psi3*)
	dcomplex result;
	double S=0.0,dS_dPsi2=0.0,dS_dF2=0.0;
	double AbsPsi,AbsF,FVec[3]; // FVec=(<Fx>,<Fy>,<Fz>)
	dcomplex tmp_dcomplex;
	int j_AbsPsi,j_AbsF,j_SqrtMinusDmu;
	double tmp_AbsPsi,tmp_AbsF,tmp_SqrtMinusDmu;
	double u_AbsPsi,u_AbsF,u_SqrtMinusDmu;
	int i2;
	
	// Calculate AbsPsi, FVec and AbsF
	calcAbsPsiAndAbsF(PsiVec,AbsPsi,AbsF,FVec,sysparams);
	
	// Fit spline
	tmp_AbsPsi=__min(AbsPsi/sysparams.D_AbsPsi,2147483646.0);
	tmp_AbsF  =__min(AbsF  /sysparams.D_AbsF  ,2147483646.0);
	tmp_SqrtMinusDmu=__min(SqrtMinusDmu/sysparams.D_SqrtMinusDmu,2147483646.0);
	j_AbsPsi=floor(tmp_AbsPsi);
	j_AbsF  =floor(tmp_AbsF);
	j_SqrtMinusDmu=floor(tmp_SqrtMinusDmu);
	u_AbsPsi=tmp_AbsPsi-j_AbsPsi;
	u_AbsF  =tmp_AbsF  -j_AbsF;
	u_SqrtMinusDmu=tmp_SqrtMinusDmu-j_SqrtMinusDmu;
	
	// What if the value of AbsPsi is out of range.
	if(j_AbsPsi>=(sysparams.AbsPsi_Steps-1)) {
		#if SPLINE_ERR
			cerr<<"AbsPsi value out of range in function evaluateSpline()."<<endl;
		#endif
		u_AbsPsi=1.0;
		j_AbsPsi=sysparams.AbsPsi_Steps-2;
		if(Option==0) {
			S+=fabs((AbsPsi-sysparams.AbsPsi_Max)*
			        (gpe.SeffLDA[__min(j_SqrtMinusDmu,sysparams.SqrtMinusDmu_Steps-1)][sysparams.AbsPsi_Steps-1][__min(j_AbsF,sysparams.AbsF_Steps-1)]
			        -gpe.SeffLDA[__min(j_SqrtMinusDmu,sysparams.SqrtMinusDmu_Steps-1)][sysparams.AbsPsi_Steps-2][__min(j_AbsF,sysparams.AbsF_Steps-1)])
			        /sysparams.D_AbsPsi
			       );
		}
		if(Option==1 || Option==2 || Option==3) {
			dS_dPsi2+=fabs((AbsPsi-sysparams.AbsPsi_Max)*
			               (gpe.dSeff_dPsi2[__min(j_SqrtMinusDmu,sysparams.SqrtMinusDmu_Steps-1)][sysparams.AbsPsi_Steps-1][__min(j_AbsF,sysparams.AbsF_Steps-1)]
			               -gpe.dSeff_dPsi2[__min(j_SqrtMinusDmu,sysparams.SqrtMinusDmu_Steps-1)][sysparams.AbsPsi_Steps-2][__min(j_AbsF,sysparams.AbsF_Steps-1)])
			               /sysparams.D_AbsPsi
			              );
			dS_dF2+=fabs((AbsPsi-sysparams.AbsPsi_Max)*
			             (gpe.dSeff_dF2[__min(j_SqrtMinusDmu,sysparams.SqrtMinusDmu_Steps-1)][sysparams.AbsPsi_Steps-1][__min(j_AbsF,sysparams.AbsF_Steps-1)]
			             -gpe.dSeff_dF2[__min(j_SqrtMinusDmu,sysparams.SqrtMinusDmu_Steps-1)][sysparams.AbsPsi_Steps-2][__min(j_AbsF,sysparams.AbsF_Steps-1)])
			             /sysparams.D_AbsPsi
			            );
		}
	}
	// What if AbsF is out of range
	if(j_AbsF>=(sysparams.AbsF_Steps-1)) {
		#if SPLINE_ERR
			cerr<<"AbsF value is out of range in function evaluateSpline()."<<endl;
		#endif
		u_AbsF=1.0;
		j_AbsF=sysparams.AbsF_Steps-2;
	}
	// What if Dmu is out of range
	if(j_SqrtMinusDmu>=(sysparams.SqrtMinusDmu_Steps-1)) {
		#if SPLINE_ERR
			cerr<<"mu value is out of range in function evaluateSpline()."<<endl;
		#endif
		u_SqrtMinusDmu=1.0;
		j_SqrtMinusDmu=sysparams.SqrtMinusDmu_Steps-2;
	}
	
	// Evaluate spline
	if(Option==0) {
		result=S;
		result.re+=(((1.0-u_SqrtMinusDmu)*gpe.SeffLDA[j_SqrtMinusDmu][j_AbsPsi  ][j_AbsF  ]+u_SqrtMinusDmu*gpe.SeffLDA[j_SqrtMinusDmu+1][j_AbsPsi  ][j_AbsF  ])*(1.0-u_AbsPsi)
		           +((1.0-u_SqrtMinusDmu)*gpe.SeffLDA[j_SqrtMinusDmu][j_AbsPsi+1][j_AbsF  ]+u_SqrtMinusDmu*gpe.SeffLDA[j_SqrtMinusDmu+1][j_AbsPsi+1][j_AbsF  ])*u_AbsPsi      )*(1.0-u_AbsF)
		          +(((1.0-u_SqrtMinusDmu)*gpe.SeffLDA[j_SqrtMinusDmu][j_AbsPsi  ][j_AbsF+1]+u_SqrtMinusDmu*gpe.SeffLDA[j_SqrtMinusDmu+1][j_AbsPsi  ][j_AbsF+1])*(1.0-u_AbsPsi)
		           +((1.0-u_SqrtMinusDmu)*gpe.SeffLDA[j_SqrtMinusDmu][j_AbsPsi+1][j_AbsF+1]+u_SqrtMinusDmu*gpe.SeffLDA[j_SqrtMinusDmu+1][j_AbsPsi+1][j_AbsF+1])*u_AbsPsi      )*u_AbsF;
		return result;
	}
	else if(Option==1 || Option==2 || Option==3) {
		if(AbsPsi>1e-8) {
			dS_dPsi2+=(((1.0-u_SqrtMinusDmu)*gpe.dSeff_dPsi2[j_SqrtMinusDmu][j_AbsPsi  ][j_AbsF  ]+u_SqrtMinusDmu*gpe.dSeff_dPsi2[j_SqrtMinusDmu+1][j_AbsPsi  ][j_AbsF  ])*(1.0-u_AbsPsi)
			          +((1.0-u_SqrtMinusDmu)*gpe.dSeff_dPsi2[j_SqrtMinusDmu][j_AbsPsi+1][j_AbsF  ]+u_SqrtMinusDmu*gpe.dSeff_dPsi2[j_SqrtMinusDmu+1][j_AbsPsi+1][j_AbsF  ])*u_AbsPsi      )*(1.0-u_AbsF)
			         +(((1.0-u_SqrtMinusDmu)*gpe.dSeff_dPsi2[j_SqrtMinusDmu][j_AbsPsi  ][j_AbsF+1]+u_SqrtMinusDmu*gpe.dSeff_dPsi2[j_SqrtMinusDmu+1][j_AbsPsi  ][j_AbsF+1])*(1.0-u_AbsPsi)
			          +((1.0-u_SqrtMinusDmu)*gpe.dSeff_dPsi2[j_SqrtMinusDmu][j_AbsPsi+1][j_AbsF+1]+u_SqrtMinusDmu*gpe.dSeff_dPsi2[j_SqrtMinusDmu+1][j_AbsPsi+1][j_AbsF+1])*u_AbsPsi      )*u_AbsF;
			dS_dF2+=(((1.0-u_SqrtMinusDmu)*gpe.dSeff_dF2[j_SqrtMinusDmu][j_AbsPsi  ][j_AbsF  ]+u_SqrtMinusDmu*gpe.dSeff_dF2[j_SqrtMinusDmu+1][j_AbsPsi  ][j_AbsF  ])*(1.0-u_AbsPsi)
			        +((1.0-u_SqrtMinusDmu)*gpe.dSeff_dF2[j_SqrtMinusDmu][j_AbsPsi+1][j_AbsF  ]+u_SqrtMinusDmu*gpe.dSeff_dF2[j_SqrtMinusDmu+1][j_AbsPsi+1][j_AbsF  ])*u_AbsPsi      )*(1.0-u_AbsF)
			       +(((1.0-u_SqrtMinusDmu)*gpe.dSeff_dF2[j_SqrtMinusDmu][j_AbsPsi  ][j_AbsF+1]+u_SqrtMinusDmu*gpe.dSeff_dF2[j_SqrtMinusDmu+1][j_AbsPsi  ][j_AbsF+1])*(1.0-u_AbsPsi)
			        +((1.0-u_SqrtMinusDmu)*gpe.dSeff_dF2[j_SqrtMinusDmu][j_AbsPsi+1][j_AbsF+1]+u_SqrtMinusDmu*gpe.dSeff_dF2[j_SqrtMinusDmu+1][j_AbsPsi+1][j_AbsF+1])*u_AbsPsi      )*u_AbsF;
			tmp_dcomplex=0.0;
			for(i2=0;i2<3;i2++) {
				tmp_dcomplex.re+=2.0*((sysparams.Fx[Option-1][i2].re*PsiVec[i2].re-sysparams.Fx[Option-1][i2].im*PsiVec[i2].im)*FVec[0]
				                     +(sysparams.Fy[Option-1][i2].re*PsiVec[i2].re-sysparams.Fy[Option-1][i2].im*PsiVec[i2].im)*FVec[1]
				                     +(sysparams.Fz[Option-1][i2].re*PsiVec[i2].re-sysparams.Fz[Option-1][i2].im*PsiVec[i2].im)*FVec[2]);
				tmp_dcomplex.im+=2.0*((sysparams.Fx[Option-1][i2].re*PsiVec[i2].im+sysparams.Fx[Option-1][i2].im*PsiVec[i2].re)*FVec[0]
				                     +(sysparams.Fy[Option-1][i2].re*PsiVec[i2].im+sysparams.Fy[Option-1][i2].im*PsiVec[i2].re)*FVec[1]
				                     +(sysparams.Fz[Option-1][i2].re*PsiVec[i2].im+sysparams.Fz[Option-1][i2].im*PsiVec[i2].re)*FVec[2]);
			};
			result.re=dS_dPsi2*PsiVec[Option-1].re+dS_dF2*(tmp_dcomplex.re/(AbsPsi*AbsPsi)-2.0*AbsF*AbsF/(AbsPsi*AbsPsi)*PsiVec[Option-1].re);
			result.im=dS_dPsi2*PsiVec[Option-1].im+dS_dF2*(tmp_dcomplex.im/(AbsPsi*AbsPsi)-2.0*AbsF*AbsF/(AbsPsi*AbsPsi)*PsiVec[Option-1].im);
			return result;
		}
		else {
			result=0.0;
			return result;
		}
	}
	else {
		cerr<<"Invalid option in function evaluateSpline()."<<endl;
		result=0.0;
		return result;
	}
	result=0.0;
	return result;
}

cmatrix evaluateSplineDensities(double SqrtMinusDmu,dcomplex PsiVec[3],GPE_Parameters gpe,System_Parameters sysparams) {
	cmatrix n(3,3);
	double alpha_p,alpha_m;
	register int j_alpha_p,j_alpha_m,j_SqrtMinusDmu;
	double tmp_alpha_p,tmp_alpha_m,tmp_SqrtMinusDmu;
	register double u_alpha_p,u_alpha_m,u_SqrtMinusDmu;
	register int i1,i2;
	dmatrix U(3,3);
	
	#if SAVE_DENSITIES
		// Dmu
		tmp_SqrtMinusDmu=__min(SqrtMinusDmu/sysparams.D_SqrtMinusDmu,2147483646.0);
		j_SqrtMinusDmu=floor(tmp_SqrtMinusDmu);
		u_SqrtMinusDmu=tmp_SqrtMinusDmu-j_SqrtMinusDmu;
		// What if Dmu is out of range
		if(j_SqrtMinusDmu>=(sysparams.SqrtMinusDmu_Steps-1)) {
			#if SPLINE_ERR
				cerr<<"mu value is out of range in function evaluateSplineDensities()."<<endl;
			#endif
			u_SqrtMinusDmu=1.0;
			j_SqrtMinusDmu=sysparams.SqrtMinusDmu_Steps-2;
		}
		
		// alpha_p,alpha_m
		transformPsi(PsiVec,alpha_p,alpha_m,U);
		tmp_alpha_p=__min(alpha_p/sysparams.D_AbsPsi,2147483646.0);
		tmp_alpha_m=__min(alpha_m/sysparams.D_AbsPsi,2147483646.0);
		j_alpha_p=floor(tmp_alpha_p);
		j_alpha_m=floor(tmp_alpha_m);
		u_alpha_p=tmp_alpha_p-j_alpha_p;
		u_alpha_m=tmp_alpha_m-j_alpha_m;
		// What if alpha_p is out of range
		if(j_alpha_p>=(sysparams.AbsPsi_Steps-1)) {
			#if SPLINE_ERR
				cerr<<"alpha_p value is out of range in function evaluateSplineDensities()."<<endl;
			#endif
			u_alpha_p=1.0;
			j_alpha_p=sysparams.AbsPsi_Steps-2;
		}
		// What if alpha_m is out of range
		if(j_alpha_m>=(sysparams.AbsPsi_Steps-1)) {
			#if SPLINE_ERR
				cerr<<"alpha_m value is out of range in function evaluateSplineDensities()."<<endl;
			#endif
			u_alpha_m=1.0;
			j_alpha_m=sysparams.AbsPsi_Steps-2;
		}
		// Fit spline to the transformed matrix of densities
		for(i1=0;i1<3;i1++) {
			for(i2=0;i2<3;i2++) {
				n[i1][i2]=(((1.0-u_SqrtMinusDmu)*gpe.n[j_SqrtMinusDmu][j_alpha_p  ][j_alpha_m  ][i1][i2]+u_SqrtMinusDmu*gpe.n[j_SqrtMinusDmu+1][j_alpha_p  ][j_alpha_m  ][i1][i2])*(1.0-u_alpha_p)
				          +((1.0-u_SqrtMinusDmu)*gpe.n[j_SqrtMinusDmu][j_alpha_p+1][j_alpha_m  ][i1][i2]+u_SqrtMinusDmu*gpe.n[j_SqrtMinusDmu+1][j_alpha_p+1][j_alpha_m  ][i1][i2])*u_alpha_p      )*(1.0-u_alpha_m)
				         +(((1.0-u_SqrtMinusDmu)*gpe.n[j_SqrtMinusDmu][j_alpha_p  ][j_alpha_m+1][i1][i2]+u_SqrtMinusDmu*gpe.n[j_SqrtMinusDmu+1][j_alpha_p  ][j_alpha_m+1][i1][i2])*(1.0-u_alpha_p)
				          +((1.0-u_SqrtMinusDmu)*gpe.n[j_SqrtMinusDmu][j_alpha_p+1][j_alpha_m+1][i1][i2]+u_SqrtMinusDmu*gpe.n[j_SqrtMinusDmu+1][j_alpha_p+1][j_alpha_m+1][i1][i2])*u_alpha_p      )*u_alpha_m;
			};
		};
		n*=U;
		U.transpose();
		n=U*n;
		return n;
	#else
		n*=0.0;
		return n;
	#endif
}

cmatrix evaluateSplinePhi(double SqrtMinusDmu,dcomplex PsiVec[3],GPE_Parameters gpe,System_Parameters sysparams) {
	cmatrix phi_average(3,1);
	double alpha_p,alpha_m;
	register int j_alpha_p,j_alpha_m,j_SqrtMinusDmu;
	double tmp_alpha_p,tmp_alpha_m,tmp_SqrtMinusDmu;
	register double u_alpha_p,u_alpha_m,u_SqrtMinusDmu;
	register int j;
	dmatrix U(3,3);
	
	// Dmu
	tmp_SqrtMinusDmu=__min(SqrtMinusDmu/sysparams.D_SqrtMinusDmu,2147483646.0);
	j_SqrtMinusDmu=floor(tmp_SqrtMinusDmu);
	u_SqrtMinusDmu=tmp_SqrtMinusDmu-j_SqrtMinusDmu;
	// What if Dmu is out of range
	if(j_SqrtMinusDmu>=(sysparams.SqrtMinusDmu_Steps-1)) {
		#if SPLINE_ERR
			cerr<<"mu value is out of range in function evaluateSplinePhi()."<<endl;
		#endif
		u_SqrtMinusDmu=1.0;
		j_SqrtMinusDmu=sysparams.SqrtMinusDmu_Steps-2;
	}
	
	// alpha_p,alpha_m
	transformPsi(PsiVec,alpha_p,alpha_m,U);
	tmp_alpha_p=__min(alpha_p/sysparams.D_AbsPsi,2147483646.0);
	tmp_alpha_m=__min(alpha_m/sysparams.D_AbsPsi,2147483646.0);
	j_alpha_p=floor(tmp_alpha_p);
	j_alpha_m=floor(tmp_alpha_m);
	u_alpha_p=tmp_alpha_p-j_alpha_p;
	u_alpha_m=tmp_alpha_m-j_alpha_m;
	// What if alpha_p is out of range
	if(j_alpha_p>=(sysparams.AbsPsi_Steps-1)) {
		#if SPLINE_ERR
			cerr<<"alpha_p value is out of range in function evaluateSplinePhi()."<<endl;
		#endif
		u_alpha_p=1.0;
		j_alpha_p=sysparams.AbsPsi_Steps-2;
	}
	// What if alpha_m is out of range
	if(j_alpha_m>=(sysparams.AbsPsi_Steps-1)) {
		#if SPLINE_ERR
			cerr<<"alpha_m value is out of range in function evaluateSplinePhi()."<<endl;
		#endif
		u_alpha_m=1.0;
		j_alpha_m=sysparams.AbsPsi_Steps-2;
	}
	// Fit spline to the transformed matrix of densities
	phi_average*=0.0;
	for(j=0;j<3;j++) {
		phi_average[j][0]=(((1.0-u_SqrtMinusDmu)*gpe.phi_avg[j_SqrtMinusDmu][j_alpha_p  ][j_alpha_m  ][j]+u_SqrtMinusDmu*gpe.phi_avg[j_SqrtMinusDmu+1][j_alpha_p  ][j_alpha_m  ][j])*(1.0-u_alpha_p)
			          +((1.0-u_SqrtMinusDmu)*gpe.phi_avg[j_SqrtMinusDmu][j_alpha_p+1][j_alpha_m  ][j]+u_SqrtMinusDmu*gpe.phi_avg[j_SqrtMinusDmu+1][j_alpha_p+1][j_alpha_m  ][j])*u_alpha_p      )*(1.0-u_alpha_m)
			         +(((1.0-u_SqrtMinusDmu)*gpe.phi_avg[j_SqrtMinusDmu][j_alpha_p  ][j_alpha_m+1][j]+u_SqrtMinusDmu*gpe.phi_avg[j_SqrtMinusDmu+1][j_alpha_p  ][j_alpha_m+1][j])*(1.0-u_alpha_p)
			          +((1.0-u_SqrtMinusDmu)*gpe.phi_avg[j_SqrtMinusDmu][j_alpha_p+1][j_alpha_m+1][j]+u_SqrtMinusDmu*gpe.phi_avg[j_SqrtMinusDmu+1][j_alpha_p+1][j_alpha_m+1][j])*u_alpha_p      )*u_alpha_m;
	};
	U.transpose();
	phi_average=U*phi_average;
	return phi_average;
}
#undef SPLINE_ERR

void initializeSpline(GPE_Parameters gpe,System_Parameters sysparams) {
	register double SqrtMinusDmu,AbsPsi,AbsF;
	register double Delta_AbsPsi2,Delta_AbsF2; // derivation steps
	register double alpha_p,alpha_m;
	register dcomplex v_p[3],v_m[3];
	dcomplex ii(0.0,1.0);
	register dcomplex PsiVec[sysparams.nthreads][3];
	register int i1,i2,i3,j1,j2,j;
	register int threadID;
	register cmatrix ntmp(3,3),phi_tmp(3,1);
	
	// Eigenvectors correcponding to the +1 and -1 eigenvalue of Fz.
	v_p[0]=1.0/sqrt(2.0); v_p[1]=-1.0*ii/sqrt(2.0); v_p[2]=0.0;
	v_m[0]=1.0/sqrt(2.0); v_m[1]=     ii/sqrt(2.0); v_m[2]=0.0;
	
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
		#pragma omp for schedule(dynamic) nowait private(i1,i2,i3,j,j1,j2,SqrtMinusDmu,AbsPsi,AbsF,Delta_AbsPsi2,Delta_AbsF2,alpha_p,alpha_m,ntmp,phi_tmp)
		for(i1=0;i1<sysparams.SqrtMinusDmu_Steps;i1++) {
			SqrtMinusDmu=i1*sysparams.D_SqrtMinusDmu;
			for(i2=0;i2<sysparams.AbsPsi_Steps;i2++) {
				AbsPsi=i2*sysparams.D_AbsPsi;
				Delta_AbsPsi2=__max(2.0*AbsPsi*sysparams.D_AbsPsi*0.01,2.0*sysparams.D_AbsPsi*sysparams.D_AbsPsi*0.01);
				for(i3=0;i3<sysparams.AbsF_Steps;i3++) {
					AbsF=i3*sysparams.D_AbsF;
					Delta_AbsF2=__max(2.0*AbsF*sysparams.D_AbsF*0.01,2.0*sysparams.D_AbsF*sysparams.D_AbsF*0.01);
					if(i3>sysparams.AbsF_Steps/2)
						Delta_AbsF2*=-1.0;
					setPsi(AbsPsi,AbsF,PsiVec[threadID]);
					gpe.SeffLDA[i1][i2][i3]=calculateSeff(SqrtMinusDmu,PsiVec[threadID],sysparams);
					setPsi(sqrt(AbsPsi*AbsPsi+Delta_AbsPsi2),AbsF,PsiVec[threadID]);
					gpe.dSeff_dPsi2[i1][i2][i3]=(calculateSeff(SqrtMinusDmu,PsiVec[threadID],sysparams)
					                            -gpe.SeffLDA[i1][i2][i3])/Delta_AbsPsi2;
					setPsi(AbsPsi,sqrt(AbsF*AbsF+Delta_AbsF2),PsiVec[threadID]);
					gpe.dSeff_dF2  [i1][i2][i3]=(calculateSeff(SqrtMinusDmu,PsiVec[threadID],sysparams)
					                            -gpe.SeffLDA[i1][i2][i3])/Delta_AbsF2;
				};
			};
			// Densities
			#if SAVE_DENSITIES
			for(i2=0;i2<sysparams.AbsPsi_Steps;i2++) {
				alpha_p=i2*sysparams.D_AbsPsi;
				for(i3=0;i3<sysparams.AbsPsi_Steps;i3++) {
					alpha_m=i3*sysparams.D_AbsPsi;
					for(j=0;j<3;j++)
						PsiVec[threadID][j]=alpha_p*v_p[j]+alpha_m*v_m[j];
					ntmp=calculateDensitiesAndPhi(0,SqrtMinusDmu,PsiVec[threadID],sysparams);
					for(j1=0;j1<3;j1++) {
						for(j2=0;j2<3;j2++)
							gpe.n[i1][i2][i3][j1][j2]=ntmp[j1][j2];
					};
				};
			};
			#endif
			for(i2=0;i2<sysparams.AbsPsi_Steps;i2++) {
				alpha_p=i2*sysparams.D_AbsPsi;
				for(i3=0;i3<sysparams.AbsPsi_Steps;i3++) {
					alpha_m=i3*sysparams.D_AbsPsi;
					for(j=0;j<3;j++)
						PsiVec[threadID][j]=alpha_p*v_p[j]+alpha_m*v_m[j];
					phi_tmp=calculateDensitiesAndPhi(1,SqrtMinusDmu,PsiVec[threadID],sysparams);
					for(j1=0;j1<3;j1++)
							gpe.phi_avg[i1][i2][i3][j1]=phi_tmp[j1][0];
				};
			};
		};
	}
}

void calcAbsPsiAndAbsF(dcomplex PsiVec[3],double& AbsPsi,double& AbsF,double FVec[3],System_Parameters sysparams) {
	register int i1,i2;
	register double tmp1,tmp2;
	AbsPsi=0.0; AbsF=0.0; FVec[0]=0.0; FVec[1]=0.0; FVec[2]=0.0;
	for(i1=0;i1<3;i1++) {
		AbsPsi+=PsiVec[i1].re*PsiVec[i1].re+PsiVec[i1].im*PsiVec[i1].im;
		for(i2=0;i2<3;i2++) {
			tmp1=PsiVec[i1].re*PsiVec[i2].re+PsiVec[i1].im*PsiVec[i2].im;
			tmp2=PsiVec[i1].re*PsiVec[i2].im-PsiVec[i1].im*PsiVec[i2].re;
			FVec[0]+=sysparams.Fx[i1][i2].re*tmp1-sysparams.Fx[i1][i2].im*tmp2;
			FVec[1]+=sysparams.Fy[i1][i2].re*tmp1-sysparams.Fy[i1][i2].im*tmp2;
			FVec[2]+=sysparams.Fz[i1][i2].re*tmp1-sysparams.Fz[i1][i2].im*tmp2;
		};
	};
	tmp1=1.0/AbsPsi;
	for(i1=0;i1<3;i1++)
		FVec[i1]*=tmp1;
	AbsPsi=sqrt(AbsPsi);
	AbsF=sqrt(FVec[0]*FVec[0]+FVec[1]*FVec[1]+FVec[2]*FVec[2]);
}

void setPsi(double AbsPsi,double AbsF,dcomplex PsiVec[3]) {
	// Sets a Psi configuration with PsiVec'*PsiVec=AbsPsi^2 
	// and ((PsiVec'*Fx*PsiVec)^2+(PsiVec'*Fy*PsiVec)^2+(PsiVec'*Fz*PsiVec)^2)/(PsiVec'*PsiVec)^2=AbsF^2
	double alpha;
	double tmp;
	alpha=((1.0+AbsF)+sqrt(1.0-AbsF*AbsF))/(2.0*(1.0+AbsF));
	tmp=1.0/sqrt(alpha*alpha+(1.0-alpha)*(1.0-alpha));
	PsiVec[0]=AbsPsi*(1.0-alpha)*tmp;
	PsiVec[1].re=0.0;
	PsiVec[1].im=AbsPsi*alpha*tmp;
	PsiVec[2]=0.0;
}

void transformPsi(dcomplex PsiVec[3],double& alpha_p,double& alpha_m,dmatrix& U) {
	// Calculates the appropriate SO(3) transformation that brings PsiVec into 
	// the subspace spanned by v_p and v_m, the +1 and -1 eigenvectors of Fz.
	// U*PsiVec=(alpha_p*v_p+alpha_m*v_m)*exp(i*tau), where alpha_p and alpha_m
	// are positive numbers.
	int j;
	dcomplex ii(0.0,1.0);
	dcomplex alpha_p_cmp,alpha_m_cmp;
	dcomplex v_p_conj[3],v_m_conj[3];
	double tmp1,tmp2,tmp3;
	double x,y,z,rxy,ryz,rxyz,rxy_,ryz_,rxyz_;
	register dmatrix U1(3,3),U2(3,3),U3(3,3);
	register cmatrix Psi(3,1);
	Psi[0][0]=PsiVec[0]; Psi[1][0]=PsiVec[1]; Psi[2][0]=PsiVec[2];
	tmp1=1.0/sqrt(2.0);
	v_p_conj[0]=tmp1; v_p_conj[1]= tmp1*ii; v_p_conj[2]=0.0;
	v_m_conj[0]=tmp1; v_m_conj[1]=-tmp1*ii; v_m_conj[2]=0.0;
	
	// U1
	x=Psi[0][0].re;
	y=Psi[1][0].re;
	z=Psi[2][0].re;
	rxy=sqrt(x*x+y*y);
	rxyz=sqrt(x*x+y*y+z*z);
	if(rxy && rxyz) {
		rxy_=1.0/rxy;
		rxyz_=1.0/rxyz;
		U1[0][0]=x*rxyz_;            U1[0][1]=y*rxyz_;            U1[0][2]=z*rxyz_;
		U1[1][0]=-y*rxy_;            U1[1][1]=x*rxy_;             U1[1][2]=0.0;
		U1[2][0]=-x*z*rxy_*rxyz_;    U1[2][1]=-y*z*rxy_*rxyz_;    U1[2][2]=rxy*rxyz_;
	}
	else if(!rxy) {
		U1[0][0]=0.0;    U1[0][1]=0.0;     U1[0][2]=1.0;
		U1[1][0]=0.0;    U1[1][1]=1.0;     U1[1][2]=0.0;
		U1[2][0]=1.0;    U1[2][1]=0.0;     U1[2][2]=0.0;
	}
	else if(!rxyz) {
		U1[0][0]=1.0;    U1[0][1]=0.0;     U1[0][2]=0.0;
		U1[1][0]=0.0;    U1[1][1]=1.0;     U1[1][2]=0.0;
		U1[2][0]=0.0;    U1[2][1]=0.0;     U1[2][2]=1.0;
	}
	Psi=U1*Psi;
	
	// U2
	y=Psi[1][0].im;
	z=Psi[2][0].im;
	ryz=sqrt(y*y+z*z);
	if(ryz) {
		ryz_=1.0/ryz;
		tmp2=y*ryz_;
		tmp3=z*ryz_;
		U2[0][0]=1.0;   U2[0][1]=0.0;    U2[0][2]=0.0;
		U2[1][0]=0.0;   U2[1][1]=tmp2;   U2[1][2]=tmp3;
		U2[2][0]=0.0;   U2[2][1]=-tmp3;  U2[2][2]=tmp2;
	}
	else {
		U2[0][0]=1.0;   U2[0][1]=0.0;    U2[0][2]=0.0;
		U2[1][0]=0.0;   U2[1][1]=1.0;    U2[1][2]=0.0;
		U2[2][0]=0.0;   U2[2][1]=0.0;    U2[2][2]=1.0;
	}
	Psi=U2*Psi;
	
	// U3
	for(alpha_p_cmp=0.0,alpha_m_cmp=0.0,j=0;j<3;j++) {
		alpha_p_cmp+=(v_p_conj[j]*Psi[j][0]);
		alpha_m_cmp+=(v_m_conj[j]*Psi[j][0]);
	};
	alpha_p=alpha_p_cmp.abs();
	alpha_m=alpha_m_cmp.abs();
	if(alpha_m) {
		// Auxiliary quantities
		dcomplex tmp_cmp=alpha_p_cmp/alpha_m_cmp;
		double Cos2Khi=tmp_cmp.re/tmp_cmp.abs();
		double sgn=(tmp_cmp.im>=0.0 ? 1.0 : -1.0);
		tmp1=sqrt(0.5*(1.0+Cos2Khi))*sgn;
		tmp2=sqrt(0.5*(1.0-Cos2Khi));
		U3[0][0]=tmp1;   U3[0][1]=tmp2;    U3[0][2]=0.0;
		U3[1][0]=-tmp2;  U3[1][1]=tmp1;    U3[1][2]=0.0;
		U3[2][0]=0.0;    U3[2][1]=0.0;     U3[2][2]=1.0;
	}
	else {
		U3[0][0]=1.0;   U3[0][1]=0.0;    U3[0][2]=0.0;
		U3[1][0]=0.0;   U3[1][1]=1.0;    U3[1][2]=0.0;
		U3[2][0]=0.0;   U3[2][1]=0.0;    U3[2][2]=1.0;
	}
	
	// U
	U3*=U2;
	U3*=U1;
	U=U3;
}

double calculateSeff(double SqrtMinusDmu,dcomplex PsiVec[3],System_Parameters sysparams) {
	// The inclusion of the 3J terms in the redefinition of the chemical potentials is a matter of taste here.
	// Either you choose (-mu+3J)(N1+N2+N3), or just -mu(N1+N2+N3).
	
	double result;
	double mu=sysparams.mu_0-SqrtMinusDmu*SqrtMinusDmu;
	
	// Set up the Hamiltonian
	register int i,i1,i2,j;
	register double beta=1.0/sysparams.T;
	register double coeff;
	register int m[3],n[3];
	register dcomplex **Hamiltonian;
	double* eigenvalues;
	double Ztmp;
	
	Hamiltonian=new dcomplex* [sysparams.sizeOfHamiltonian];
	for(j=0;j<sysparams.sizeOfHamiltonian;j++)
		Hamiltonian[j]=new (nothrow) dcomplex[sysparams.sizeOfHamiltonian];
	eigenvalues=new (nothrow) double[sysparams.sizeOfHamiltonian];
	
	for(i1=0;i1<sysparams.sizeOfHamiltonian;i1++) {
		for(i2=0;i2<sysparams.sizeOfHamiltonian;i2++) {
			Hamiltonian[i1][i2]=sysparams.HInt[i1][i2];
		};
		Hamiltonian[i1][i1].re+=(-mu+3.0*sysparams.J)*(sysparams.N[0][0][i1][i1].re+sysparams.N[1][1][i1][i1].re+sysparams.N[2][2][i1][i1].re);
	};
	for(i=0;i<3;i++) {
		for(i2=0;i2<sysparams.sizeOfHamiltonian;i2++) {
			inv_index(n[0],n[1],n[2],i2,sysparams);
			// Psi_i Phi_i^+
			coeff=sysparams.V0;
			m[0]=n[0]; m[1]=n[1]; m[2]=n[2];
			create(i,m,coeff,sysparams);
			if(coeff!=0.0 && index(m[0],m[1],m[2],i1,sysparams)) {
				Hamiltonian[i1][i2].re+=coeff*PsiVec[i].re;
				Hamiltonian[i1][i2].im+=coeff*PsiVec[i].im;
			}
			// Psi_i^* Phi_i
			coeff=sysparams.V0;
			m[0]=n[0]; m[1]=n[1]; m[2]=n[2];
			annihilate(i,m,coeff,sysparams);
			if(coeff!=0 && index(m[0],m[1],m[2],i1,sysparams)) {
				Hamiltonian[i1][i2].re+=coeff*PsiVec[i].re;
				Hamiltonian[i1][i2].im-=coeff*PsiVec[i].im;
			}
		};
	};
	
	// Diagonalization and calculation of Z and densities
	diagonalize(Hamiltonian,eigenvalues,sysparams.sizeOfHamiltonian);
	for(Ztmp=0.0,j=0;j<sysparams.sizeOfHamiltonian;j++) {
		eigenvalues[j]=exp(-beta*eigenvalues[j]);
		Ztmp+=eigenvalues[j];
	};
	result=-log(Ztmp)+((PsiVec[0].re*PsiVec[0].re+PsiVec[0].im*PsiVec[0].im)
		                  +(PsiVec[1].re*PsiVec[1].re+PsiVec[1].im*PsiVec[1].im)
		                  +(PsiVec[2].re*PsiVec[2].re+PsiVec[2].im*PsiVec[2].im)
		                  )*beta*sysparams.V0*sysparams.V0/(6.0*sysparams.J);
	// Free Hamiltonian
	for(j=0;j<sysparams.sizeOfHamiltonian;j++)
		delete[] Hamiltonian[j];
	delete[] Hamiltonian;
	delete[] eigenvalues;
	
	return result;
}

cmatrix calculateDensitiesAndPhi(int Option,double SqrtMinusDmu,dcomplex PsiVec[3],System_Parameters sysparams) {
	// Option==0: densities
	// Option==1: phi_average
	// The inclusion of the 3J terms in the redefinition of the chemical potentials is a matter of taste here.
	// Either you choose (-mu+3J)(N1+N2+N3), or just -mu(N1+N2+N3).
	
	cmatrix densities(3,3);
	cmatrix phi_average(3,1);
	
	double mu=sysparams.mu_0-SqrtMinusDmu*SqrtMinusDmu;
	
	// Set up the Hamiltonian
	register int i,i1,i2,i3,j,j1,j2;
	register double beta=1.0/sysparams.T;
	register double coeff;
	register int m[3],n[3];
	register dcomplex **Hamiltonian;
	register dcomplex **expH;
	double* eigenvalues;
	double Ztmp;
	dcomplex ctmp;
	
	Hamiltonian=new dcomplex* [sysparams.sizeOfHamiltonian];
	for(j=0;j<sysparams.sizeOfHamiltonian;j++)
		Hamiltonian[j]=new (nothrow) dcomplex[sysparams.sizeOfHamiltonian];
	expH=new dcomplex* [sysparams.sizeOfHamiltonian];
	for(j=0;j<sysparams.sizeOfHamiltonian;j++)
		expH[j]=new (nothrow) dcomplex[sysparams.sizeOfHamiltonian];
	eigenvalues=new (nothrow) double[sysparams.sizeOfHamiltonian];
	
	for(i1=0;i1<sysparams.sizeOfHamiltonian;i1++) {
		for(i2=0;i2<sysparams.sizeOfHamiltonian;i2++) {
			Hamiltonian[i1][i2]=sysparams.HInt[i1][i2];
		};
		Hamiltonian[i1][i1].re+=(-mu+3.0*sysparams.J)*(sysparams.N[0][0][i1][i1].re+sysparams.N[1][1][i1][i1].re+sysparams.N[2][2][i1][i1].re);
	};
	for(i=0;i<3;i++) {
		for(i2=0;i2<sysparams.sizeOfHamiltonian;i2++) {
			inv_index(n[0],n[1],n[2],i2,sysparams);
			// Psi_i Phi_i^+
			coeff=sysparams.V0;
			m[0]=n[0]; m[1]=n[1]; m[2]=n[2];
			create(i,m,coeff,sysparams);
			if(coeff!=0.0 && index(m[0],m[1],m[2],i1,sysparams)) {
				Hamiltonian[i1][i2].re+=coeff*PsiVec[i].re;
				Hamiltonian[i1][i2].im+=coeff*PsiVec[i].im;
			}
			// Psi_i^* Phi_i
			coeff=sysparams.V0;
			m[0]=n[0]; m[1]=n[1]; m[2]=n[2];
			annihilate(i,m,coeff,sysparams);
			if(coeff!=0 && index(m[0],m[1],m[2],i1,sysparams)) {
				Hamiltonian[i1][i2].re+=coeff*PsiVec[i].re;
				Hamiltonian[i1][i2].im-=coeff*PsiVec[i].im;
			}
		};
	};
	
	// Diagonalization and calculation of Z and densities
	diagonalize(Hamiltonian,eigenvalues,sysparams.sizeOfHamiltonian);
	for(Ztmp=0.0,j=0;j<sysparams.sizeOfHamiltonian;j++) {
		eigenvalues[j]=exp(-beta*eigenvalues[j]);
		Ztmp+=eigenvalues[j];
	};
	
	// Calculate the exponentiated beta*Hamiltonian, 'expH'
	for(i1=0;i1<sysparams.sizeOfHamiltonian;i1++) {
		for(i2=0;i2<sysparams.sizeOfHamiltonian;i2++) {
			expH[i1][i2]=0.0;
			for(i3=0;i3<sysparams.sizeOfHamiltonian;i3++) {
				ctmp.re=Hamiltonian[i2][i3].re;
				ctmp.im=-Hamiltonian[i2][i3].im;
				expH[i1][i2]+=Hamiltonian[i1][i3]*eigenvalues[i3]*ctmp;
			};
		};
	};
	if(!Option) {
		densities*=0.0;
		for(j1=0;j1<3;j1++) {
			for(j2=0;j2<3;j2++) {
				for(i1=0;i1<sysparams.sizeOfHamiltonian;i1++) {
					for(i2=0;i2<sysparams.sizeOfHamiltonian;i2++)
						densities[j1][j2]+=sysparams.N[j1][j2][i1][i2]*expH[i2][i1];
				};
			};
		};
		densities*=1.0/Ztmp;
	}
	else {
		phi_average*=0.0;
		for(j1=0;j1<3;j1++) {
			for(i1=0;i1<sysparams.sizeOfHamiltonian;i1++) {
				for(i2=0;i2<sysparams.sizeOfHamiltonian;i2++)
					phi_average[j1][0]+=sysparams.phi[j1][i1][i2]*expH[i2][i1];
			};
		};
		phi_average*=1.0/Ztmp;
	}

	// Free Hamiltonian
	for(j=0;j<sysparams.sizeOfHamiltonian;j++) {
		delete[] Hamiltonian[j];
		delete[] expH[j];
	}
	delete[] Hamiltonian;
	delete[] expH;
	delete[] eigenvalues;
	
	if(!Option) {
		return densities;
	}
	else {
		return phi_average;
	}
}

void phiAvereages(fftw_complex* phiAverage1,fftw_complex* phiAverage2,fftw_complex* phiAverage3,
		  GPE_Psis Psis,GPE_Parameters gpe,System_Parameters sysparams) {
	// Calculates the expectation values of the annihilation operators phi1, phi2 and phi3,
	// and puts the results in phiAverage1, phiAverage2 etc.
	// These arrays shall be of the size sysparams.M.
	register int j,jx,jy,jz;
	cmatrix phi_avg(3,1);;
	dcomplex PsiVec[3];
	register double diffjx,diffjy,diffjz;
	register double SqrtMinusDmu;
	
	for(j=0,jx=0;jx<sysparams.Mx;jx++) {
		diffjx=jx-0.5*(sysparams.Mx-1.0);
		for(jy=0;jy<sysparams.My;jy++) {
			diffjy=jy-0.5*(sysparams.My-1.0);
			for(jz=0;jz<sysparams.Mz;jz++,j++) {
				diffjz=jz-0.5*(sysparams.Mz-1.0);
				SqrtMinusDmu=sqrt(sysparams.trapx*diffjx*diffjx
			                         +sysparams.trapy*diffjy*diffjy
				                 +sysparams.trapz*diffjz*diffjz);
				PsiVec[0]=Psis.Psi1[j];
				PsiVec[1]=Psis.Psi2[j];
				PsiVec[2]=Psis.Psi3[j];
				phi_avg=evaluateSplinePhi(SqrtMinusDmu,PsiVec,gpe,sysparams);
				phiAverage1[j][0]=phi_avg[0][0].re; phiAverage1[j][1]=phi_avg[0][0].im;
				phiAverage2[j][0]=phi_avg[1][0].re; phiAverage2[j][1]=phi_avg[1][0].im;
				phiAverage3[j][0]=phi_avg[2][0].re; phiAverage3[j][1]=phi_avg[2][0].im;
			};
		};
	};
}


void densityCorrelations(fftw_complex* n11,fftw_complex* n12,fftw_complex* n13,fftw_complex* n22,fftw_complex* n23,fftw_complex* n33,
                         GPE_Psis Psis,GPE_Parameters gpe,System_Parameters sysparams,FFTW_Struct FFTW_Variables) {
	// The arrays n11, n12 etc. shall be of the size sysparams.M.
	register int j,jx,jy,jz;
	cmatrix phi_avg(3,1),n(3,3);
	dcomplex phi_avg_conj[3];
	dcomplex PsiVec[3];
	register double diffjx,diffjy,diffjz;
	register double SqrtMinusDmu;
	dcomplex constMx[3][3];
	register double tmp;
	fftw_complex* tmpArray;
	tmpArray=new (nothrow) fftw_complex [sysparams.M];
	// Calculate the expectation values of the annihilation operators
	for(jx=0;jx<3;jx++) {
		for(jy=0;jy<3;jy++) {
			constMx[jx][jy]=0.0;
		};
	};
	for(j=0,jx=0;jx<sysparams.Mx;jx++) {
		diffjx=jx-0.5*(sysparams.Mx-1.0);
		for(jy=0;jy<sysparams.My;jy++) {
			diffjy=jy-0.5*(sysparams.My-1.0);
			for(jz=0;jz<sysparams.Mz;jz++,j++) {
				diffjz=jz-0.5*(sysparams.Mz-1.0);
				SqrtMinusDmu=sqrt(sysparams.trapx*diffjx*diffjx
			                         +sysparams.trapy*diffjy*diffjy
				                 +sysparams.trapz*diffjz*diffjz);
				PsiVec[0]=Psis.Psi1[j];
				PsiVec[1]=Psis.Psi2[j];
				PsiVec[2]=Psis.Psi3[j];
				phi_avg=evaluateSplinePhi(SqrtMinusDmu,PsiVec,gpe,sysparams);
				n11[j][0]=phi_avg[0][0].re; n11[j][1]=phi_avg[0][0].im;
				n22[j][0]=phi_avg[1][0].re; n22[j][1]=phi_avg[1][0].im;
				n33[j][0]=phi_avg[2][0].re; n33[j][1]=phi_avg[2][0].im;
				phi_avg_conj[0].re=phi_avg[0][0].re; phi_avg_conj[0].im=-phi_avg[0][0].im; 
				phi_avg_conj[1].re=phi_avg[1][0].re; phi_avg_conj[1].im=-phi_avg[1][0].im; 
				phi_avg_conj[2].re=phi_avg[2][0].re; phi_avg_conj[2].im=-phi_avg[2][0].im; 
				
				// Calculate the constant part of the n(k): constMx
				n=evaluateSplineDensities(SqrtMinusDmu,PsiVec,gpe,sysparams);
				constMx[0][0]+=n[0][0]-phi_avg_conj[0]*phi_avg[0][0];
				constMx[0][1]+=n[0][1]-phi_avg_conj[0]*phi_avg[1][0];
				constMx[0][2]+=n[0][2]-phi_avg_conj[0]*phi_avg[2][0];
				constMx[1][1]+=n[1][1]-phi_avg_conj[1]*phi_avg[1][0];
				constMx[1][2]+=n[1][2]-phi_avg_conj[1]*phi_avg[2][0];
				constMx[2][2]+=n[2][2]-phi_avg_conj[2]*phi_avg[2][0];
			};
		};
	};
	
	//Calculate the 'Fourier transform' x 'sqrt(Mx*My*Mz)' of the expectation values of the annihilation operators 
	//stored currently in n11 etc. and overwrite n11, n22 and n33 with them.
	fftw_execute_dft(FFTW_Variables.plan,n11,n11);
	fftw_execute_dft(FFTW_Variables.plan,n22,n22);
	fftw_execute_dft(FFTW_Variables.plan,n33,n33);
	tmp=1.0/sysparams.M;
	for(j=0;j<sysparams.M;j++) {
		n12[j][0]=(n11[j][0]*n22[j][0]+n11[j][1]*n22[j][1]+constMx[0][1].re)*tmp;
		n12[j][1]=(n11[j][0]*n22[j][1]-n11[j][1]*n22[j][0]+constMx[0][1].im)*tmp;
		n13[j][0]=(n11[j][0]*n33[j][0]+n11[j][1]*n33[j][1]+constMx[0][2].re)*tmp;
		n13[j][1]=(n11[j][0]*n33[j][1]-n11[j][1]*n33[j][0]+constMx[0][2].im)*tmp;
		n23[j][0]=(n22[j][0]*n33[j][0]+n22[j][1]*n33[j][1]+constMx[1][2].re)*tmp;
		n23[j][1]=(n22[j][0]*n33[j][1]-n22[j][1]*n33[j][0]+constMx[1][2].im)*tmp;
		
		n11[j][0]=(n11[j][0]*n11[j][0]+n11[j][1]*n11[j][1]+constMx[0][0].re)*tmp;
		n11[j][1]=0.0;
		n22[j][0]=(n22[j][0]*n22[j][0]+n22[j][1]*n22[j][1]+constMx[1][1].re)*tmp;
		n22[j][1]=0.0;
		n33[j][0]=(n33[j][0]*n33[j][0]+n33[j][1]*n33[j][1]+constMx[2][2].re)*tmp;
		n33[j][1]=0.0;
	};
	
	// Shift the functions in Fourier space
	register int jxOld,jyOld,jzOld,jOld;
	register int jxNew,jyNew,jzNew,jNew;
	register int tmp_int_x=sysparams.Mx/2,tmp_int_y=sysparams.My/2,tmp_int_z=sysparams.Mz/2;
	register fftw_complex* Array;
	for (j=0;j<6;j++) {
		switch (j) {
		  case 0:
			Array=n11;
			break;
		  case 1:
			Array=n12;
			break;
		  case 2:
			Array=n13;
			break;
		  case 3:
			Array=n22;
			break;
		  case 4:
			Array=n23;
			break;
		  case 5:
			Array=n33;
			break;
		}
		for(jOld=0;jOld<sysparams.M;jOld++) {
			tmpArray[jOld][0]=Array[jOld][0];
			tmpArray[jOld][1]=Array[jOld][1];
		};
		for(jOld=0,jxOld=0;jxOld<sysparams.Mx;jxOld++) {
			jxNew=(jxOld<tmp_int_x ? jxOld+tmp_int_x : jxOld-tmp_int_x);
			for(jyOld=0;jyOld<sysparams.My;jyOld++) {
				jyNew=(jyOld<tmp_int_y ? jyOld+tmp_int_y : jyOld-tmp_int_y);
				for(jzOld=0;jzOld<sysparams.Mz;jzOld++,jOld++) {
					jzNew=(jzOld<tmp_int_z ? jzOld+tmp_int_z : jzOld-tmp_int_z);
					jNew=jzNew+sysparams.Mz*(jyNew+sysparams.My*jxNew);
					Array[jNew][0]=tmpArray[jOld][0];
					Array[jNew][1]=tmpArray[jOld][1];
				};
			};
		};
	};
	delete [] tmpArray;
}

void resetSplineParameters(System_Parameters& sysparams) {
	// Sets the boundaries of the spline to appropriate values.
	// Finds the minimum of Seff in terms of AbsPsi and AbsF, and 
	// sets sysparams.AbsPsi_Max to be twice the value of this.
	// It also sets sysparams.SqrtMinusDmu_Max to be the maximal 
	// possible value at the edge of the cloud.
	
	// Number of grid points to search for the minimum.
	int AbsPsi_SearchSteps=31;
	int AbsF_SearchSteps=26;
	int SqrtMinusDmu_SearchSteps=51;
	
	// Variables
	double SqrtMinusDmu,AbsPsi,AbsF;
	double D_SqrtMinusDmu_Search,D_AbsPsi_Search,D_AbsF_Search; // stepsizes in the search
	dcomplex PsiVec[sysparams.nthreads][3];
	double greatestAbsPsiMinPoint[sysparams.nthreads];
	double tmp1,tmp2;
	
	// Counters
	int i1,i2,i3,jMin;
	int threadID;
	
	// Reset SqrtMinusDmu_Max to the optimal value
	sysparams.SqrtMinusDmu_Max=
	         (sysparams.SqrtMinusDmu_Steps+1.0)/(sysparams.SqrtMinusDmu_Steps-1.0)
	         *sqrt(sysparams.trapx*0.25*(sysparams.Mx-1.0)*(sysparams.Mx-1.0)+
	               sysparams.trapy*0.25*(sysparams.My-1.0)*(sysparams.My-1.0)+
	               sysparams.trapz*0.25*(sysparams.Mz-1.0)*(sysparams.Mz-1.0));
	sysparams.D_SqrtMinusDmu=sysparams.SqrtMinusDmu_Max/(sysparams.SqrtMinusDmu_Steps-1.0);
	
	// Reset sysparams.D_AbsF
	sysparams.D_AbsF=1.0/(sysparams.AbsF_Steps-1.0);
	
	// Set sysparams.AbsPsi_Max to be twice the value of the greatest AbsPsi at the LDA minimum.
	D_SqrtMinusDmu_Search=sysparams.SqrtMinusDmu_Max/(SqrtMinusDmu_SearchSteps-1.0);
	D_AbsPsi_Search      =sysparams.AbsPsi_Max/(AbsPsi_SearchSteps-1.0);
	D_AbsF_Search        =1.0/(AbsF_SearchSteps-1.0);
	for(i1=0;i1<sysparams.nthreads;i1++)
		greatestAbsPsiMinPoint[i1]=0.0;
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
		#pragma omp for schedule(dynamic) nowait private(i1,i2,i3,jMin,SqrtMinusDmu,AbsPsi,AbsF,tmp1,tmp2)
		for(i1=0;i1<SqrtMinusDmu_SearchSteps;i1++){
			SqrtMinusDmu=i1*D_SqrtMinusDmu_Search;
			for(i3=0;i3<AbsF_SearchSteps;i3++) {
				AbsF=i3*D_AbsF_Search;
				PsiVec[threadID][0]=0.0;
				PsiVec[threadID][1]=0.0;
				PsiVec[threadID][2]=0.0;
				tmp1=calculateSeff(SqrtMinusDmu,PsiVec[threadID],sysparams);
				for(jMin=0,i2=1;i2<AbsPsi_SearchSteps;i2++) {
					AbsPsi=i2*D_AbsPsi_Search;
					setPsi(AbsPsi,AbsF,PsiVec[threadID]);
					tmp2=calculateSeff(SqrtMinusDmu,PsiVec[threadID],sysparams);
					if(isnan(tmp2) || isinf(tmp2)) {tmp2=1e10*(1.0+AbsPsi);}
					if(tmp2<tmp1) { jMin=i2; tmp1=tmp2;}
				};
				if(jMin==(AbsPsi_SearchSteps-1)) {
					do {
						tmp1=tmp2;
						AbsPsi+=D_AbsPsi_Search;
						setPsi(AbsPsi,AbsF,PsiVec[threadID]);
						tmp2=calculateSeff(SqrtMinusDmu,PsiVec[threadID],sysparams);
						if(isnan(tmp2) || isinf(tmp2)) {tmp2=1e10*(1.0+AbsPsi);}
					}while (tmp2<tmp1);
				}
				AbsPsi=jMin*D_AbsPsi_Search;
				if(AbsPsi>greatestAbsPsiMinPoint[threadID]) {greatestAbsPsiMinPoint[threadID]=AbsPsi;}
			};
		};
	}
	for(i1=1;i1<sysparams.nthreads;i1++) {
		greatestAbsPsiMinPoint[0]=__max(greatestAbsPsiMinPoint[0],greatestAbsPsiMinPoint[i1]);
	};
	if(greatestAbsPsiMinPoint[0]>0.0) { sysparams.AbsPsi_Max=2.0*greatestAbsPsiMinPoint[0]; }
	else { sysparams.AbsPsi_Max=0.05; }
	sysparams.D_AbsPsi=sysparams.AbsPsi_Max/(sysparams.AbsPsi_Steps-1.0);
}

void createHIntAndNs(System_Parameters& sysparams) {
	register int i,i1,i2,i3,i4,j,j1,j2;
	register int m[3],n[3];
	register double coeff;
	register double tmp; 
	register dcomplex tmp_dcomplex;
	register dcomplex **tmpMx1,**tmpMx2,**tmpMx3;
	
	// Calculate the size of these matrices
	for(sysparams.sizeOfHamiltonian=0,i=sysparams.Nmax+2;i>=2;i--) {
		sysparams.sizeOfHamiltonian+=(i*(i-1))/2;
	};
	
	// Allocate HInt, N and phi
	sysparams.N=new (nothrow) dcomplex*** [3];
	for(j1=0;j1<3;j1++) {
		sysparams.N[j1]=new (nothrow) dcomplex** [3];
		for(j2=0;j2<3;j2++) {
			sysparams.N[j1][j2]=new (nothrow) dcomplex* [sysparams.sizeOfHamiltonian];
			for(i1=0;i1<sysparams.sizeOfHamiltonian;i1++)
				sysparams.N[j1][j2][i1]=new (nothrow) dcomplex [sysparams.sizeOfHamiltonian];
		};
	};
	sysparams.phi=new (nothrow) dcomplex** [3];
	for(j1=0;j1<3;j1++) {
		sysparams.phi[j1]=new (nothrow) dcomplex* [sysparams.sizeOfHamiltonian];
		for(j2=0;j2<sysparams.sizeOfHamiltonian;j2++)
			sysparams.phi[j1][j2]=new (nothrow) dcomplex [sysparams.sizeOfHamiltonian];
	};
	sysparams.HInt=new (nothrow) dcomplex *[sysparams.sizeOfHamiltonian];
	for(i=0;i<sysparams.sizeOfHamiltonian;i++)
		sysparams.HInt[i]=new (nothrow) dcomplex[sysparams.sizeOfHamiltonian];
	
	tmpMx1=new (nothrow) dcomplex *[sysparams.sizeOfHamiltonian];
	for(i=0;i<sysparams.sizeOfHamiltonian;i++)
		tmpMx1[i]=new (nothrow) dcomplex[sysparams.sizeOfHamiltonian];
	tmpMx2=new (nothrow) dcomplex *[sysparams.sizeOfHamiltonian];
	for(i=0;i<sysparams.sizeOfHamiltonian;i++)
		tmpMx2[i]=new (nothrow) dcomplex[sysparams.sizeOfHamiltonian];
	tmpMx3=new (nothrow) dcomplex *[sysparams.sizeOfHamiltonian];
	for(i=0;i<sysparams.sizeOfHamiltonian;i++)
		tmpMx3[i]=new (nothrow) dcomplex[sysparams.sizeOfHamiltonian];
	
	// Generate N
	for(j1=0;j1<3;j1++) {
		for(j2=0;j2<3;j2++) {
			for(i1=0;i1<sysparams.sizeOfHamiltonian;i1++){
				for(i2=0;i2<sysparams.sizeOfHamiltonian;i2++)
					sysparams.N[j1][j2][i1][i2]=0.0;
			};
			for(i2=0;i2<sysparams.sizeOfHamiltonian;i2++){
				coeff=1.0;
				inv_index(n[0],n[1],n[2],i2,sysparams);
				annihilate(j2,n,coeff,sysparams);
				create(j1,n,coeff,sysparams);
				if(coeff && index(n[0],n[1],n[2],i1,sysparams)) {
					sysparams.N[j1][j2][i1][i2].re=coeff;
				}
			};
		};
	};
	
	// Generate phi
	for(j1=0;j1<3;j1++) {
		for(i1=0;i1<sysparams.sizeOfHamiltonian;i1++){
			for(i2=0;i2<sysparams.sizeOfHamiltonian;i2++)
				sysparams.phi[j1][i1][i2]=0.0;
		};
		for(i2=0;i2<sysparams.sizeOfHamiltonian;i2++){
			coeff=1.0;
			inv_index(n[0],n[1],n[2],i2,sysparams);
			annihilate(j1,n,coeff,sysparams);
			if(coeff && index(n[0],n[1],n[2],i1,sysparams)) {
				sysparams.phi[j1][i1][i2].re=coeff;
			}
		};
	};
	
	
	// Generate interaction part of the Hamiltonian
	for(i=0;i<sysparams.sizeOfHamiltonian;i++) {
		for(j=0;j<sysparams.sizeOfHamiltonian;j++) {
			sysparams.HInt[i][j]=0.0;
		};
	};
	for(i1=0;i1<3;i1++) {
		for(i2=0;i2<3;i2++) {
			// coefficient of (phi^+_i1 phi^+_i2 phi_i2 phi_i1).
			tmp=0.5*sysparams.V0;
			MxAdjoint(sysparams.phi[i1],tmpMx1,sysparams.sizeOfHamiltonian);
			MxAdjoint(sysparams.phi[i2],tmpMx2,sysparams.sizeOfHamiltonian);
			multMatrices(tmpMx1,tmpMx2,tmpMx3,sysparams.sizeOfHamiltonian);
			multMatrices(tmpMx3,sysparams.phi[i2],tmpMx1,sysparams.sizeOfHamiltonian);
			multMatrices(tmpMx1,sysparams.phi[i1],tmpMx2,sysparams.sizeOfHamiltonian);
			for(j1=0;j1<sysparams.sizeOfHamiltonian;j1++) {
				for(j2=0;j2<sysparams.sizeOfHamiltonian;j2++) {
					sysparams.HInt[j1][j2]+=tmp*tmpMx2[j1][j2];
				};
			};
		};
	};
	for(i1=0;i1<3;i1++) {
		for(i2=0;i2<3;i2++) {
			for(i3=0;i3<3;i3++) {
				for(i4=0;i4<3;i4++) {
					// coefficient of (phi^+_i1 phi^+_i2 phi_i3 phi_i4).
					tmp_dcomplex=0.5*sysparams.V2*(sysparams.Fx[i1][i4]*sysparams.Fx[i2][i3]+
					                               sysparams.Fy[i1][i4]*sysparams.Fy[i2][i3]+
					                               sysparams.Fz[i1][i4]*sysparams.Fz[i2][i3]);
					MxAdjoint(sysparams.phi[i1],tmpMx1,sysparams.sizeOfHamiltonian);
					MxAdjoint(sysparams.phi[i2],tmpMx2,sysparams.sizeOfHamiltonian);
					multMatrices(tmpMx1,tmpMx2,tmpMx3,sysparams.sizeOfHamiltonian);
					multMatrices(tmpMx3,sysparams.phi[i3],tmpMx1,sysparams.sizeOfHamiltonian);
					multMatrices(tmpMx1,sysparams.phi[i4],tmpMx2,sysparams.sizeOfHamiltonian);
					for(j1=0;j1<sysparams.sizeOfHamiltonian;j1++) {
						for(j2=0;j2<sysparams.sizeOfHamiltonian;j2++) {
							sysparams.HInt[j1][j2]+=tmp_dcomplex*tmpMx2[j1][j2];
						};
					};
				};
			};
		};
	};
	
	for(j=0;j<sysparams.sizeOfHamiltonian;j++) {
		delete [] tmpMx1[j];
		delete [] tmpMx2[j];
		delete [] tmpMx3[j];
	};
	delete [] tmpMx1;
	delete [] tmpMx2;
	delete [] tmpMx3;
}

void multMatrices(dcomplex** A,dcomplex** B,dcomplex** C,int sizeMx) {
	// Performs the matrix multiplication A*B=C on the matrices of size 'sizeMx',
	register int j1,j2,j3;
	for(j1=0;j1<sizeMx;j1++) {
		for(j2=0;j2<sizeMx;j2++) {
			for(C[j1][j2]=0,j3=0;j3<sizeMx;j3++) {
				C[j1][j2]+=A[j1][j3]*B[j3][j2];
			};
		};
	};
}

void MxAdjoint(dcomplex** A,dcomplex** A_dagger,int sizeMx) {
	register int j1,j2;
	for(j1=0;j1<sizeMx;j1++) {
		for(j2=0;j2<sizeMx;j2++) {
			A_dagger[j1][j2].re=A[j2][j1].re;
			A_dagger[j1][j2].im=-A[j2][j1].im;
		};
	};
}

void freeHIntAndNs(System_Parameters sysparams) {
	register int i1,i2,i3;
	for(i1=0;i1<3;i1++){
		for(i2=0;i2<3;i2++){
			for(i3=0;i3<sysparams.sizeOfHamiltonian;i3++)
				delete [] sysparams.N[i1][i2][i3];
			delete [] sysparams.N[i1][i2];
		};
		delete[] sysparams.N[i1];
	};
	delete[] sysparams.N;
	
	for(i1=0;i1<3;i1++){
		for(i2=0;i2<sysparams.sizeOfHamiltonian;i2++)
			delete [] sysparams.phi[i1][i2];
		delete[] sysparams.phi[i1];
	};
	delete[] sysparams.phi;
	
	for(i1=0;i1<sysparams.sizeOfHamiltonian;i1++){
		delete[] sysparams.HInt[i1];
	};
	delete[] sysparams.HInt;
}

void annihilate(int flavor,int n[3],double& coefficient,System_Parameters sysparams) {
	if(n[flavor]>0 && n[flavor]<=sysparams.Nmax) {
		coefficient*=sqrt((double) n[flavor]);
		n[flavor]--;
	}
	else {
		coefficient=0.0;
		n[0]=0; n[1]=0; n[2]=0;
	}
}

void create(int flavor,int n[3],double& coefficient,System_Parameters sysparams) {
	if(n[flavor]>=0 && n[flavor]<sysparams.Nmax) {
		coefficient*=sqrt(n[flavor]+1.0);
		n[flavor]++;
	}
	else {
		coefficient=0.0;
		n[0]=0; n[1]=0; n[2]=0;
	}
}

int index(int n1,int n2,int n3,int& n,System_Parameters sysparams) {
	// (n1,n2,n3) --> n
	int i1,i2,i3;
	for(n=0,i1=0;i1<=sysparams.Nmax;i1++) {
		for(i2=0;i2<=sysparams.Nmax-i1;i2++) {
			for(i3=0;i3<=sysparams.Nmax-(i1+i2);i3++) {
				if(i1==n1 && i2==n2 && i3==n3) {
					return 1;
				}
				n++;
			};
		};
	};
	n=-1;
	return 0;
}

int inv_index(int& n1,int& n2,int& n3,int n,System_Parameters sysparams) {
	// n --> (n1,n2,n3)
	int i;
	for(i=n,n1=0;n1<=sysparams.Nmax;n1++) {
		for(n2=0;n2<=sysparams.Nmax-n1;n2++) {
			for(n3=0;n3<=sysparams.Nmax-(n1+n2);n3++) {
				if(i==0) {
					return 1;
				}
				i--;
			};
		};
	};
	n1=-1; n2=-1; n3=-1;
	return 0;
}

// ZHEEV prototype 
extern "C" void zheev_( char* jobz, char* uplo, int* n, dcomplex* a, int* lda,
			double* w, dcomplex* work, int* lwork, double* rwork, int* info );

int diagonalize(dcomplex **HermitianMx,double* eigenvalues,int sizeMx) {
	// Diagonalizes HermitianMx. The eigenvectors will be put into the array HermitianMx.
	// HermitianMx=U*eigenvalues*U'
	// HermitianMx <-- U
	int returnValue=1;
	dcomplex *A;
	A=new (nothrow) dcomplex[sizeMx*sizeMx];
	int n=sizeMx; 
	int lda=sizeMx; 		/* The leading of the array A. lda >= max(1, N).*/
	int info;			/* info == 0:  successful exit
				        < 0:  if INFO == -i, the i-th argument had an illegal value
				        > 0:  if INFO == i, the algorithm failed to converge; i off-diagonal elements of an intermediate
							   tridiagonal form did not converge to zero.*/
	int lwork;			/* Length of the array work. */
	dcomplex wkopt;			/* The variable used to get the size of the optimal workspace array. */
	dcomplex* work;			/* Workspace */
	double rwork[3*sizeMx-2]; 	/* rwork dimension should be at least max(1, 3*n-2). */

	int i1,i2;
	for(i1=0;i1<sizeMx;i1++) {
		for(i2=i1;i2<sizeMx;i2++)
			A[i1*sizeMx+i2]=HermitianMx[i2][i1];
		for(i2=0;i2<i1;i2++)
			A[i1*sizeMx+i2]=0.0;
	};
	// Calculate the optimal workspace to be allocated.
	lwork = -1;
	zheev_("Vectors","Lower",&n,A,&lda,eigenvalues,&wkopt,&lwork,rwork,&info);
	lwork=(int)wkopt.re;
	work=new (nothrow) dcomplex[lwork];
	// Diagonalize
	zheev_("Vectors","Lower",&n,A,&lda,eigenvalues,work,&lwork,rwork,&info);
	if(info>0) {
		cerr<<"zheev_ error.";
		returnValue=0;
	}
	// Put eigenvectors into the array HermitianMx.
	for(i1=0;i1<sizeMx;i1++) {
		for(i2=0;i2<sizeMx;i2++) {
			HermitianMx[i2][i1]=A[i1*sizeMx+i2];
		};
	};
	// Free workspace
	delete[] A;
	delete[] work;
	
	return returnValue;
}

/*
____________________________zheev documentation____________________________

      SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
*  Purpose
*  =======
*
*  ZHEEV computes all eigenvalues and, optionally, eigenvectors of a
*  complex Hermitian matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,2*N-1).
*          For optimal efficiency, LWORK >= (NB+1)*N,
*          where NB is the blocksize for ZHETRD returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(1, 3*N-2))
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
____________________________________________________________________________________
*/