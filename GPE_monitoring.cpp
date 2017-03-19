#include<GPE_dataStructures.h>
#include<GPE_solverFunctions.h>
#include<GPE_monitoring.h>
#include<GPE_spline.h>

void monitoring(const char filenameMonitoring[],GPE_Psis Psis,GPE_Parameters gpe,
		System_Parameters sysparams,FFTW_Struct FFTW_Variables) {
	double kinEnergy1,kinEnergy2,kinEnergy3,potEnergy;
	double Lx,Ly,Lz;
	double normPsi1,normPsi2,normPsi3;
	double tmp;
	FILE* fp;
	
	kineticEnergy(&kinEnergy1,&kinEnergy2,&kinEnergy3,Psis,sysparams,FFTW_Variables);
	potentialEnergy(&potEnergy,Psis,sysparams,gpe);
	angularMomenta(&Lx,&Ly,&Lz,Psis,sysparams);
	norms(&normPsi1,&normPsi2,&normPsi3,Psis,sysparams);
		   
	tmp=1.0/(normPsi1*normPsi1+normPsi2*normPsi2+normPsi3*normPsi3);
	
	// Print to file
	fp=fopen(filenameMonitoring,"a");
	fprintf(fp,"t=%6.5f",gpe.t);
	fprintf(fp,"\t<Energy>=%6.7f",(kinEnergy1+kinEnergy2+kinEnergy3+potEnergy)*tmp);
	fprintf(fp,"\t<Ekin>=%6.5f",(kinEnergy1+kinEnergy2+kinEnergy3)*tmp);
	fprintf(fp,"\t<Epot>=%6.5f",potEnergy*tmp);
	fprintf(fp,"\t<Lx>=%6.5f\t<Ly>=%6.5f\t<Lz>=%6.5f",Lx*tmp,Ly*tmp,Lz*tmp);
	fprintf(fp,"\t|Psi1|_2=%6.5f\t|Psi2|_2=%6.5f\t|Psi3|_2=%6.5f\n",normPsi1,normPsi2,normPsi3);
	fclose(fp);
	
	// Print to the screen
	printf("t=%6.5f",gpe.t);
	printf("\t<Energy>=%6.7f",(kinEnergy1+kinEnergy2+kinEnergy3+potEnergy)*tmp);
	printf("\t<Ekin>=%6.5f",(kinEnergy1+kinEnergy2+kinEnergy3)*tmp);
	printf("\t<Epot>=%6.5f",potEnergy*tmp);
	printf("\t<Lx>=%6.5f\t<Ly>=%6.5f\t<Lz>=%6.5f",Lx*tmp,Ly*tmp,Lz*tmp);
	printf("\t|Psi1|_2=%6.5f\t|Psi2|_2=%6.5f\t|Psi3|_2=%6.5f\n",normPsi1,normPsi2,normPsi3);
}


#define ContinuumLaplacianSpectrum
// #define DiscreteLaplacianSpectrum
void kineticEnergy(double *kinEnergy1,double *kinEnergy2,double *kinEnergy3,GPE_Psis Psis,
	           System_Parameters sysparams,FFTW_Struct FFTW_Variables) {
	register int jloc,jx,jy,jz;
	register double MASS_x,MASS_y,MASS_z;
	double spectx[sysparams.Mx],specty[sysparams.My],spectz[sysparams.Mz];
	register double tmp1,tmp2;
	MASS_x=6.0*6.0*sysparams.J*sysparams.T/(2.0*sysparams.ax*sysparams.ax*sysparams.V0*sysparams.V0);
	MASS_y=6.0*6.0*sysparams.J*sysparams.T/(2.0*sysparams.ay*sysparams.ay*sysparams.V0*sysparams.V0);
	MASS_z=6.0*6.0*sysparams.J*sysparams.T/(2.0*sysparams.az*sysparams.az*sysparams.V0*sysparams.V0);
	
	#ifdef ContinuumLaplacianSpectrum
		// Continuum Laplacian spectrum
		for(jloc=0;jloc<(sysparams.Mx/2);jloc++){
			tmp1=PI*jloc/((sysparams.Mx-1)*sysparams.ax);
			spectx[jloc]=2.0/MASS_x*tmp1*tmp1;
		};
		for(jloc=(sysparams.Mx/2);jloc<sysparams.Mx;jloc++) {
			tmp1=PI*(jloc-sysparams.Mx)/((sysparams.Mx-1)*sysparams.ax);
			spectx[jloc]=2.0/MASS_x*tmp1*tmp1;
		};
		for(jloc=0;jloc<(sysparams.My/2);jloc++){
			tmp1=PI*jloc/((sysparams.My-1)*sysparams.ay);
			specty[jloc]=2.0/MASS_y*tmp1*tmp1;
		};
		for(jloc=(sysparams.My/2);jloc<sysparams.My;jloc++) {
			tmp1=PI*(jloc-sysparams.My)/((sysparams.My-1)*sysparams.ay);
			specty[jloc]=2.0/MASS_y*tmp1*tmp1;
		};
		for(jloc=0;jloc<(sysparams.Mz/2);jloc++){
			tmp1=PI*jloc/((sysparams.Mz-1)*sysparams.az);
			spectz[jloc]=2.0/MASS_z*tmp1*tmp1;
		};
		for(jloc=(sysparams.Mz/2);jloc<sysparams.Mz;jloc++) {
			tmp1=PI*(jloc-sysparams.Mz)/((sysparams.Mz-1)*sysparams.az);
			spectz[jloc]=2.0/MASS_z*tmp1*tmp1;
		};
	#endif
	#ifdef DiscreteLaplacianSpectrum
		// Discrete Laplacian spectrum
		tmp1=1.0/(MASS_x*sysparams.ax*sysparams.ax);
		tmp2=2.0*PI/(sysparams.Mx);
		for(jloc=0;jloc<sysparams.Mx;jloc++){
			spectx[jloc]=tmp1*(1.0-cos(tmp2*jloc));
		};
		tmp1=1.0/(MASS_y*sysparams.ay*sysparams.ay);
		tmp2=2.0*PI/(sysparams.My);
		for(jloc=0;jloc<sysparams.My;jloc++){
			specty[jloc]=tmp1*(1.0-cos(tmp2*jloc));
		};
		tmp1=1.0/(MASS_z*sysparams.az*sysparams.az);
		tmp2=2.0*PI/(sysparams.Mz);
		for(jloc=0;jloc<sysparams.Mz;jloc++){
			spectz[jloc]=tmp1*(1.0-cos(tmp2*jloc));
		};
	#endif
	
	tmp1=1.0/((double) sysparams.M);
	fftw_execute_dft(FFTW_Variables.plan,Psis.Psi1,Psis.Psi1);
	fftw_execute_dft(FFTW_Variables.plan,Psis.Psi2,Psis.Psi2);
	fftw_execute_dft(FFTW_Variables.plan,Psis.Psi3,Psis.Psi3);
	for(jx=0,*kinEnergy1=0.0,*kinEnergy2=0.0,*kinEnergy3=0.0;jx<sysparams.Mx;jx++){
		for(jy=0;jy<sysparams.My;jy++) {
			for(jz=0;jz<sysparams.Mz;jz++){
				jloc=jz+sysparams.Mz*(jy+sysparams.My*jx);
				*kinEnergy1+=sqr_norm(Psis.Psi1[jloc])*(spectx[jx]+specty[jy]+spectz[jz])*tmp1;
				*kinEnergy2+=sqr_norm(Psis.Psi2[jloc])*(spectx[jx]+specty[jy]+spectz[jz])*tmp1;
				*kinEnergy3+=sqr_norm(Psis.Psi3[jloc])*(spectx[jx]+specty[jy]+spectz[jz])*tmp1;
				Psis.Psi1[jloc][0]*=tmp1;
				Psis.Psi1[jloc][1]*=tmp1;
				Psis.Psi2[jloc][0]*=tmp1;
				Psis.Psi2[jloc][1]*=tmp1;
				Psis.Psi3[jloc][0]*=tmp1;
				Psis.Psi3[jloc][1]*=tmp1;
			};
		};
	};
	fftw_execute_dft(FFTW_Variables.inverseplan,Psis.Psi1,Psis.Psi1);
	fftw_execute_dft(FFTW_Variables.inverseplan,Psis.Psi2,Psis.Psi2);
	fftw_execute_dft(FFTW_Variables.inverseplan,Psis.Psi3,Psis.Psi3);
}

void potentialEnergy(double *potEnergy,GPE_Psis Psis,System_Parameters sysparams,GPE_Parameters gpe) {
	register int jloc,jx,jy,jz;
	register double SqrtMinusDmu;
	register dcomplex PsiVec[3];
	for(jx=0,*potEnergy=0.0;jx<sysparams.Mx;jx++){
		for(jy=0;jy<sysparams.My;jy++) {
			for(jz=0;jz<sysparams.Mz;jz++){
				jloc=jz+sysparams.Mz*(jy+sysparams.My*jx);
				SqrtMinusDmu=sqrt(sysparams.trapx*(jx-0.5*(sysparams.Mx-1.0))*(jx-0.5*(sysparams.Mx-1.0))
				                 +sysparams.trapy*(jy-0.5*(sysparams.My-1.0))*(jy-0.5*(sysparams.My-1.0))
				                 +sysparams.trapz*(jz-0.5*(sysparams.Mz-1.0))*(jz-0.5*(sysparams.Mz-1.0)));
				PsiVec[0]=Psis.Psi1[jloc];
				PsiVec[1]=Psis.Psi2[jloc];
				PsiVec[2]=Psis.Psi3[jloc];
				*potEnergy+=(evaluateSpline(0,SqrtMinusDmu,PsiVec,gpe,sysparams)).re;
			};
		};
	};
}

void angularMomenta(double *Lx,double *Ly,double *Lz,GPE_Psis Psis,System_Parameters sysparams) {
	register int i;
	register int jloc,jx,jy,jz;
	register int Dx=sysparams.My*sysparams.Mz,Dy=sysparams.Mz,Dz=1;
	double x[sysparams.Mx],y[sysparams.My],z[sysparams.Mz];
	register double tmpax,tmpay,tmpaz;
	fftw_complex *whichPsi[3];
	
	whichPsi[0]=Psis.Psi1; whichPsi[1]=Psis.Psi2; whichPsi[2]=Psis.Psi3;
	tmpax=1.0/sysparams.ax; tmpay=1.0/sysparams.ay; tmpaz=1.0/sysparams.az; 
	for(jloc=0;jloc<sysparams.Mx;jloc++)
		x[jloc]=(((double) jloc)-0.5*(sysparams.Mx-1.0))*sysparams.ax;
	for(jloc=0;jloc<sysparams.My;jloc++)
		y[jloc]=(((double) jloc)-0.5*(sysparams.My-1.0))*sysparams.ay;
	for(jloc=0;jloc<sysparams.Mz;jloc++)
		z[jloc]=(((double) jloc)-0.5*(sysparams.Mz-1.0))*sysparams.az;
	for(jx=0,*Lx=0.0;jx<sysparams.Mx;jx++) {
		for(jy=0;jy<sysparams.My-1;jy++) {
			for(jz=0;jz<sysparams.Mz-1;jz++) {
				jloc=jz+sysparams.Mz*(jy+jx*sysparams.My);
				for(i=0;i<3;i++) {
					*Lx+=y[jy]*((whichPsi[i][jloc+Dz][1]-whichPsi[i][jloc][1])*(whichPsi[i][jloc+Dz][0]+whichPsi[i][jloc][0])-
					            (whichPsi[i][jloc+Dz][0]-whichPsi[i][jloc][0])*(whichPsi[i][jloc+Dz][1]+whichPsi[i][jloc][1]))*tmpaz*0.5-
					     z[jz]*((whichPsi[i][jloc+Dy][1]-whichPsi[i][jloc][1])*(whichPsi[i][jloc+Dy][0]+whichPsi[i][jloc][0])-
					            (whichPsi[i][jloc+Dy][0]-whichPsi[i][jloc][0])*(whichPsi[i][jloc+Dy][1]+whichPsi[i][jloc][1]))*tmpay*0.5;
				};
			};
		};
	};
	for(jy=0,*Ly=0.0;jy<sysparams.My;jy++) {
		for(jz=0;jz<sysparams.Mz-1;jz++) {
			for(jx=0;jx<sysparams.Mx-1;jx++) {
				jloc=jz+sysparams.Mz*(jy+jx*sysparams.My);
				for(i=0;i<3;i++) {
					*Ly+=z[jz]*((whichPsi[i][jloc+Dx][1]-whichPsi[i][jloc][1])*(whichPsi[i][jloc+Dx][0]+whichPsi[i][jloc][0])-
					            (whichPsi[i][jloc+Dx][0]-whichPsi[i][jloc][0])*(whichPsi[i][jloc+Dx][1]+whichPsi[i][jloc][1]))*tmpax*0.5-
					     x[jx]*((whichPsi[i][jloc+Dz][1]-whichPsi[i][jloc][1])*(whichPsi[i][jloc+Dz][0]+whichPsi[i][jloc][0])-
					            (whichPsi[i][jloc+Dz][0]-whichPsi[i][jloc][0])*(whichPsi[i][jloc+Dz][1]+whichPsi[i][jloc][1]))*tmpaz*0.5;
				};
			};
		};
	};
	for(jz=0,*Lz=0.0;jz<sysparams.Mz;jz++) {
		for(jx=0;jx<sysparams.Mx-1;jx++) {
			for(jy=0;jy<sysparams.My-1;jy++) {
				jloc=jz+sysparams.Mz*(jy+jx*sysparams.My);
				for(i=0;i<3;i++) {
					*Lz+=x[jx]*((whichPsi[i][jloc+Dy][1]-whichPsi[i][jloc][1])*(whichPsi[i][jloc+Dy][0]+whichPsi[i][jloc][0])-
					            (whichPsi[i][jloc+Dy][0]-whichPsi[i][jloc][0])*(whichPsi[i][jloc+Dy][1]+whichPsi[i][jloc][1]))*tmpay*0.5-
					     y[jy]*((whichPsi[i][jloc+Dx][1]-whichPsi[i][jloc][1])*(whichPsi[i][jloc+Dx][0]+whichPsi[i][jloc][0])-
					            (whichPsi[i][jloc+Dx][0]-whichPsi[i][jloc][0])*(whichPsi[i][jloc+Dx][1]+whichPsi[i][jloc][1]))*tmpax*0.5;
				};
			};
		};
	};
}

void norms(double *normPsi1,double *normPsi2,double *normPsi3,GPE_Psis Psis,System_Parameters sysparams) {
	register int j;  
	for(*normPsi1=0.0,*normPsi2=0.0,*normPsi3=0.0,j=0;j<sysparams.M;j++) {
		*normPsi1+=sqr_norm(Psis.Psi1[j]);
		*normPsi2+=sqr_norm(Psis.Psi2[j]);
		*normPsi3+=sqr_norm(Psis.Psi3[j]);
	};
	*normPsi1=sqrt(*normPsi1);
	*normPsi2=sqrt(*normPsi2);
	*normPsi3=sqrt(*normPsi3);
}

dcomplex PontrjaginIndex(GPE_Psis Psis,System_Parameters sysparams) {
	register int jx,jy,jz,jloc;
	register int i1,i2,i3;
	register int Dx=sysparams.My*sysparams.Mz,Dy=sysparams.Mz,Dz=1;
	register dcomplex W=0.0;
	register dcomplex derivatives[3][3];
	register double epsilon[3][3][3];
	register double normPsi1,normPsi2,normPsi3;
	register double normalize;
	register dcomplex tmp;
	
	norms(&normPsi1,&normPsi2,&normPsi3,Psis,sysparams);
	normPsi1/=sqrt(sysparams.M);
	normPsi2/=sqrt(sysparams.M);
	normPsi3/=sqrt(sysparams.M);
	cout<<normPsi1+normPsi2+normPsi3<<endl;
	
	// Fill up Levi-Civita symbol
	for(i1=0;i1<3;i1++) {
		for(i2=0;i2<3;i2++) {
			for(i3=0;i3<3;i3++) {
				epsilon[i1][i2][i3]=0.0;
			};
		};
	};
	epsilon[0][1][2]= 1.0; epsilon[0][2][1]=-1.0;
	epsilon[1][2][0]= 1.0; epsilon[1][0][2]=-1.0;
	epsilon[2][0][1]= 1.0; epsilon[2][1][0]=-1.0;
	
	// Calculate Pontrjagin index
	for(jx=1;jx<(sysparams.Mx-1);jx++) {
		for(jy=1;jy<(sysparams.My-1);jy++) {
			for(jz=1;jz<(sysparams.Mz-1);jz++) {
				jloc=jz+sysparams.Mz*(jy+jx*sysparams.My);
				normalize=sqrt(sqr_norm(Psis.Psi1[jloc])+sqr_norm(Psis.Psi2[jloc])+sqr_norm(Psis.Psi3[jloc]));
				if(normalize>0.00001) {
					normalize=1.0/normalize;
					derivatives[0][0].re=(Psis.Psi1[jloc+Dx][0]-Psis.Psi1[jloc][0])*normalize;
					derivatives[0][0].im=(Psis.Psi1[jloc+Dx][1]-Psis.Psi1[jloc][1])*normalize;
					derivatives[1][0].re=(Psis.Psi1[jloc+Dy][0]-Psis.Psi1[jloc][0])*normalize;
					derivatives[1][0].im=(Psis.Psi1[jloc+Dy][1]-Psis.Psi1[jloc][1])*normalize;
					derivatives[2][0].re=(Psis.Psi1[jloc+Dz][0]-Psis.Psi1[jloc][0])*normalize;
					derivatives[2][0].im=(Psis.Psi1[jloc+Dz][1]-Psis.Psi1[jloc][1])*normalize;
					derivatives[0][1].re=(Psis.Psi2[jloc+Dx][0]-Psis.Psi2[jloc][0])*normalize;
					derivatives[0][1].im=(Psis.Psi2[jloc+Dx][1]-Psis.Psi2[jloc][1])*normalize;
					derivatives[1][1].re=(Psis.Psi2[jloc+Dy][0]-Psis.Psi2[jloc][0])*normalize;
					derivatives[1][1].im=(Psis.Psi2[jloc+Dy][1]-Psis.Psi2[jloc][1])*normalize;
					derivatives[2][1].re=(Psis.Psi2[jloc+Dz][0]-Psis.Psi2[jloc][0])*normalize;
					derivatives[2][1].im=(Psis.Psi2[jloc+Dz][1]-Psis.Psi2[jloc][1])*normalize;
					derivatives[0][2].re=(Psis.Psi3[jloc+Dx][0]-Psis.Psi3[jloc][0])*normalize;
					derivatives[0][2].im=(Psis.Psi3[jloc+Dx][1]-Psis.Psi3[jloc][1])*normalize;
					derivatives[1][2].re=(Psis.Psi3[jloc+Dy][0]-Psis.Psi3[jloc][0])*normalize;
					derivatives[1][2].im=(Psis.Psi3[jloc+Dy][1]-Psis.Psi3[jloc][1])*normalize;
					derivatives[2][2].re=(Psis.Psi3[jloc+Dz][0]-Psis.Psi3[jloc][0])*normalize;
					derivatives[2][2].im=(Psis.Psi3[jloc+Dz][1]-Psis.Psi3[jloc][1])*normalize;
					tmp=PontrAuxiliaryFunction(epsilon,derivatives);
					W.re+=tmp.re; W.im+=tmp.im;
				}
			};
		};
	};
	W.re/=(8.0*PI); W.im/=(8.0*PI);
	return W;
}

dcomplex mult_dcomplex(dcomplex a,dcomplex b) {
	return a*b;
}

dcomplex PontrAuxiliaryFunction(double epsilon[3][3][3],dcomplex derivatives[3][3]) {
	register int i1,i2,i3;
	register int alpha1,alpha2,alpha3;
	register dcomplex tmp;
	register dcomplex Wtmp=0.0;
	for(i1=0;i1<3;i1++) {
		for(i2=0;i2<3;i2++) {
			for(i3=0;i3<3;i3++) {
				for(alpha1=0;alpha1<3;alpha1++) {
					for(alpha2=0;alpha2<3;alpha2++) {
						for(alpha3=0;alpha3<3;alpha3++) {
							tmp=mult_dcomplex(mult_dcomplex(derivatives[i1][alpha1],derivatives[i2][alpha2]),derivatives[i3][alpha3]);
							tmp.re*=epsilon[i1][i2][i3]*epsilon[alpha1][alpha2][alpha3];
							tmp.im*=epsilon[i1][i2][i3]*epsilon[alpha1][alpha2][alpha3];
							Wtmp.re+=tmp.re; Wtmp.im+=tmp.im;
						};
					};
				};
			};
		};
	};
	return Wtmp;
}