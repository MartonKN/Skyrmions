#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include "dcomplex.h"
#include "matrix.h"
using namespace std;

// random number generator seed
long idum_matrix=-13514;



// CLASS DMATRIX: REAL MATRICES
// Constructors
dmatrix::dmatrix(void) {
	m1=0;
	m2=0;
	pm=NULL;
}

dmatrix::dmatrix(int n1,int n2) {
	register int i1,i2;
	m1=n1;
	m2=n2;
	pm=new double* [n1];
	for(i1=0;i1<n1;i1++)
		pm[i1]=new double [n2];
}

dmatrix::dmatrix(dmatrix& A) {
	register int i1,i2;
	m1=A.m1;
	m2=A.m2;
	pm=new double* [m1];
	for(i1=0;i1<m1;i1++) {
		pm[i1]=new double [m2];
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]=A.pm[i1][i2];
	};
};

dmatrix::dmatrix(double& d) {
	m1=m2=1;
	pm=new double* [1];
	pm[0]=new double[1];
	pm[0][0]=d;
}

dmatrix::dmatrix(int n1,int n2,double** A){
	register int i1,i2;
	m1=n1; m2=n2;
	pm=new double* [m1];
	for(i1=0;i1<m1;i1++) {
		pm[i1]=new double [m2];
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]=A[i1][i2];
	};
}


// Functions
void dmatrix::size(int& n1,int& n2) {
	n1=m1; n2=m2;
}

double dmatrix::max(int& i1,int& i2) {
	if(!pm) {
		cerr<<endl<<"You wanted to take the greatest element of an empty matrix."<<endl;
		return NAN;
	}
	register double maxValue;
	register int j1,j2;
	i1=0;i2=0; maxValue=pm[0][0];
	for(j1=0;j1<m1;j1++) {
		for(j2=0;j2<m2;j2++) {
			if (maxValue<pm[j1][j2]) {
				maxValue=pm[j1][j2];
				i1=j1; i2=j2;
			}
		};
	};
	return maxValue;
}
  
double dmatrix::min(int& i1,int& i2) {
	if(!pm) {
		cerr<<endl<<"You wanted to take the smallest element of an empty matrix."<<endl;
		return NAN;
	}
	register double minValue;
	register int j1,j2;
	i1=0;i2=0; minValue=pm[0][0];
	for(j1=0;j1<m1;j1++) {
		for(j2=0;j2<m2;j2++) {
			if (minValue>pm[j1][j2]) {
				minValue=pm[j1][j2];
				i1=j1; i2=j2;
			}
		};
	};
	return minValue;
}

double dmatrix::maxAbs(int& i1,int& i2) {
	if(!pm) {
		cerr<<endl<<"You wanted to take the greatest absolute value element of an empty matrix."<<endl;
		return NAN;
	}
	register double maxValue;
	register int j1,j2;
	i1=0;i2=0; maxValue=fabs(pm[0][0]);
	for(j1=0;j1<m1;j1++) {
		for(j2=0;j2<m2;j2++) {
			if (maxValue<fabs(pm[j1][j2])) {
				maxValue=fabs(pm[j1][j2]);
				i1=j1; i2=j2;
			}
		};
	};
	return maxValue;
}

double dmatrix::minAbs(int& i1,int& i2) {
	if(!pm) {
		cerr<<endl<<"You wanted to take the smallest absolute value element of an empty matrix."<<endl;
		return NAN;
	}
	register double minValue;
	register int j1,j2;
	i1=0;i2=0; minValue=fabs(pm[0][0]);
	for(j1=0;j1<m1;j1++) {
		for(j2=0;j2<m2;j2++) {
			if (minValue>fabs(pm[j1][j2])) {
				minValue=fabs(pm[j1][j2]);
				i1=j1; i2=j2;
			}
		};
	};
	return minValue;
}

double dmatrix::max(void) {
	int i1,i2;
	return (*this).max(i1,i2);
}

double dmatrix::min(void) {
	int i1,i2;
	return (*this).min(i1,i2);
}

double dmatrix::maxAbs(void) {
	int i1,i2;
	return (*this).maxAbs(i1,i2);
}

double dmatrix::minAbs(void) {
	int i1,i2;
	return (*this).minAbs(i1,i2);
}

double dmatrix::norm(void) {
	register double norm=0.0;
	register int j1,j2;
	for(j1=0;j1<m1;j1++) {
		for(j2=0;j2<m2;j2++)
			norm+=pm[j1][j2]*pm[j1][j2];
	};
	norm=sqrt(norm);
	return norm;
}

dmatrix dmatrix::transpose(void) {
	register int i1,i2;
	if(m1!=m2) {
		double** p=pm;
		i1=m1; m1=m2; m2=i1;
		pm=new double* [m1];
		for(i1=0;i1<m1;i1++) {
			pm[i1]=new double [m2];
			for(i2=0;i2<m2;i2++)
				pm[i1][i2]=p[i2][i1];
		};
		for(i1=0;i1<m2;i1++)
			delete [] p[i1];
		delete [] p;
	}
	else {
		register double tmp;
		for(i1=0;i1<m1;i1++) {
			for(i2=0;i2<i1;i2++) {
				tmp=pm[i1][i2];
				pm[i1][i2]=pm[i2][i1];
				pm[i2][i1]=tmp;
			};
		};
	}
	return *this;
}

double dmatrix::trace(void) {
	register int i;
	register double tr=0.0;
	if(m1!=m2) {
		cerr<<endl<<"You took the trace of a non-rectangular matrix!"<<endl;
	}
	for(i=0;i<(m1<m2?m1:m2);i++)
		tr+=pm[i][i];
	return tr;
}

dmatrix dmatrix::diagonal(void) {
	register int j1,j2;
	if(m1==1) {
		dmatrix result(m2,m2);
		for(j1=0;j1<m2;j1++) {
			for(j2=0;j2<m2;j2++) {
				result.pm[j1][j2]=0.0;
			};
			result.pm[j1][j1]=pm[0][j1];
		};
		return result;
	}
	else if(m2==1) {
		dmatrix result(m1,m1);
		for(j1=0;j1<m1;j1++) {
			for(j2=0;j2<m1;j2++) {
				result.pm[j1][j2]=0.0;
			};
			result.pm[j1][j1]=pm[j1][0];
		};
		return result;
	}
	else {
		int m=(m1<m2 ? m1 : m2);
		dmatrix result(m,1);
		for(j1=0;j1<m;j1++)
			result.pm[j1][0]=pm[j1][j1];
		return result;
	}
}

dmatrix dmatrix::randomize(void) {
	register double ran;
	register int i1,i2;
	
	// initialize the random number generator
	time_t now=time(0);
	tm *ltime=localtime(&now);
	int secs=1+ltime->tm_sec;
	for(i1=0;i1<secs;i1++)
		ran=ran_matrix(&idum_matrix);
	
	// fill matrix with random numbers
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]=ran_matrix(&idum_matrix);
	};
	return *this;
}


// Operators
double* dmatrix::operator [] (int i) {
	return pm[i];
}

dmatrix dmatrix::operator = (const dmatrix& M) {
	register int i1,i2;
	if(M.m1!=m1 || M.m2!=m2) {
		for(i1=0;i1<m1;i1++) {
			delete [] pm[i1];
		};
		delete [] pm;
		m1=M.m1; m2=M.m2;
		pm=new double* [m1];
		for(i1=0;i1<m1;i1++) {
			pm[i1]=new double [m2];
		};
	}
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++) {
			pm[i1][i2]=M.pm[i1][i2];
		};
	};
	return *this;
}

dmatrix dmatrix::operator = (const double& d) {
	int i1;
	if(m1!=1 || m2!=1) {
		for(i1=0;i1<m1;i1++) {
			delete [] pm[i1];
		};
		delete [] pm;
		m1=1; m2=1;
		pm=new double* [1];
		pm[0]=new double [1];
	}
	pm[0][0]=d;
	return *this;  
}

dmatrix dmatrix::operator += (const dmatrix& B) {
	register int i1,i2;
	if(B.m1==1 && B.m2==1) {
		(*this)+=B.pm[0][0];
		return *this;
	}
	else if(m1==1 && m2==1) {
		double tmp=pm[0][0];
		*this=B;
		*this+=tmp;
		return *this;
	}
	if(m1!=B.m1 || m2!=B.m2) {
		cerr<<endl<<"You can't add two matrices of non-matching sizes. Nothing has been done."<<endl;
		return *this;
	}
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]+=B.pm[i1][i2];
	};
	return *this;
}

dmatrix dmatrix::operator += (const double& d) {
	register int i1,i2;
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]+=d;
	};
	return *this;
}

dmatrix dmatrix::operator -= (const dmatrix& B) {
	register int i1,i2;
	if(B.m1==1 && B.m2==1) {
		(*this)-=B.pm[0][0];
		return *this;
	}
	else if(m1==1 && m2==1) {
		double tmp=pm[0][0];
		*this=B;
		*this-=tmp;
		return *this;
	}
	if(m1!=B.m1 || m2!=B.m2) {
		cerr<<endl<<"You can't subtract two matrices of non-matching sizes. Nothing has been done."<<endl;
		return *this;
	}
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]-=B.pm[i1][i2];
	};
	return *this;
}

dmatrix dmatrix::operator -= (const double& d) {
	register int i1,i2;
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]-=d;
	};
	return *this;
}

dmatrix dmatrix::operator *= (const dmatrix& B) {
	register int i1,i2,i3;
	if(B.m1==1 && B.m2==1) {
		(*this)*=B.pm[0][0];
		return *this;
	}
	else if(m1==1 && m2==1) {
		double tmp=pm[0][0];
		*this=B;
		(*this)*=tmp;
		return *this;
	}
	else {
		if(m2!=B.m1) {
			cerr<<endl<<"You can't multiply two matrices of non-matching sizes. Nothing has been done."<<endl;
			return *this;
		}
		double** p;
		p=new double* [m1];
		for(i1=0;i1<m1;i1++)
			p[i1]=new double [B.m2];
		for(i1=0;i1<m1;i1++) {
			for(i2=0;i2<B.m2;i2++) {
				for(i3=0,p[i1][i2]=0.0;i3<m2;i3++)
					p[i1][i2]+=(pm[i1][i3]*B.pm[i3][i2]);
			};
		};
		for(i1=0;i1<m1;i1++)
			delete [] pm[i1];
		delete [] pm;
		pm=p;
		m2=B.m2;
		return *this;
	}
}

dmatrix dmatrix::operator *= (const double& d) {
	register int i1,i2;
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]*=d;
	};
	return *this;
}

dmatrix operator + (const dmatrix& A,const dmatrix& B) {
	register int i1,i2;
	dmatrix result;
	if(A.m1==1 && A.m2==1) {
		result.m1=B.m1;
		result.m2=B.m2;
		result.pm=new double* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new double [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=B.pm[i1][i2]+A.pm[0][0];
		};
	}
	else if(B.m1==1 && B.m2==1) {
		result.m1=A.m1;
		result.m2=A.m2;
		result.pm=new double* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new double [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A.pm[i1][i2]+B.pm[0][0];
		};
	}
	else {
		if(A.m1!=B.m1 || A.m2!=B.m2) {
			cerr<<endl<<"You can't add two matrices of different size."<<endl;
			return result;
		}
		else {
			result.m1=A.m1;
			result.m2=A.m2;
			result.pm=new double* [result.m1];
			for(i1=0;i1<result.m1;i1++) {
				result.pm[i1]=new double [result.m2];
				for(i2=0;i2<result.m2;i2++)
					result.pm[i1][i2]=A.pm[i1][i2]+B.pm[i1][i2];
			};
		}
	}
	return result;
}

dmatrix operator - (const dmatrix& A,const dmatrix& B) {
	register int i1,i2;
	dmatrix result;
	if(A.m1==1 && A.m2==1) {
		result.m1=B.m1;
		result.m2=B.m2;
		result.pm=new double* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new double [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A.pm[0][0]-B.pm[i1][i2];
		};
	}
	else if(B.m1==1 && B.m2==1) {
		result.m1=A.m1;
		result.m2=A.m2;
		result.pm=new double* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new double [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A.pm[i1][i2]-B.pm[0][0];
		};
	}
	else {
		if(A.m1!=B.m1 || A.m2!=B.m2) {
			cerr<<endl<<"You can't add two matrices of different size."<<endl;
			return result;
		}
		else {
			result.m1=A.m1;
			result.m2=A.m2;
			result.pm=new double* [result.m1];
			for(i1=0;i1<result.m1;i1++) {
				result.pm[i1]=new double [result.m2];
				for(i2=0;i2<result.m2;i2++)
					result.pm[i1][i2]=A.pm[i1][i2]-B.pm[i1][i2];
			};
		}
	}
	return result;
}

dmatrix operator - (const dmatrix& A,const double& d) {
	register int i1,i2;
	dmatrix result(A.m1,A.m2);
	for(i1=0;i1<A.m1;i1++) {
		for(i2=0;i2<A.m2;i2++)
			result.pm[i1][i2]=A.pm[i1][i2]-d;
	};
	return result;
}

dmatrix operator - (const double& d,const dmatrix& A) {
	register int i1,i2;
	dmatrix result(A.m1,A.m2);
	for(i1=0;i1<A.m1;i1++) {
		for(i2=0;i2<A.m2;i2++)
			result.pm[i1][i2]=d-A.pm[i1][i2];
	};
	return result;
}

dmatrix operator * (const dmatrix& A,const dmatrix& B) {
	register int i1,i2;
	dmatrix result;
	if(A.m1==1 && A.m2==1) {
		result.m1=B.m1;
		result.m2=B.m2;
		result.pm=new double* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new double [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=B.pm[i1][i2]*A.pm[0][0];
		};
	}
	else if(B.m1==1 && B.m2==1) {
		result.m1=A.m1;
		result.m2=A.m2;
		result.pm=new double* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new double [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A.pm[i1][i2]*B.pm[0][0];
		};
	}
	else {
		if(A.m2!=B.m1) {
			cerr<<endl<<"You can't multiply two matrices of non-matching sizes."<<endl;
			return result;
		}
		else {
			register int i3;
			result.m1=A.m1;
			result.m2=B.m2;
			result.pm=new double* [result.m1];
			for(i1=0;i1<result.m1;i1++) {
				result.pm[i1]=new double [result.m2];
				for(i2=0;i2<result.m2;i2++) {
					for(i3=0,result.pm[i1][i2]=0.0;i3<A.m2;i3++)
						result.pm[i1][i2]+=A.pm[i1][i3]*B.pm[i3][i2];
				};
			};
		}
	}
	return result;
}

dmatrix operator * (const dmatrix& A,const double& d) {
	register int i1,i2;
	dmatrix result(A.m1,A.m2);
	for(i1=0;i1<A.m1;i1++) {
		for(i2=0;i2<A.m2;i2++) {
			result.pm[i1][i2]=A.pm[i1][i2]*d;
		};
	};
	return result;
}

dmatrix operator * (const double& d,const dmatrix& A) {
	register int i1,i2;
	dmatrix result(A.m1,A.m2);
	for(i1=0;i1<A.m1;i1++) {
		for(i2=0;i2<A.m2;i2++) {
			result.pm[i1][i2]=A.pm[i1][i2]*d;
		};
	};
	return result;
}

// IO
ostream& operator << (ostream& os,const dmatrix& M) {
	int array_width=12;
	int precision=5;
	int i1,i2;
	cout<<setiosflags(ios::left);
	cout<<setprecision(precision);
	cout<<fixed;
	cout<<endl;
	for(i1=0;i1<M.m1;i1++) {
		for(i2=0;i2<M.m2;i2++)
			if(fabs(M.pm[i1][i2])>=1e5) {cout<<scientific;}
	};
	for(i1=0;i1<M.m1;i1++) {
		for(i2=0;i2<M.m2;i2++)
			cout<<"\t"<<(M.pm[i1][i2]<0 ? '-' : ' ')<<setw(array_width)<<fabs(M.pm[i1][i2]);
		cout<<endl;
	};
	cout.unsetf(ios::floatfield);
	return os;
}

// LU decomposition (See Numerical Recipies in C, p. 46)
dmatrix dmatrix::LU(int *indx,double *d)
// Gives back the LU decomposition of the matix. (The matrix itself will remain untouched,) 
// indx and d are auxiliary variables used in calculation of the inverse, the determinant and solving linear equations.
{
	if(m1!=m2) {
		cerr<<endl<<"You can't make the LU decomposition of a non-square matrix. Nothing has been done."<<endl;
		return *this;
	}
	else {
		int n=m1;
		register int j1,j2;
		register double** tmpMx;
		tmpMx=new double* [n+1];
		tmpMx[0]=new double [n+1];
		for(j1=1;j1<=n;j1++) {
			tmpMx[j1]=new double [n+1];
			for(j2=1;j2<=n;j2++) {
				tmpMx[j1][j2]=pm[j1-1][j2-1];
			};
		};
		
		ludcmp(tmpMx,n,indx,d);
		
		dmatrix result(n,n);
		delete [] tmpMx[0];
		for(j1=1;j1<=n;j1++) {
			for(j2=1;j2<=n;j2++) {
				result.pm[j1-1][j2-1]=tmpMx[j1][j2];
			};
			delete [] tmpMx[j1];
		};
		delete [] tmpMx;
		return result;
	}
}

dmatrix dmatrix::LU(void) 
// Same as the previous one. It just doesn't keep track of indx and d.
{
	int* indx;
	indx=new int [m1];
	double d;
	dmatrix result;
	result=(*this).LU(indx,&d);
	delete [] indx;
	return result;
}

dmatrix dmatrix::LU_backsub(int *indx, dmatrix b)
{
	if(m1!=m2 || m2!=b.m1 || b.m2!=1) {
		cerr<<endl<<"Inappropriate size of the matrix or the vector in function LU_backsub.";
		cerr<<"Nothing has been done."<<endl;
		dmatrix emptyMx;
		return emptyMx;
	}
	else {
		int n=m1;
		register double* b_;
		b_=new double [m1+1];
		register int j1,j2;
		for(j1=1;j1<=n;j1++)
			b_[j1]=b.pm[j1-1][0];
		register double** tmpMx;
		tmpMx=new double* [n+1];
		tmpMx[0]=new double [n+1];
		for(j1=1;j1<=n;j1++) {
			tmpMx[j1]=new double [n+1];
			for(j2=1;j2<=n;j2++) {
				tmpMx[j1][j2]=pm[j1-1][j2-1];
			};
		};
		
		lubksb(tmpMx,n,indx,b_);
		
		dmatrix result(b.m1,b.m2);
		for(j1=1;j1<=n;j1++) {
			result.pm[j1-1][0]=b_[j1];
			delete [] tmpMx[j1];
		};
		delete [] tmpMx[0];
		delete [] tmpMx;
		delete [] b_;
		return result;
	}
}





// Linear equation solver
dmatrix linearSolve(dmatrix& A,dmatrix& b) {
	// Solves the linear equation A*x=b;
	if(A.m1!=A.m2 || A.m2!=b.m1 || b.m2!=1) {
		cerr<<endl<<"Inappropriate size of the matrix or the vector in function linearSolve. ";
		cerr<<"Nothing has been done."<<endl;
		dmatrix x;
		return x;
	}
	else {
		dmatrix x(b.m1,b.m2);
		dmatrix tmpA(A.m1,A.m2);
		int* indx;
		indx=new int [A.m1];
		double d;
		tmpA=A.LU(indx,&d);
		x=tmpA.LU_backsub(indx,b);
		delete [] indx;
		return x;
	}
}

dmatrix dmatrix::invert(void) {
	if(m1!=m2) {
		cerr<<endl<<"You can't take the inverse of a non-square matrix. Nothing has been done."<<endl;
		return *this;
	}
	else {
		int n=m1;
		register int j1,j2;
		register double** tmpMx;
		tmpMx=new double* [n+1];
		tmpMx[0]=new double [n+1];
		for(j1=1;j1<=n;j1++) {
			tmpMx[j1]=new double [n+1];
			for(j2=1;j2<=n;j2++)
				tmpMx[j1][j2]=pm[j1-1][j2-1];
		};
		int* indx; indx=new int [n+1];
		double* col; col=new double [n+1];
		double d;
		
		ludcmp(tmpMx,n,indx,&d);
		for(j1=1;j1<=n;j1++) {
			for(j2=1;j2<=n;j2++) 
				col[j2]=0.0;
			col[j1]=1.0;
			lubksb(tmpMx,n,indx,col);
			for(j2=1;j2<=n;j2++) 
				pm[j2-1][j1-1]=col[j2];
		};
		delete [] indx;
		delete [] col;
		for(j1=0;j1<=n;j1++)
			delete [] tmpMx[j1];
		delete [] tmpMx;
		return *this;
	}
}

double dmatrix::determinant(void) {
	if(m1!=m2) {
		cerr<<endl<<"You can't take the determinant of a non-square matrix. Nothing has been done."<<endl;
		return NAN;
	}
	else {
		register double det;
		register int* indx;
		indx=new int [m1];
		dmatrix tmpMx(m1,m2);
		register int j;
		
		tmpMx=(*this).LU(indx,&det);
		for(j=0;j<m1;j++)
			det*=tmpMx[j][j];
		delete [] indx;
		return det;
	}
}


// CLASS CMATRIX: COMPLEX MATRICES
// Constructors
cmatrix::cmatrix(void) {
	m1=0;
	m2=0;
	pm=NULL;
}

cmatrix::cmatrix(int n1,int n2) {
	register int i1,i2;
	m1=n1;
	m2=n2;
	pm=new dcomplex* [n1];
	for(i1=0;i1<n1;i1++)
		pm[i1]=new dcomplex [n2];
}

cmatrix::cmatrix(cmatrix& A) {
	register int i1,i2;
	m1=A.m1;
	m2=A.m2;
	pm=new dcomplex* [m1];
	for(i1=0;i1<m1;i1++) {
		pm[i1]=new dcomplex [m2];
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]=A.pm[i1][i2];
	};
};

cmatrix::cmatrix(dmatrix& A) {
	register int i1,i2;
	A.size(m1,m2);
	pm=new dcomplex* [m1];
	for(i1=0;i1<m1;i1++) {
		pm[i1]=new dcomplex [m2];
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]=A[i1][i2];
	};
};

cmatrix::cmatrix(dcomplex& c) {
	m1=m2=1;
	pm=new dcomplex* [1];
	pm[0]=new dcomplex[1];
	pm[0][0]=c;
}

cmatrix::cmatrix(double& d) {
	m1=m2=1;
	pm=new dcomplex* [1];
	pm[0]=new dcomplex[1];
	pm[0][0]=d;
}

cmatrix::cmatrix(int n1,int n2,dcomplex** A){
	register int i1,i2;
	m1=n1; m2=n2;
	pm=new dcomplex* [m1];
	for(i1=0;i1<m1;i1++) {
		pm[i1]=new dcomplex [m2];
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]=A[i1][i2];
	};
}

cmatrix::cmatrix(int n1,int n2,double** A){
	register int i1,i2;
	m1=n1; m2=n2;
	pm=new dcomplex* [m1];
	for(i1=0;i1<m1;i1++) {
		pm[i1]=new dcomplex [m2];
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]=A[i1][i2];
	};
}

// Functions
void cmatrix::size(int& n1,int& n2) {
	n1=m1; n2=m2;
}

cmatrix cmatrix::re(void) {
	register int j1,j2;
	for(j1=0;j1<m1;j1++) {
		for(j2=0;j2<m2;j2++)
			pm[j1][j2]=pm[j1][j2].re;
	};
	return *this;
}

cmatrix cmatrix::im(void) {
	register int j1,j2;
	for(j1=0;j1<m1;j1++) {
		for(j2=0;j2<m2;j2++)
			pm[j1][j2]=pm[j1][j2].im;
	};
	return *this;
}

double cmatrix::maxAbs(int& i1,int& i2) {
	if(!pm) {
		cerr<<endl<<"You wanted to take the greatest absolute value element of an empty matrix."<<endl;
		return NAN;
	}
	register double maxValue,tmp;
	register int j1,j2;
	i1=0;i2=0; maxValue=pm[0][0].re*pm[0][0].re+pm[0][0].im*pm[0][0].im;
	for(j1=0;j1<m1;j1++) {
		for(j2=0;j2<m2;j2++) {
			tmp=pm[j1][j2].re*pm[j1][j2].re+pm[j1][j2].im*pm[j1][j2].im;
			if (maxValue<tmp) {
				maxValue=tmp;
				i1=j1; i2=j2;
			}
		};
	};
	maxValue=sqrt(maxValue);
	return maxValue;
}

double cmatrix::minAbs(int& i1,int& i2) {
	if(!pm) {
		cerr<<endl<<"You wanted to take the smallest absolute value element of an empty matrix."<<endl;
		return NAN;
	}
	register double minValue,tmp;
	register int j1,j2;
	i1=0;i2=0; minValue=pm[0][0].re*pm[0][0].re+pm[0][0].im*pm[0][0].im;
	for(j1=0;j1<m1;j1++) {
		for(j2=0;j2<m2;j2++) {
			tmp=pm[j1][j2].re*pm[j1][j2].re+pm[j1][j2].im*pm[j1][j2].im;
			if (minValue>tmp) {
				minValue=tmp;
				i1=j1; i2=j2;
			}
		};
	};
	minValue=sqrt(minValue);
	return minValue;
}

double cmatrix::maxAbs(void) {
	int i1,i2;
	return (*this).maxAbs(i1,i2);
}

double cmatrix::minAbs(void) {
	int i1,i2;
	return (*this).minAbs(i1,i2);
}

double cmatrix::norm(void) {
	register double norm=0.0;
	register int j1,j2;
	for(j1=0;j1<m1;j1++) {
		for(j2=0;j2<m2;j2++)
			norm+=pm[j1][j2].re*pm[j1][j2].re+pm[j1][j2].im*pm[j1][j2].im;
	};
	norm=sqrt(norm);
	return norm;
}

cmatrix cmatrix::conjugate(void) {
	register int j1,j2;
	for(j1=0;j1<m1;j1++) {
		for(j2=0;j2<m2;j2++)
			pm[j1][j2].conjugate();
	};
	return *this;
}

cmatrix cmatrix::transpose(void) {
	register int i1,i2;
	if(m1!=m2) {
		dcomplex** p=pm;
		i1=m1; m1=m2; m2=i1;
		pm=new dcomplex* [m1];
		for(i1=0;i1<m1;i1++) {
			pm[i1]=new dcomplex [m2];
			for(i2=0;i2<m2;i2++)
				pm[i1][i2]=p[i2][i1];
		};
		for(i1=0;i1<m2;i1++)
			delete [] p[i1];
		delete [] p;
	}
	else {
		register dcomplex tmp;
		for(i1=0;i1<m1;i1++) {
			for(i2=0;i2<i1;i2++) {
				tmp=pm[i1][i2];
				pm[i1][i2]=pm[i2][i1];
				pm[i2][i1]=tmp;
			};
		};
	}
	return *this;
}

cmatrix cmatrix::adjoint(void) {
	register int i1,i2;
	if(m1!=m2) {
		dcomplex** p=pm;
		i1=m1; m1=m2; m2=i1;
		pm=new dcomplex* [m1];
		for(i1=0;i1<m1;i1++) {
			pm[i1]=new dcomplex [m2];
			for(i2=0;i2<m2;i2++)
				pm[i1][i2]=p[i2][i1].conjugate();
		};
		for(i1=0;i1<m2;i1++)
			delete [] p[i1];
		delete [] p;
	}
	else {
		register dcomplex tmp;
		for(i1=0;i1<m1;i1++) {
			for(i2=0;i2<i1;i2++) {
				tmp=pm[i1][i2];
				pm[i1][i2]=pm[i2][i1].conjugate();
				pm[i2][i1]=tmp.conjugate();
			};
			pm[i1][i1]=pm[i1][i1].conjugate();
		};
	}
	return *this;
}

cmatrix cmatrix::diagonal(void) {
	register int j1,j2;
	if(m1==1) {
		cmatrix result(m2,m2);
		for(j1=0;j1<m2;j1++) {
			for(j2=0;j2<m2;j2++) {
				result.pm[j1][j2]=0.0;
			};
			result.pm[j1][j1]=pm[0][j1];
		};
		return result;
	}
	else if(m2==1) {
		cmatrix result(m1,m1);
		for(j1=0;j1<m1;j1++) {
			for(j2=0;j2<m1;j2++) {
				result.pm[j1][j2]=0.0;
			};
			result.pm[j1][j1]=pm[j1][0];
		};
		return result;
	}
	else {
		int m=(m1<m2 ? m1 : m2);
		cmatrix result(m,1);
		for(j1=0;j1<m;j1++)
			result.pm[j1][0]=pm[j1][j1];
		return result;
	}
}

cmatrix cmatrix::randomize(void) {
	register double ran;
	register int i1,i2;
	
	// initialize the random number generator
	time_t now=time(0);
	tm *ltime=localtime(&now);
	int secs=1+ltime->tm_sec;
	for(i1=0;i1<secs;i1++)
		ran=ran_matrix(&idum_matrix);
	
	// fill matrix with random numbers
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++) {
			pm[i1][i2].re=ran_matrix(&idum_matrix);
			pm[i1][i2].im=ran_matrix(&idum_matrix);
		};
	};
	return *this;
}

dcomplex cmatrix::trace(void) {
	register int i;
	register dcomplex tr=0.0;
	if(m1!=m2) {
		cerr<<endl<<"You took the trace of a non-rectangular matrix!"<<endl;
	}
	for(i=0;i<(m1<m2?m1:m2);i++)
		tr+=pm[i][i];
	return tr;
}


// Operators
dcomplex* cmatrix::operator [] (int i) {
	return pm[i];
}

cmatrix cmatrix::operator = (const cmatrix& M) {
	register int i1,i2;
	if(M.m1!=m1 || M.m2!=m2) {
		for(i1=0;i1<m1;i1++) {
			delete [] pm[i1];
		};
		delete [] pm;
		m1=M.m1; m2=M.m2;
		pm=new dcomplex* [m1];
		for(i1=0;i1<m1;i1++) {
			pm[i1]=new dcomplex [m2];
		};
	}
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++) {
			pm[i1][i2]=M.pm[i1][i2];
		};
	};
	return *this;
}

cmatrix cmatrix::operator = (dmatrix& M) {
	register int i1,i2;
	int n1,n2; M.size(n1,n2);
	if(n1!=m1 || n2!=m2) {
		for(i1=0;i1<m1;i1++) {
			delete [] pm[i1];
		};
		delete [] pm;
		m1=n1; m2=n2;
		pm=new dcomplex* [m1];
		for(i1=0;i1<m1;i1++) {
			pm[i1]=new dcomplex [m2];
		};
	}
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++) {
			pm[i1][i2]=M[i1][i2];
		};
	};
	return *this;
}

cmatrix cmatrix::operator = (const double& d) {
	int i1;
	if(m1!=1 || m2!=1) {
		for(i1=0;i1<m1;i1++) {
			delete [] pm[i1];
		};
		delete [] pm;
		m1=1; m2=1;
		pm=new dcomplex* [1];
		pm[0]=new dcomplex [1];
	}
	pm[0][0]=d;
	return *this;  
}

cmatrix cmatrix::operator = (const dcomplex& c) {
	int i1;
	if(m1!=1 || m2!=1) {
		for(i1=0;i1<m1;i1++) {
			delete [] pm[i1];
		};
		delete [] pm;
		m1=1; m2=1;
		pm=new dcomplex* [1];
		pm[0]=new dcomplex [1];
	}
	pm[0][0]=c;
	return *this;  
}

cmatrix cmatrix::operator += (const cmatrix& B) {
	register int i1,i2;
	if(B.m1==1 && B.m2==1) {
		(*this)+=B.pm[0][0];
		return *this;
	}
	else if(m1==1 && m2==1) {
		dcomplex tmp=pm[0][0];
		*this=B;
		*this+=tmp;
		return *this;
	}
	if(m1!=B.m1 || m2!=B.m2) {
		cerr<<endl<<"You can't add two matrices of non-matching sizes. Nothing has been done."<<endl;
		return *this;
	}
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]+=B.pm[i1][i2];
	};
	return *this;
}

cmatrix cmatrix::operator += (dmatrix& B) {
	register int i1,i2;
	int Bm1,Bm2;
	B.size(Bm1,Bm2);
	if(Bm1==1 && Bm2==1) {
		(*this)+=B[0][0];
		return *this;
	}
	else if(m1==1 && m2==1) {
		dcomplex tmp=pm[0][0];
		*this=B;
		*this+=tmp;
		return *this;
	}
	if(m1!=Bm1 || m2!=Bm2) {
		cerr<<endl<<"You can't add two matrices of non-matching sizes. Nothing has been done."<<endl;
		return *this;
	}
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]+=B[i1][i2];
	};
	return *this;
}

cmatrix cmatrix::operator += (const dcomplex& d) {
	register int i1,i2;
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]+=d;
	};
	return *this;
}

cmatrix cmatrix::operator += (const double& d) {
	register int i1,i2;
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]+=d;
	};
	return *this;
}

cmatrix cmatrix::operator -= (const cmatrix& B) {
	register int i1,i2;
	if(B.m1==1 && B.m2==1) {
		(*this)-=B.pm[0][0];
		return *this;
	}
	else if(m1==1 && m2==1) {
		dcomplex tmp=pm[0][0];
		*this=B;
		*this-=tmp;
		return *this;
	}
	if(m1!=B.m1 || m2!=B.m2) {
		cerr<<endl<<"You can't subtract two matrices of non-matching sizes. Nothing has been done."<<endl;
		return *this;
	}
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]-=B.pm[i1][i2];
	};
	return *this;
}

cmatrix cmatrix::operator -= (dmatrix& B) {
	register int i1,i2;
	int Bm1,Bm2;
	B.size(Bm1,Bm2);
	if(Bm1==1 && Bm2==1) {
		(*this)-=B[0][0];
		return *this;
	}
	else if(m1==1 && m2==1) {
		dcomplex tmp=pm[0][0];
		*this=B;
		*this-=tmp;
		return *this;
	}
	if(m1!=Bm1 || m2!=Bm2) {
		cerr<<endl<<"You can't subtract two matrices of non-matching sizes. Nothing has been done."<<endl;
		return *this;
	}
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]-=B[i1][i2];
	};
	return *this;
}

cmatrix cmatrix::operator -= (const dcomplex& d) {
	register int i1,i2;
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]-=d;
	};
	return *this;
}

cmatrix cmatrix::operator -= (const double& d) {
	register int i1,i2;
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]-=d;
	};
	return *this;
}

cmatrix cmatrix::operator *= (const cmatrix& B) {
	register int i1,i2,i3;
	if(B.m1==1 && B.m2==1) {
		(*this)*=B.pm[0][0];
		return *this;
	}
	else if(m1==1 && m2==1) {
		dcomplex tmp=pm[0][0];
		*this=B;
		(*this)*=tmp;
		return *this;
	}
	else {
		if(m2!=B.m1) {
			cerr<<endl<<"You can't multiply two matrices of non-matching sizes. Nothing has been done."<<endl;
			return *this;
		}
		dcomplex** p;
		p=new dcomplex* [m1];
		for(i1=0;i1<m1;i1++)
			p[i1]=new dcomplex [B.m2];
		for(i1=0;i1<m1;i1++) {
			for(i2=0;i2<B.m2;i2++) {
				for(i3=0,p[i1][i2]=0.0;i3<m2;i3++)
					p[i1][i2]+=(pm[i1][i3]*B.pm[i3][i2]);
			};
		};
		for(i1=0;i1<m1;i1++)
			delete [] pm[i1];
		delete [] pm;
		pm=p;
		m2=B.m2;
		return *this;
	}
}

cmatrix cmatrix::operator *= (dmatrix& B) {
	register int i1,i2,i3;
	int Bm1,Bm2;
	B.size(Bm1,Bm2);
	if(Bm1==1 && Bm2==1) {
		(*this)*=B[0][0];
		return *this;
	}
	else if(m1==1 && m2==1) {
		dcomplex tmp=pm[0][0];
		*this=B;
		(*this)*=tmp;
		return *this;
	}
	else {
		if(m2!=Bm1) {
			cerr<<endl<<"You can't multiply two matrices of non-matching sizes. Nothing has been done."<<endl;
			return *this;
		}
		dcomplex** p;
		p=new dcomplex* [m1];
		for(i1=0;i1<m1;i1++)
			p[i1]=new dcomplex [Bm2];
		for(i1=0;i1<m1;i1++) {
			for(i2=0;i2<Bm2;i2++) {
				for(i3=0,p[i1][i2]=0.0;i3<m2;i3++)
					p[i1][i2]+=(pm[i1][i3]*B[i3][i2]);
			};
		};
		for(i1=0;i1<m1;i1++)
			delete [] pm[i1];
		delete [] pm;
		pm=p;
		m2=Bm2;
		return *this;
	}
}

cmatrix cmatrix::operator *= (const dcomplex& d) {
	register int i1,i2;
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]*=d;
	};
	return *this;
}

cmatrix cmatrix::operator *= (const double& d) {
	register int i1,i2;
	for(i1=0;i1<m1;i1++) {
		for(i2=0;i2<m2;i2++)
			pm[i1][i2]*=d;
	};
	return *this;
}

cmatrix operator + (const cmatrix& A,const cmatrix& B) {
	register int i1,i2;
	cmatrix result;
	if(A.m1==1 && A.m2==1) {
		result.m1=B.m1;
		result.m2=B.m2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=B.pm[i1][i2]+A.pm[0][0];
		};
	}
	else if(B.m1==1 && B.m2==1) {
		result.m1=A.m1;
		result.m2=A.m2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A.pm[i1][i2]+B.pm[0][0];
		};
	}
	else {
		if(A.m1!=B.m1 || A.m2!=B.m2) {
			cerr<<endl<<"You can't add two matrices of different size."<<endl;
			return result;
		}
		else {
			result.m1=A.m1;
			result.m2=A.m2;
			result.pm=new dcomplex* [result.m1];
			for(i1=0;i1<result.m1;i1++) {
				result.pm[i1]=new dcomplex [result.m2];
				for(i2=0;i2<result.m2;i2++)
					result.pm[i1][i2]=A.pm[i1][i2]+B.pm[i1][i2];
			};
		}
	}
	return result;
}

cmatrix operator + (const cmatrix& A,dmatrix& B) {
	register int i1,i2;
	int Bm1,Bm2;
	B.size(Bm1,Bm2);
	cmatrix result;
	if(A.m1==1 && A.m2==1) {
		result.m1=Bm1;
		result.m2=Bm2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=B[i1][i2]+A.pm[0][0];
		};
	}
	else if(Bm1==1 && Bm2==1) {
		result.m1=A.m1;
		result.m2=A.m2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A.pm[i1][i2]+B[0][0];
		};
	}
	else {
		if(A.m1!=Bm1 || A.m2!=Bm2) {
			cerr<<endl<<"You can't add two matrices of different size."<<endl;
			return result;
		}
		else {
			result.m1=A.m1;
			result.m2=A.m2;
			result.pm=new dcomplex* [result.m1];
			for(i1=0;i1<result.m1;i1++) {
				result.pm[i1]=new dcomplex [result.m2];
				for(i2=0;i2<result.m2;i2++)
					result.pm[i1][i2]=A.pm[i1][i2]+B[i1][i2];
			};
		}
	}
	return result;
}

cmatrix operator + (dmatrix& A,const cmatrix& B) {
	cmatrix result;
	result=B+A;
	return result;
}

cmatrix operator + (const cmatrix& A,const dcomplex& d) {
	register int i1,i2;
	cmatrix result(A.m1,A.m2);
	for(i1=0;i1<A.m1;i1++) {
		for(i2=0;i2<A.m2;i2++)
			result.pm[i1][i2]=A.pm[i1][i2]+d;
	};
	return result;
}

cmatrix operator + (const dcomplex& d,const cmatrix& A) {
	cmatrix result;
	result=A+d;
	return result;
}

cmatrix operator + (const cmatrix& A,const double& d) {
	register int i1,i2;
	cmatrix result(A.m1,A.m2);
	for(i1=0;i1<A.m1;i1++) {
		for(i2=0;i2<A.m2;i2++)
			result.pm[i1][i2]=A.pm[i1][i2]+d;
	};
	return result;
}

cmatrix operator + (const double& d,const cmatrix& A) {
	cmatrix result;
	result=A+d;
	return result;
}

cmatrix operator - (const cmatrix& A,const cmatrix& B) {
	register int i1,i2;
	cmatrix result;
	if(A.m1==1 && A.m2==1) {
		result.m1=B.m1;
		result.m2=B.m2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A.pm[0][0]-B.pm[i1][i2];
		};
	}
	else if(B.m1==1 && B.m2==1) {
		result.m1=A.m1;
		result.m2=A.m2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A.pm[i1][i2]-B.pm[0][0];
		};
	}
	else {
		if(A.m1!=B.m1 || A.m2!=B.m2) {
			cerr<<endl<<"You can't add two matrices of different size."<<endl;
			return result;
		}
		else {
			result.m1=A.m1;
			result.m2=A.m2;
			result.pm=new dcomplex* [result.m1];
			for(i1=0;i1<result.m1;i1++) {
				result.pm[i1]=new dcomplex [result.m2];
				for(i2=0;i2<result.m2;i2++)
					result.pm[i1][i2]=A.pm[i1][i2]-B.pm[i1][i2];
			};
		}
	}
	return result;
}

cmatrix operator - (const cmatrix& A,dmatrix& B) {
	register int i1,i2;
	int Bm1,Bm2;
	B.size(Bm1,Bm2);
	cmatrix result;
	if(A.m1==1 && A.m2==1) {
		result.m1=Bm1;
		result.m2=Bm2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A.pm[0][0]-B[i1][i2];
		};
	}
	else if(Bm1==1 && Bm2==1) {
		result.m1=A.m1;
		result.m2=A.m2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A.pm[i1][i2]-B[0][0];
		};
	}
	else {
		if(A.m1!=Bm1 || A.m2!=Bm2) {
			cerr<<endl<<"You can't add two matrices of different size."<<endl;
			return result;
		}
		else {
			result.m1=A.m1;
			result.m2=A.m2;
			result.pm=new dcomplex* [result.m1];
			for(i1=0;i1<result.m1;i1++) {
				result.pm[i1]=new dcomplex [result.m2];
				for(i2=0;i2<result.m2;i2++)
					result.pm[i1][i2]=A.pm[i1][i2]-B[i1][i2];
			};
		}
	}
	return result;
}

cmatrix operator - (dmatrix& A,const cmatrix& B) {
	register int i1,i2;
	int Am1,Am2;
	A.size(Am1,Am2);
	cmatrix result;
	if(Am1==1 && Am2==1) {
		result.m1=B.m1;
		result.m2=B.m2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A[0][0]-B.pm[i1][i2];
		};
	}
	else if(B.m1==1 && B.m2==1) {
		result.m1=Am1;
		result.m2=Am2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A[i1][i2]-B.pm[0][0];
		};
	}
	else {
		if(Am1!=B.m1 || Am2!=B.m2) {
			cerr<<endl<<"You can't add two matrices of different size."<<endl;
			return result;
		}
		else {
			result.m1=Am1;
			result.m2=Am2;
			result.pm=new dcomplex* [result.m1];
			for(i1=0;i1<result.m1;i1++) {
				result.pm[i1]=new dcomplex [result.m2];
				for(i2=0;i2<result.m2;i2++)
					result.pm[i1][i2]=A[i1][i2]-B.pm[i1][i2];
			};
		}
	}
	return result;
}

cmatrix operator - (const cmatrix& A,const dcomplex& d) {
	register int i1,i2;
	cmatrix result(A.m1,A.m2);
	for(i1=0;i1<A.m1;i1++) {
		for(i2=0;i2<A.m2;i2++)
			result.pm[i1][i2]=A.pm[i1][i2]-d;
	};
	return result;
}

cmatrix operator - (const dcomplex& d,const cmatrix& A) {
	register int i1,i2;
	cmatrix result(A.m1,A.m2);
	for(i1=0;i1<A.m1;i1++) {
		for(i2=0;i2<A.m2;i2++)
			result.pm[i1][i2]=d-A.pm[i1][i2];
	};
	return result;
}

cmatrix operator - (const cmatrix& A,const double& d) {
	register int i1,i2;
	cmatrix result(A.m1,A.m2);
	for(i1=0;i1<A.m1;i1++) {
		for(i2=0;i2<A.m2;i2++)
			result.pm[i1][i2]=A.pm[i1][i2]-d;
	};
	return result;
}

cmatrix operator - (const double& d,const cmatrix& A) {
	register int i1,i2;
	cmatrix result(A.m1,A.m2);
	for(i1=0;i1<A.m1;i1++) {
		for(i2=0;i2<A.m2;i2++)
			result.pm[i1][i2]=d-A.pm[i1][i2];
	};
	return result;
}

cmatrix operator * (const cmatrix& A,const cmatrix& B) {
	register int i1,i2;
	cmatrix result;
	if(A.m1==1 && A.m2==1) {
		result.m1=B.m1;
		result.m2=B.m2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=B.pm[i1][i2]*A.pm[0][0];
		};
	}
	else if(B.m1==1 && B.m2==1) {
		result.m1=A.m1;
		result.m2=A.m2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A.pm[i1][i2]*B.pm[0][0];
		};
	}
	else {
		if(A.m2!=B.m1) {
			cerr<<endl<<"You can't multiply two matrices of non-matching sizes."<<endl;
			return result;
		}
		else {
			register int i3;
			result.m1=A.m1;
			result.m2=B.m2;
			result.pm=new dcomplex* [result.m1];
			for(i1=0;i1<result.m1;i1++) {
				result.pm[i1]=new dcomplex [result.m2];
				for(i2=0;i2<result.m2;i2++) {
					for(i3=0,result.pm[i1][i2]=0.0;i3<A.m2;i3++)
						result.pm[i1][i2]+=A.pm[i1][i3]*B.pm[i3][i2];
				};
			};
		}
	}
	return result;
}

cmatrix operator * (const cmatrix& A,dmatrix& B) {
	register int i1,i2;
	cmatrix result;
	int Am1,Am2,Bm1,Bm2;
	Am1=A.m1; Am2=A.m2;
	B.size(Bm1,Bm2);
	if(A.m1==1 && A.m2==1) {
		result.m1=Bm1;
		result.m2=Bm2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=B[i1][i2]*A.pm[0][0];
		};
	}
	else if(Bm1==1 && Bm2==1) {
		result.m1=Am1;
		result.m2=Am2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A.pm[i1][i2]*B[0][0];
		};
	}
	else {
		if(Am2!=Bm1) {
			cerr<<endl<<"You can't multiply two matrices of non-matching sizes."<<endl;
			return result;
		}
		else {
			register int i3;
			result.m1=Am1;
			result.m2=Bm2;
			result.pm=new dcomplex* [result.m1];
			for(i1=0;i1<result.m1;i1++) {
				result.pm[i1]=new dcomplex [result.m2];
				for(i2=0;i2<result.m2;i2++) {
					for(i3=0,result.pm[i1][i2]=0.0;i3<Am2;i3++)
						result.pm[i1][i2]+=A.pm[i1][i3]*B[i3][i2];
				};
			};
		}
	}
	return result;
}

cmatrix operator * (dmatrix& A,const cmatrix& B) {
	register int i1,i2;
	cmatrix result;
	int Am1,Am2,Bm1,Bm2;
	A.size(Am1,Am2); Bm1=B.m1; Bm2=B.m2;
	if(Am1==1 && Am2==1) {
		result.m1=Bm1;
		result.m2=Bm2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=B.pm[i1][i2]*A[0][0];
		};
	}
	else if(Bm1==1 && Bm2==1) {
		result.m1=Am1;
		result.m2=Am2;
		result.pm=new dcomplex* [result.m1];
		for(i1=0;i1<result.m1;i1++) {
			result.pm[i1]=new dcomplex [result.m2];
			for(i2=0;i2<result.m2;i2++)
				result.pm[i1][i2]=A[i1][i2]*B.pm[0][0];
		};
	}
	else {
		if(Am2!=Bm1) {
			cerr<<endl<<"You can't multiply two matrices of non-matching sizes."<<endl;
			return result;
		}
		else {
			register int i3;
			result.m1=Am1;
			result.m2=Bm2;
			result.pm=new dcomplex* [result.m1];
			for(i1=0;i1<result.m1;i1++) {
				result.pm[i1]=new dcomplex [result.m2];
				for(i2=0;i2<result.m2;i2++) {
					for(i3=0,result.pm[i1][i2]=0.0;i3<Am2;i3++)
						result.pm[i1][i2]+=A[i1][i3]*B.pm[i3][i2];
				};
			};
		}
	}
	return result;
}

cmatrix operator * (const cmatrix& A,const double& d) {
	register int i1,i2;
	cmatrix result(A.m1,A.m2);
	for(i1=0;i1<A.m1;i1++) {
		for(i2=0;i2<A.m2;i2++) {
			result.pm[i1][i2]=A.pm[i1][i2]*d;
		};
	};
	return result;
}

cmatrix operator * (const double& d,const cmatrix& A) {
	cmatrix result;
	result=A*d;
	return result;
}

cmatrix operator * (const cmatrix& A,const dcomplex& d) {
	register int i1,i2;
	cmatrix result(A.m1,A.m2);
	for(i1=0;i1<A.m1;i1++) {
		for(i2=0;i2<A.m2;i2++) {
			result.pm[i1][i2]=A.pm[i1][i2]*d;
		};
	};
	return result;
}

cmatrix operator * (const dcomplex& d,const cmatrix& A) {
	cmatrix result;
	result=A*d;
	return result;
}

// IO
ostream& operator << (ostream& os,const cmatrix& M) {
	int array_width=12;
	int precision=5;
	int i1,i2;
	cout<<setiosflags(ios::left);
	cout<<setprecision(precision);
	cout<<fixed;
	cout<<endl;
	for(i1=0;i1<M.m1;i1++) {
		for(i2=0;i2<M.m2;i2++)
			if(fabs(M.pm[i1][i2].re)>=1e5 || fabs(M.pm[i1][i2].im)>=1e5) {cout<<scientific;}
	};
	for(i1=0;i1<M.m1;i1++) {
		for(i2=0;i2<M.m2;i2++) {
			cout<<"\t  "<<(M.pm[i1][i2].re<0 ? '-' : ' ')<<setw(array_width)<<fabs(M.pm[i1][i2].re);
			cout<<(M.pm[i1][i2].im<0 ? " - i*" : " + i*")<<setw(array_width)<<fabs(M.pm[i1][i2].im);
		};
		cout<<endl;
	};
	cout.unsetf(ios::floatfield);
	return os;
}


// LU decomposition (See Numerical Recipies in C, p. 46)
cmatrix cmatrix::LU(int *indx,dcomplex *d)
// Gives back the LU decomposition of the matix. (The matrix itself will remain untouched,) 
// indx and d are auxiliary variables used in calculation of the inverse, the determinant and solving linear equations.
{
	if(m1!=m2) {
		cerr<<endl<<"You can't make the LU decomposition of a non-square matrix. Nothing has been done."<<endl;
		return *this;
	}
	else {
		int n=m1;
		register int j1,j2;
		register dcomplex** tmpMx;
		tmpMx=new dcomplex* [n+1];
		tmpMx[0]=new dcomplex [n+1];
		for(j1=1;j1<=n;j1++) {
			tmpMx[j1]=new dcomplex [n+1];
			for(j2=1;j2<=n;j2++) {
				tmpMx[j1][j2]=pm[j1-1][j2-1];
			};
		};
		
		ludcmp_dcomplex(tmpMx,n,indx,d);
		
		cmatrix result(n,n);
		delete [] tmpMx[0];
		for(j1=1;j1<=n;j1++) {
			for(j2=1;j2<=n;j2++) {
				result.pm[j1-1][j2-1]=tmpMx[j1][j2];
			};
			delete [] tmpMx[j1];
		};
		delete [] tmpMx;
		return result;
	}
}

cmatrix cmatrix::LU(void) 
// Same as the previous one. It just doesn't keep track of indx and d.
{
	int* indx;
	indx=new int [m1];
	dcomplex d;
	cmatrix result;
	result=(*this).LU(indx,&d);
	delete [] indx;
	return result;
}

cmatrix cmatrix::LU_backsub(int *indx, cmatrix b)
{
	if(m1!=m2 || m2!=b.m1 || b.m2!=1) {
		cerr<<endl<<"Inappropriate size of the matrix or the vector in function LU_backsub.";
		cerr<<"Nothing has been done."<<endl;
		cmatrix emptyMx;
		return emptyMx;
	}
	else {
		int n=m1;
		register dcomplex* b_;
		b_=new dcomplex [m1+1];
		register int j1,j2;
		for(j1=1;j1<=n;j1++)
			b_[j1]=b.pm[j1-1][0];
		register dcomplex** tmpMx;
		tmpMx=new dcomplex* [n+1];
		tmpMx[0]=new dcomplex [n+1];
		for(j1=1;j1<=n;j1++) {
			tmpMx[j1]=new dcomplex [n+1];
			for(j2=1;j2<=n;j2++) {
				tmpMx[j1][j2]=pm[j1-1][j2-1];
			};
		};
		
		lubksb_dcomplex(tmpMx,n,indx,b_);
		
		cmatrix result(b.m1,b.m2);
		for(j1=1;j1<=n;j1++) {
			result.pm[j1-1][0]=b_[j1];
			delete [] tmpMx[j1];
		};
		delete [] tmpMx[0];
		delete [] tmpMx;
		delete [] b_;
		return result;
	}
}

// Linear equation solver
cmatrix linearSolve(cmatrix& A,cmatrix& b) {
	// Solves the linear equation A*x=b;
	if(A.m1!=A.m2 || A.m2!=b.m1 || b.m2!=1) {
		cerr<<endl<<"Inappropriate size of the matrix or the vector in function linearSolve. ";
		cerr<<"Nothing has been done."<<endl;
		cmatrix x;
		return x;
	}
	else {
		cmatrix x(b.m1,b.m2);
		cmatrix tmpA(A.m1,A.m2);
		int* indx;
		indx=new int [A.m1];
		dcomplex d;
		tmpA=A.LU(indx,&d);
		x=tmpA.LU_backsub(indx,b);
		delete [] indx;
		return x;
	}
}

cmatrix cmatrix::invert(void) {
	if(m1!=m2) {
		cerr<<endl<<"You can't take the inverse of a non-square matrix. Nothing has been done."<<endl;
		return *this;
	}
	else {
		int n=m1;
		register int j1,j2;
		register dcomplex** tmpMx;
		tmpMx=new dcomplex* [n+1];
		tmpMx[0]=new dcomplex [n+1];
		for(j1=1;j1<=n;j1++) {
			tmpMx[j1]=new dcomplex [n+1];
			for(j2=1;j2<=n;j2++)
				tmpMx[j1][j2]=pm[j1-1][j2-1];
		};
		int* indx; indx=new int [n+1];
		dcomplex* col; col=new dcomplex [n+1];
		dcomplex d;
		
		ludcmp_dcomplex(tmpMx,n,indx,&d);
		for(j1=1;j1<=n;j1++) {
			for(j2=1;j2<=n;j2++) 
				col[j2]=0.0;
			col[j1]=1.0;
			lubksb_dcomplex(tmpMx,n,indx,col);
			for(j2=1;j2<=n;j2++) 
				pm[j2-1][j1-1]=col[j2];
		};
		delete [] indx;
		delete [] col;
		for(j1=0;j1<=n;j1++)
			delete [] tmpMx[j1];
		delete [] tmpMx;
		return *this;
	}
}

dcomplex cmatrix::determinant(void) {
	if(m1!=m2) {
		cerr<<endl<<"You can't take the determinant of a non-square matrix. Nothing has been done."<<endl;
		dcomplex det;
		det.re=NAN; det.im=NAN;
		return det;
	}
	else {
		register dcomplex det;
		register int* indx;
		indx=new int [m1];
		cmatrix tmpMx(m1,m2);
		register int j;
		
		tmpMx=(*this).LU(indx,&det);
		for(j=0;j<m1;j++)
			det*=tmpMx[j][j];
		delete [] indx;
		return det;
	}
}

// External Lapack routine for diagonalization of a Hermitian matrix
extern "C" void zheev_( char* jobz, char* uplo, int* n, dcomplex* a, int* lda,
			double* w, dcomplex* work, int* lwork, double* rwork, int* info );
void cmatrix::diagonalize_hermitian(cmatrix& eigs,cmatrix& U) {
	register int i1,i2;
	if(m1!=m2) {
		cerr<<endl<<"You can't diagonalize a non-square matrix."<<endl;
	}
	else {
		// Reset the size of eigs and U if inappropriate
		if(eigs.m1!=m1 || eigs.m2!=1) {
			for(i1=0;i1<eigs.m1;i1++)
				delete [] eigs.pm[i1];
			delete [] eigs.pm;
			eigs.m1=m1;
			eigs.m2=1;
			eigs.pm=new dcomplex* [m1];
			for(i1=0;i1<eigs.m1;i1++)
				eigs.pm[i1]=new dcomplex [1];
		}
		if(U.m1!=m1 || U.m2!=m1) {
			for(i1=0;i1<U.m1;i1++)
				delete [] U.pm[i1];
			delete [] U.pm;
			U.m1=m1;
			U.m2=m1;
			U.pm=new dcomplex* [m1];
			for(i1=0;i1<U.m1;i1++)
				U.pm[i1]=new dcomplex [m1];
		}
		// Calculate diagonalization
		register dcomplex* A;
		A=new (nothrow) dcomplex[m1*m1];
		register double* eigenvalues;
		eigenvalues=new (nothrow) double [m1];
		int n=m1; 
		int lda=m1; 			/* The leading of the array A. lda >= max(1, N).*/
		int info;			/* info == 0:  successful exit
					        < 0:  if INFO == -i, the i-th argument had an illegal value
					        > 0:  if INFO == i, the algorithm failed to converge; i off-diagonal elements of an intermediate
								   tridiagonal form did not converge to zero.*/
		int lwork;			/* Length of the array work. */
		dcomplex wkopt;			/* The variable used to get the size of the optimal workspace array. */
		dcomplex* work;			/* Workspace */
		double rwork[3*m1-2];   	/* rwork dimension should be at least max(1, 3*n-2). */
		
		for(i1=0;i1<m1;i1++) {
			for(i2=i1;i2<m1;i2++)
				A[i1*m1+i2]=pm[i2][i1];
			for(i2=0;i2<i1;i2++)
				A[i1*m1+i2]=0.0;
		};
		// Calculate the optimal workspace to be allocated.
		lwork = -1;
		zheev_("Vectors","Lower",&n,A,&lda,eigenvalues,&wkopt,&lwork,rwork,&info);
		lwork=(int)wkopt.re;
		work=new (nothrow) dcomplex[lwork];
		// Diagonalize
		zheev_("Vectors","Lower",&n,A,&lda,eigenvalues,work,&lwork,rwork,&info);
		if(info>0) {
			cerr<<endl<<"zheev_ error."<<endl;
			eigs=NAN; U=NAN;
		}
		else {
			for(i1=0;i1<m1;i1++) {
				for(i2=0;i2<m1;i2++) {
					U.pm[i2][i1]=A[i1*m1+i2];
				};
				eigs.pm[i1][0]=eigenvalues[i1];
			};
		}
		// Free workspace
		delete[] A;
		delete[] eigenvalues;
		delete[] work;
	}
}



// AUXILIARY ROUTINES

// random number generator
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
// Produces random number between 0 and 1 using multiplicative method 
// (see Numerical Recipies p. 280). Cycle around 1e8.
// Recommended instead of the built in random number generators of 
// the C language implementations.
// idum: shall be initialized as a negative long variable once in the 
// program and not touched afterwards. (This is the seed of ran1.)
double ran_matrix(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


// Numerical recipies routine for LU decomposition
#define NRANSI
#define TINY 1.0e-20;
void ludcmp(double **a, int n, int *indx, double *d)
// Given a matrix a[0..n-1][0..n-1], this routine replaces it by the LU decomposition of a rowwise
// permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
// indx[1..n] is an output vector that records the row permutation effected by the partial
// pivoting; d is output as ±1 depending on whether the number of row interchanges was even
// or odd, respectively. This routine is used in combination with lubksb to solve linear equations
// or invert a matrix.
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=new double [n+1];
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) {cerr<<endl<<"Singular matrix in routine ludcmp"<<endl;}
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	delete [] vv;
}
#undef TINY
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software */

#define NRANSI
#define TINY 1.0e-20;
void ludcmp_dcomplex(dcomplex **a, int n, int *indx, dcomplex *d)
// dcomplex version of ludcmp
{
	int i,imax,j,k;
	double temp,big,dum;
	dcomplex sum,dum_complex;
	dcomplex *vv;

	vv=new dcomplex [n+1];
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=a[i][j].abs()) > big) big=temp;
		if (big == 0.0) {cerr<<endl<<"Singular matrix in routine ludcmp"<<endl;}
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=(vv[i]*sum).abs()) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum_complex=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum_complex;
			}
			(*d).re = -((*d).re);
			(*d).im = -((*d).im);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j].abs() == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum_complex=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum_complex;
		}
	}
	delete [] vv;
}
#undef TINY
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software */


void lubksb(double **a, int n, int *indx, double* b)
// Solves the set of n linear equations A·X = B. Here a[0..n-1][0..n-1] is input, not as the matrix
// A but rather as its LU decomposition, determined by the routine ludcmp. indx[0..n-1] is input
// as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
// B, and returns with the solution vector X. a, n, and indx are not modified by this routine
// and can be left in place for successive calls with different right-hand sides b. This routine takes
// into account the possibility that b will begin with many zero elements, so it is efficient for use
// in matrix inversion.
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software */

void lubksb_dcomplex(dcomplex **a, int n, int *indx, dcomplex* b)
{
	int i,ii=0,ip,j;
	dcomplex sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum.abs()) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

