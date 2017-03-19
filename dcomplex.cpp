#include "dcomplex.h"
using namespace std;

// Constructors and copy functions
dcomplex::dcomplex(void) {
  re=0.0;
  im=0.0;
}
dcomplex::dcomplex(double r) {
  re=r;
  im=0.0;
}
dcomplex::dcomplex(double r,double i) {
  re=r;
  im=i;
}
dcomplex::dcomplex(const dcomplex& c) {
  re=c.re;
  im=c.im;
}
dcomplex::dcomplex(fftw_complex fc) {
  re=fc[0];
  im=fc[1];
}

// = operator
dcomplex dcomplex::operator= (const dcomplex& c) {
  re=c.re;
  im=c.im;
  return *this;
}
dcomplex dcomplex::operator= (const fftw_complex& fc) {
  re=fc[0];
  im=fc[1];
  return *this;
}
dcomplex dcomplex::operator= (const double& d) {
  re=d;
  im=0.0;
  return *this;
}

// += operator
dcomplex dcomplex::operator+= (const dcomplex& c) {
  re+=c.re;
  im+=c.im;
  return *this;
}
dcomplex dcomplex::operator+= (const fftw_complex& fc) {
  re+=fc[0];
  im+=fc[1];
  return *this;
}
dcomplex dcomplex::operator+= (const double& d) {
  re+=d;
  return *this;
}

// -= operator
dcomplex dcomplex::operator-= (const dcomplex& c) {
  re-=c.re;
  im-=c.im;
  return *this;
}
dcomplex dcomplex::operator-= (const fftw_complex& fc) {
  re-=fc[0];
  im-=fc[1];
  return *this;
}
dcomplex dcomplex::operator-= (const double& d) {
  re-=d;
  return *this;
}

// + operator
dcomplex operator + (const dcomplex& c1, const dcomplex& c2) {
  dcomplex c(c1);
  c+=c2;
  return c;
}
dcomplex operator+ (const dcomplex& c1, const fftw_complex& fc2) {
  dcomplex c(c1);
  c+=fc2;
  return c;
}
dcomplex operator+ (const dcomplex& c1, const double& d2) {
  dcomplex c(c1);
  c.re+=d2;
  return c;
}
dcomplex operator+ (const fftw_complex& fc1, const dcomplex& c2) {
  dcomplex c(c2);
  c+=fc1;
  return c;
}
dcomplex operator+ (const double& d1, const dcomplex& c2) {
  dcomplex c(c2);
  c.re+=d1;
  return c;
}

// - operator
dcomplex operator - (const dcomplex& c1, const dcomplex& c2) {
  dcomplex c(c1);
  c-=c2;
  return c;
}
dcomplex operator - (const dcomplex& c1, const fftw_complex& fc2) {
  dcomplex c(c1);
  c-=fc2;
  return c;
}
dcomplex operator - (const dcomplex& c1, const double& d2) {
  dcomplex c(c1);
  c.re-=d2;
  return c;
}
dcomplex operator - (const fftw_complex& fc1, const dcomplex& c2) {
  dcomplex c;
  c.re=fc1[0]-c2.re;
  c.im=fc1[1]-c2.im;
  return c;
}
dcomplex operator - (const double& d1, const dcomplex& c2) {
  dcomplex c;
  c.re=d1-c2.re;
  c.im=-c2.im;
  return c;
}

// *= operator
dcomplex dcomplex::operator*= (const dcomplex& c) {
  double tmp;
  tmp=re*c.re-im*c.im;
  im=im*c.re+re*c.im;
  re=tmp;
  return *this;
}
dcomplex dcomplex::operator*= (const fftw_complex& fc) {
  double tmp;
  tmp=re*fc[0]-im*fc[1];
  im=im*fc[0]+re*fc[1];
  re=tmp;
  return *this;
}
dcomplex dcomplex::operator*= (const double& d) {
  re*=d;
  im*=d;
  return *this;
}

// * operator
dcomplex operator* (const dcomplex& c1, const dcomplex& c2) {
  dcomplex c(c1);
  c*=c2;
  return c;
}
dcomplex operator* (const dcomplex& c1, const fftw_complex& fc2) {
  dcomplex c(c1);
  c*=fc2;
  return c;
}
dcomplex operator* (const dcomplex& c1, const double& d2) {
  dcomplex c(c1);
  c.re*=d2;
  c.im*=d2;
  return c;
}
dcomplex operator* (const fftw_complex& fc1, const dcomplex& c2) {
  dcomplex c(c2);
  c*=fc1;
  return c;
}
dcomplex operator* (const double& d1, const dcomplex& c2) {
  dcomplex c(c2);
  c.re*=d1;
  c.im*=d1;
  return c;
}

// / operator
dcomplex operator / (const dcomplex& c1, const dcomplex& c2) {
  dcomplex c;
  double r,den;
  if (fabs(c2.re) >= fabs(c2.im)) {
    r=c2.im/c2.re;
    den=c2.re+r*c2.im;
    c.re=(c1.re+r*c1.im)/den;
    c.im=(c1.im-r*c1.re)/den;
  } else {
    r=c2.re/c2.im;
    den=c2.im+r*c2.re;
    c.re=(c1.re*r+c1.im)/den;
    c.im=(c1.im*r-c1.re)/den;
  }
  return c;
}
dcomplex operator / (const fftw_complex& fc1, const dcomplex& c2) {
  dcomplex c1; c1.re=fc1[0]; c1.im=fc1[1];
  return c1/c2;
}
dcomplex operator / (const dcomplex& c1, const fftw_complex& fc2){
  dcomplex c2; c2.re=fc2[0]; c2.im=fc2[1];
  return c1/c2;  
}
dcomplex operator / (const dcomplex& c1, const double& d){
  double tmp=1.0/d;
  dcomplex c;
  c.re=c1.re*tmp;
  c.im=c1.im*tmp;
  return c;
}
dcomplex operator / (const double& d, const dcomplex& c2){
  dcomplex c1; c1.re=d; c1.im=0.0;
  return c1/c2;
}


// IO
istream& operator >> (istream& is,dcomplex& c) {
  char str[200];
  is.getline(str,200);
  if(sscanf(str,"%lf+i*%lf",&c.re,&c.im)==2) {}
  else if(sscanf(str,"%lf-i*%lf",&c.re,&c.im)==2) {
    c.im*=-1.0;
  }  
  else {
    c=dcomplex(0.0);
  }
  return is;
}

ostream& operator << (ostream& os,const dcomplex& c) {
  os<<" "<<c.re<<" "<<(c.im>=0.0 ? '+' : '-')<<" i*"<<abs(c.im)<<" ";
  return os;
}