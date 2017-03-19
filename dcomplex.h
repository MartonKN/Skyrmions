#ifndef __DCOMPLEX
#define __DCOMPLEX

#include <cmath>
#include<iostream>
#include<string>
#include<fftw3.h>
using namespace std;

class dcomplex {
  public:
    double re,im;
    // Constructors
    dcomplex(void);
    dcomplex(double r);
    dcomplex(double r,double i);
    dcomplex(const dcomplex& c);
    explicit dcomplex(fftw_complex fc);
    // Overdefined operators
    dcomplex operator = (const dcomplex& c);
    dcomplex operator = (const fftw_complex& fc);
    dcomplex operator = (const double& d);

    dcomplex operator += (const dcomplex& c);
    dcomplex operator += (const fftw_complex& fc);
    dcomplex operator += (const double& d);

    dcomplex operator -= (const dcomplex& c);
    dcomplex operator -= (const fftw_complex& fc);
    dcomplex operator -= (const double& d);

    dcomplex operator *= (const dcomplex& c);
    dcomplex operator *= (const fftw_complex& fc);
    dcomplex operator *= (const double& d);

    friend dcomplex operator + (const dcomplex& c1, const dcomplex& c2);
    friend dcomplex operator + (const dcomplex& c1, const fftw_complex& fc2);
    friend dcomplex operator + (const dcomplex& c1, const double& d2);
    friend dcomplex operator + (const fftw_complex& fc1, const dcomplex& c2);
    friend dcomplex operator + (const double& d1, const dcomplex& c2);

    friend dcomplex operator - (const dcomplex& c1, const dcomplex& c2);
    friend dcomplex operator - (const dcomplex& c1, const fftw_complex& fc2);
    friend dcomplex operator - (const dcomplex& c1, const double& d2);
    friend dcomplex operator - (const fftw_complex& fc1, const dcomplex& c2);
    friend dcomplex operator - (const double& d1, const dcomplex& c2);

    friend dcomplex operator * (const dcomplex& c1, const dcomplex& c2);
    friend dcomplex operator * (const dcomplex& c1, const fftw_complex& fc2);
    friend dcomplex operator * (const dcomplex& c1, const double& d2);
    friend dcomplex operator * (const fftw_complex& fc1, const dcomplex& c2);
    friend dcomplex operator * (const double& d1, const dcomplex& c2);
    
    friend dcomplex operator / (const dcomplex& c1, const dcomplex& c2);
    friend dcomplex operator / (const fftw_complex& fc1, const dcomplex& c2);
    friend dcomplex operator / (const dcomplex& c1, const fftw_complex& fc2);
    friend dcomplex operator / (const dcomplex& c1, const double& d);
    friend dcomplex operator / (const double& d, const dcomplex& c2);
    
    friend istream& operator>>(istream&,dcomplex& c);
    friend ostream& operator<<(ostream&,const dcomplex& c);

    
    // Functions
    double abs(void) {
      double x,y,ans,temp;
      x=fabs(re);
      y=fabs(im);
      if (x == 0.0)
	ans=y;
      else if (y == 0.0)
	ans=x;
      else if (x > y) {
	temp=y/x;
	ans=x*sqrt(1.0+temp*temp);
      } else {
	temp=x/y;
	ans=y*sqrt(1.0+temp*temp);
      }
      return ans;
    };
    
    double sqr_abs(void) {return (re*re+im*im);};
    
    dcomplex conjugate(void) {im*=-1.0; return *this;}
    
    dcomplex csqrt(void) {
      dcomplex c;
      float x,y,w,r;
      if ((re == 0.0) && (im == 0.0)) {
	c.re=0.0;
	c.im=0.0;
	return c;
      } else {
	x=fabs(re);
	y=fabs(im);
	if (x >= y) {
	  r=y/x;
	  w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
	} else {
	  r=x/y;
	  w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
	}
	if (re >= 0.0) {
	  c.re=w;
	  c.im=im/(2.0*w);
	} else {
	  c.im=(im >= 0) ? w : -w;
	  c.re=im/(2.0*c.im);
	}
	return c;
      }
    }

};
#endif