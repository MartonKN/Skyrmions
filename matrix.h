#ifndef __MATRIXH
#define __MATRIXH

#include <cmath>
#include<iostream>
#include<string>
using namespace std;


// REAL MATRICES
class dmatrix {
  public:
    // Constructors and destructors
    dmatrix(void);
    dmatrix(int n1,int n2);
    dmatrix(dmatrix&);
    dmatrix(double&);
    dmatrix(int n1,int n2,double**);
    ~dmatrix(void) {
      register int i;
      for(i=0;i<m1;i++)
	delete[] pm[i];
      delete[] pm;
    }
    
    // Functions
    void size(int& n1,int& n2);
    // Returns the size of the matrix: n1 rows, n2 coloumns.
    double max(int& i1,int& i2);
    double min(int& i1,int& i2);
    double maxAbs(int& i1,int& i2);
    double minAbs(int& i1,int& i2);
    double max(void);
    double min(void);
    double maxAbs(void);
    double minAbs(void);

    // Return the maximal/minimal element of the matrix 
    // or that of the greatest absolute value/the smallest absolute value.
    double norm(void);
    // Returns the 2-norm of the matrix.
    dmatrix transpose(void);
    // Substitutes the matrix with its transposed.
    dmatrix diagonal(void);
    // Transforms a vector to a diagonal matrix OR takes the diagonal of a matrix.
    // This new matrix is returned, and the original one is not modified.
    dmatrix randomize(void);
    // Substitutes all the elements of the matrix with random numbers.
    double trace(void);
    // Returns the trace of the matrix.
    double determinant(void);
    // Returns the determinant of the matrix.
    dmatrix invert(void);
    // Inverts the matrix.
    dmatrix LU(int *indx,double *d);
    // Gives back the LU decomposition of the matix. (The matrix itself will remain untouched,) 
    // indx and d are auxiliary variables used in calculation of the inverse, the determinant and solving linear equations.
    dmatrix LU(void);
    // Same as the previous one. It just doesn't keep track of indx and d.
    dmatrix LU_backsub(int *indx,dmatrix b);
    // LU backsubstitution (see Numerical Recipies in C, p. 46). Auxiliary function, used
    // in functions 'invert' and 'linearSolve'.
    friend dmatrix linearSolve(dmatrix& A,dmatrix& b);
    // Solves the linear equation A*x=b and returns x. A and b aren't modified.
    
    
    // Operators
    double* operator [] (int i);
    dmatrix operator = (const dmatrix&);
    dmatrix operator = (const double&);
    dmatrix operator += (const dmatrix&);
    dmatrix operator += (const double&);
    dmatrix operator -= (const dmatrix&);
    dmatrix operator -= (const double&);
    dmatrix operator *= (const dmatrix&);
    dmatrix operator *= (const double&);
    
    friend dmatrix operator + (const dmatrix&,const dmatrix&);
    friend dmatrix operator + (const double&,const dmatrix&);
    friend dmatrix operator + (const dmatrix&,const double&);
    friend dmatrix operator - (const dmatrix&,const dmatrix&);
    friend dmatrix operator - (const double&,const dmatrix&);
    friend dmatrix operator - (const dmatrix&,const double&);
    friend dmatrix operator * (const dmatrix&,const dmatrix&);
    friend dmatrix operator * (const double&,const dmatrix&);
    friend dmatrix operator * (const dmatrix&,const double&);
    
    friend ostream& operator<<(ostream&,const dmatrix&);
  private:
    double** pm; // pointer to the matrix
    int m1,m2;     // number of rows (m1) and coloumns (m2)
};


// COMPLEX MATRICES
class cmatrix {
  public:
    // Constructors and destructors
    cmatrix(void);
    cmatrix(int n1,int n2);
    cmatrix(cmatrix&);
    cmatrix(dmatrix&);
    cmatrix(dcomplex&);
    cmatrix(double&);
    cmatrix(int n1,int n2,dcomplex**);
    cmatrix(int n1,int n2,double**);
    ~cmatrix(void) {
      register int i;
      for(i=0;i<m1;i++)
	delete[] pm[i];
      delete[] pm;
    }
    
    // Functions
    void size(int& n1,int& n2);
    // Returns the size of the matrix: n1 rows, n2 coloumns.
    cmatrix re(void);
    // Substitutes the matrix with its real part.
    cmatrix im(void);
    // Substitutes the matrix with its imaginary part.
    double maxAbs(int& i1,int& i2);
    double maxAbs(void);
    // Returns the maximal absolute value element of the matrix.
    double minAbs(int& i1,int& i2);
    double minAbs(void);
    // Returns the minimal absolute value element of the matrix.
    double norm(void);
    // Returns the 2-norm of the matrix.
    cmatrix conjugate(void);
    // Substitutes the matrix with its conjugate.
    cmatrix transpose(void);
    // Substitutes the matrix with its transposed.
    cmatrix adjoint(void);
    // Substitutes the matrix with its adjoint.
    cmatrix diagonal(void);
    // Transforms a vector to a diagonal matrix OR takes the diagonal of a matrix.
    // This new matrix is returned, and the original one is not modified.
    cmatrix randomize(void);
    // Substitutes all the elements of the matrix with random numbers.
    dcomplex trace(void);
    // Returns the trace of the matrix.
    dcomplex determinant(void);
    // Returns the determinant of the matrix.
    cmatrix invert(void);
    // Inverts the matrix.
/*    void diagonalize(cmatrix& eigs,cmatrix& U);
    // Diagonalizes the matrix. M=(U*diag(eigs))*U.adjoint().*/
    void diagonalize_hermitian(cmatrix& eigs,cmatrix& U);
    // Diagonalizes the matrix. M=(U*diag(eigs))*U.adjoint().
    
    cmatrix LU(int *indx,dcomplex *d);
    // Substitutes the matrix with its LU decomposition. indx and d are auxiliary variables
    // used in calculation of the inverse, the determinant and solving linear equations.
    cmatrix LU(void);
    // Same as the previous one. It just doesn't keep track of indx and d.
    cmatrix LU_backsub(int *indx, cmatrix b);
    // LU backsubstitution (see Numerical Recipies in C, p. 46). Auxiliary function, used
    // in functions 'invert' and 'linearSolve'.
    friend cmatrix linearSolve(cmatrix& A,cmatrix& b);
    // Solves the linear equation A*x=b and returns x. A and b aren't modified.
    
    // Operators
    dcomplex* operator [] (int i);
    cmatrix operator = (const cmatrix&);
    cmatrix operator = (dmatrix&);
    cmatrix operator = (const dcomplex&);
    cmatrix operator = (const double&);
    cmatrix operator += (const cmatrix&);
    cmatrix operator += (dmatrix&);
    cmatrix operator += (const dcomplex&);
    cmatrix operator += (const double&);
    cmatrix operator -= (const cmatrix&);
    cmatrix operator -= (dmatrix&);
    cmatrix operator -= (const dcomplex&);
    cmatrix operator -= (const double&);
    cmatrix operator *= (const cmatrix&);
    cmatrix operator *= (dmatrix&);
    cmatrix operator *= (const dcomplex&);
    cmatrix operator *= (const double&);
    
    friend cmatrix operator + (const cmatrix&,const cmatrix&);
    friend cmatrix operator + (dmatrix&,const cmatrix&);
    friend cmatrix operator + (const cmatrix&,dmatrix&);
    friend cmatrix operator + (const dcomplex&,const cmatrix& );
    friend cmatrix operator + (const cmatrix&, const dcomplex&);
    friend cmatrix operator + (const double&,const cmatrix&);
    friend cmatrix operator + (const cmatrix&,const double& );
    friend cmatrix operator - (const cmatrix&,const cmatrix&);
    friend cmatrix operator - (dmatrix&,const cmatrix&);
    friend cmatrix operator - (const cmatrix&,dmatrix&);
    friend cmatrix operator - (const dcomplex&,const cmatrix& );
    friend cmatrix operator - (const cmatrix&,const dcomplex&);
    friend cmatrix operator - (const double& ,const cmatrix&);
    friend cmatrix operator - (const cmatrix&,const double& );
    friend cmatrix operator * (const cmatrix&,const cmatrix&);
    friend cmatrix operator * (dmatrix&,const cmatrix&);
    friend cmatrix operator * (const cmatrix&,dmatrix&);
    friend cmatrix operator * (const dcomplex&,const cmatrix&);
    friend cmatrix operator * (const cmatrix&,const dcomplex&);
    friend cmatrix operator * (const double&,const cmatrix&);
    friend cmatrix operator * (const cmatrix&,const double&);
    
    friend ostream& operator<<(ostream&,const cmatrix&);
  private:
    dcomplex** pm; // pointer to the matrix
    int m1,m2;     // number of rows (m1) and coloumns (m2)
};

double ran_matrix(long *idum);
void ludcmp(double **a, int n, int *indx, double *d);
void ludcmp_dcomplex(dcomplex **a, int n, int *indx, dcomplex *d);
void lubksb(double **a, int n, int *indx, double* b);
void lubksb_dcomplex(dcomplex **a, int n, int *indx, dcomplex* b);
#endif