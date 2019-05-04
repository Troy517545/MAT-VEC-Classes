
// matrix class
#ifndef MAT_H
#define MAT_H

#include "VEC.h"

class MAT {
public:
    MAT(int dim);            // uninit constructor
    MAT(const MAT &m1);      // copy constructor
    MAT(int dim, double *v); // init constructor
    ~MAT();                  // destructor
    int dim();               // return dimension of the matrix
    MAT tpose();             // transpose
    double det();            // get determinant
    
    MAT &operator-();          // unary operator, negative value
    MAT &operator=(MAT m1);    // assignment
    MAT &operator+=(MAT &m1);  // m+=m1;
    MAT &operator-=(MAT &m1);  // m-=m1;
    MAT &operator*=(double a); // m*=dbl;
    MAT &operator/=(double a); // m/=dbl;
    MAT operator+(MAT m1);     // m1+m2
    MAT operator-(MAT m1);     // m1-m2
    MAT operator*(MAT m1);     // m1*m2
    VEC &operator[](int M);    // m'th row
    VEC operator*(VEC v1);     // m x v1
    MAT operator*(double a);   // m*dbl
    MAT operator/(double a);   // m/dbl
    
    friend MAT operator*(double a, MAT &m1); // dbl x m
    friend VEC operator*(VEC &v1, MAT &m1);  // vT x m
    void print();                            //
    int get_size();
    
    bool operator==(MAT);
    
private:
    int n;    // dedine nxn matrix
    VEC **va; // array of n pointer to vectors
};

double determinant(MAT m, int n);

MAT operator*(double a, const MAT &m1); // dbl x m
VEC operator*(VEC &v1, MAT &m1);        // vT x m

MAT &luFact(MAT &m1);

MAT &cholesky(MAT &A);
VEC fwdSubs(MAT &m1, VEC b);
VEC bckSubs(MAT &m1, VEC b);
VEC choSolve(MAT &L, VEC b);

int jacobi(MAT &A, VEC b, VEC &x, int maxIter, double tol);
int gaussSeidel(MAT &A, VEC b, VEC &x, int maxIter, double tol);
int sgs(MAT &A, VEC b, VEC &x, int maxIter, double tol);

int cg(MAT &A, VEC b, VEC &x, int maxIter, double tol);

int EVpwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter);
int EViPwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter);
int EViPwrShft(MAT &A, VEC &q0, double &lambda, double mu, double tol,
               int maxiter);

void qrFact(MAT A, MAT &Q, MAT &R);
int EVqr(MAT &A, double tol, int maxiter);
int EVqrShifted(MAT &A, double mu, double tol, int maxiter);

double Lagrange(double x, VEC &XDATA, VEC &YDATA);
double NEV(double x, VEC &XDATA, VEC &YDATA);

void splineM(int N, VEC &X, VEC &Y, VEC &M); // generate spline momentum M
double spline(double x, int N, VEC &X, VEC &Y, VEC &M); // spline interp at x
void parametric_spline(int subintervalNumber, int N, VEC &X, VEC &Y,
                       VEC &spline_x, VEC &spline_y); // p[aratric spline interp

#endif
