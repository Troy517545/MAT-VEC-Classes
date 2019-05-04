
// MAT class functions

#include "MAT.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

MAT::MAT(int dim) {
    n = dim;
    va = (VEC **)malloc(n * sizeof(VEC *));
    for (int i = 0; i < n; i++) {
        va[i] = newVEC(n);
    }
}

MAT::MAT(const MAT &m1) {
    VEC **vsrc = m1.va;
    n = m1.n;
    va = (VEC **)malloc(n * sizeof(VEC *));
    for (int i = 0; i < n; i++) {
        va[i] = newVEC(n);
        *va[i] = *vsrc[i];
    }
}

MAT::MAT(int dim, double *v) {
    n = dim;
    va = (VEC **)malloc(n * sizeof(VEC *));
    for (int i = 0; i < n; i++) {
        va[i] = newVEC(n);
        for (int j = 0; j < n; j++) {
            (*va[i])[j] = *(v++);
        }
    }
}

MAT::~MAT() {
    for (int i = n - 1; i >= 0; i--) {
        (*va[i]).~VEC();
    }
    free(va);
}

int MAT::dim() { return n; }

MAT MAT::tpose() {
    MAT mnew(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mnew[i][j] = (*va[j])[i];
        }
    }
    return mnew;
}

MAT &MAT::operator-() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            (*va[i])[j] = -(*va[i])[j];
        }
    }
    return *this;
}

MAT &MAT::operator=(MAT m1) {
    for (int i = 0; i < n; i++)
        (*va[i]) = m1[i];
    return *this;
}

MAT &MAT::operator+=(MAT &m1) {
    for (int i = 0; i < n; i++)
        (*va[i]) += m1[i];
    return *this;
}

MAT &MAT::operator-=(MAT &m1) {
    for (int i = 0; i < n; i++)
        (*va[i]) -= m1[i];
    return *this;
}

MAT &MAT::operator*=(double a) {
    for (int i = 0; i < n; i++)
        (*va[i]) *= a;
    return *this;
}

MAT &MAT::operator/=(double a) {
    for (int i = 0; i < n; i++)
        (*va[i]) /= a;
    return *this;
}

MAT MAT::operator+(MAT m1) {
    MAT s(n);
    for (int i = 0; i < n; i++)
        s[i] = (*va[i]) + m1[i];
    return s;
}

MAT MAT::operator-(MAT m1) {
    MAT s(n);
    for (int i = 0; i < n; i++)
        s[i] = (*va[i]) - m1[i];
    return s;
}

MAT MAT::operator*(MAT m1) {
    MAT z(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            z[i][j] = 0;
            for (int k = 0; k < n; k++)
                z[i][j] += ((*va[i])[k] * m1[k][j]);
        }
    }
    return z;
}

VEC &MAT::operator[](int m) { return (*va[m]); }

VEC MAT::operator*(VEC v1) {
    VEC s(n);
    for (int i = 0; i < n; i++) {
        s[i] = (*va[i]) * v1;
    }
    return s;
}

MAT MAT::operator*(double a) {
    MAT s(*this);
    for (int i = 0; i < n; i++) {
        (*va[i]) *= a;
    }
    return s;
}

MAT MAT::operator/(double a) {
    MAT s(*this);
    for (int i = 0; i < n; i++) {
        (*va[i]) /= a;
    }
    return s;
}

void MAT::print() {
    for (int i = 0; i < n; i++) {
        (*va[i]).print();
    }
    std::cout << std::endl;
}

VEC operator*(VEC &v1, MAT &m1) {
    VEC v2(m1.n);
    for (int i = 0; i < m1.n; i++) {
        v2[i] = 0;
        for (int j = 0; j < m1.n; j++) {
            v2[i] += v1[j] * m1[j][i];
        }
    }
    return v2;
}

MAT &luFact(MAT &m1) {
    int n = m1.dim();
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            m1[j][i] /= m1[i][i];
        }
        for (int j = i + 1; j < n; j++) {
            for (int k = i + 1; k < n; k++) {
                m1[j][k] -= m1[j][i] * m1[i][k];
            }
        }
    }
    return m1;
}

// this function is based on an algorithm searched from the internet
MAT &cholesky(MAT &A) {
    int n = A.dim(); // and I transferred it to c++ code from pseudo code
    MAT *L;
    L = new MAT(n);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < (i + 1); j++) {
            double s = 0;
            for (int k = 0; k < j; k++) {
                s += (*L)[i][k] * (*L)[j][k];
            }
            (*L)[i][j] =
            (i == j) ? sqrt(A[i][i] - s) : (1.0 / ((*L)[j][j]) * (A[i][j] - s));
        }
    }
    
    return *L;
}

// forward substitutions function
VEC fwdSubs(MAT &m1, VEC b) {
    int len = b.len();
    
    VEC V(len);
    
    for (int i = 0; i < len; i++) {
        double s = 0;
        for (int j = 0; j < i; j++) {
            s += m1[i][j] * V[j];
        } // s=m1[i, 0]*V[0]+m1[i, 1]*V[1]+.....+m1[i, i-1]*V[i-1]
        V[i] = (b[i] - s) / m1[i][i]; // V[i]=(b[i]-s)/m1[i][i];
    }
    
    return V;
}

// backward substitution function
VEC bckSubs(MAT &m1, VEC b) {
    int len = b.len();
    
    VEC V(len);
    
    for (int i = 0; i < len; i++) {
        V[i] = b[i];
    }
    for (int i = len - 1; i >= 0; i--) {
        V[i] = V[i] / m1[i][i];
        for (int j = i - 1; j >= 0; j--) {
            V[j] = V[j] - (m1[j][i] * V[i]);
        }
    }
    return V;
}

VEC choSolve(MAT &L, VEC b) {
    // Ax=b, L*LT*x=b. Solve Ly=b first, then solve LT*x=y. x is the result.
    
    VEC y = fwdSubs(L, b); // Solve Ly=b with forward substitution
    MAT LT = L.tpose();
    VEC x = bckSubs(LT, y); // Solve LT*x=y
    
    return x;
}

int jacobi(MAT &A, VEC b, VEC &x, int maxIter, double tol) {
    int iter, n;
    n = x.len();
    VEC x_next(n);                            // vector to store updated x
    for (iter = 1; iter <= maxIter; iter++) { // loop with each iteration
        for (int i = 0; i < n; i++) {
            double temp1 = 0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    temp1 += A[i][j] * x[j];
                }
            }
            x_next[i] = (b[i] - temp1) / A[i][i];
        }
        double error = p_norm(x, x_next, -1); // p_norm is defined in VEC.h
        
        if (error < tol) { // if error is smaller then tol, stop and return
            x = x_next;
            return iter;
        }
        x = x_next;
    }
    x = x_next;
    return iter;
}

int gaussSeidel(MAT &A, VEC b, VEC &x, int maxIter, double tol) {
    int iter, n;
    n = x.len();
    VEC x_next(n);                            // vector to store updated x
    for (iter = 1; iter <= maxIter; iter++) { // loop with each iteration
        for (int i = 0; i < n; i++) {
            double temp1 = 0, temp2 = 0;
            for (int j = 0; j < i; j++) {
                // use x_next because it is already be updated
                temp1 += A[i][j] * x_next[j];
            }
            for (int j = i + 1; j < n; j++) {
                temp2 += A[i][j] * x[j];
            }
            x_next[i] = (b[i] - temp1 - temp2) / A[i][i];
        }
        
        double error = p_norm(x_next, x, 2);
        
        if (error < tol) { // if error is smaller then tol, stop and return
            x = x_next;
            return iter;
        }
        x = x_next;
    }
    x = x_next;
    return iter;
}

int sgs(MAT &A, VEC b, VEC &x, int maxIter, double tol) {
    int iter, n;
    n = x.len();
    
    VEC x_half(n); // vector to store updated x from Forward Gauss Seidel
    VEC x_next(n); // vector to store updated x from Backward Gauss Seidel
    for (iter = 1; iter <= maxIter; iter++) { // loop with each iteration
        for (int i = 0; i < n; i++) {           // Forward Gauss Seidel
            double temp1 = 0, temp2 = 0;
            for (int j = 0; j < i; j++) {
                temp1 += A[i][j] * x_half[j];
            }
            for (int j = i + 1; j < n; j++) {
                temp2 += A[i][j] * x[j];
            }
            x_half[i] = (b[i] - temp1 - temp2) / A[i][i];
        }
        
        for (int i = n - 1; i > 0; i--) { // Backward Gauss Seidel
            double temp1 = 0, temp2 = 0;
            for (int j = n - 1; j > i; j--) {
                temp1 += A[i][j] * x_next[j];
            }
            for (int j = i - 1; j > 0; j--) {
                temp2 += A[i][j] * x_half[j];
            }
            x_next[i] = (b[i] - temp1 - temp2) / A[i][i];
        }
        
        double error = p_norm(x_next, x, 2);
        
        if (error < tol) { // if error is smaller then tol, stop and return
            x = x_next;
            return iter;
        }
        x = x_next;
    }
    x = x_next;
    return iter;
}

// Conjugate Gradient Method function
int cg(MAT &A, VEC b, VEC &x, int maxIter, double tol) {
    int iter, n;
    n = x.len();
    
    VEC r = b - A * x;
    VEC p = r;
    double alpha = 0, beta = 0, error = 0;
    
    for (iter = 1; iter <= maxIter; iter++) { // loop with each iteration
        VEC temp1 = A * p;                      // only need to calculate p*A once
        double temp2 = temp1 * p;               // only need to calculate p*A*p once
        alpha = (p * r) / temp2;
        x += p * alpha;
        r -= temp1 * alpha;
        beta = (temp1 * r) / temp2;
        p = r - p * beta;
        error = sqrt((r * r) / r.len()); // error defined in hw06
        if (error < tol) {
            return iter;
        }
    }
    return iter;
}

/*
 Solve x in Ax=b with LU method
 */
VEC LU_Solve(MAT A, VEC b) {
    int n = A.dim();
    
    A = luFact(A);
    MAT L(n);
    MAT U(n);
    
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            L[i][j] = A[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        L[i][i] = 1;
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            U[i][j] = A[i][j];
        }
    }
    
    VEC y = fwdSubs(L, b);
    VEC x = bckSubs(U, y);
    
    return x;
}

int EVpwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter) {
    double e = tol + 1;
    q0 = q0 / Lp_norm(q0, 2);
    int k = 0;
    
    VEC temp1(q0.len()); // temporary vector for optimization
    VEC temp2 = A * q0;  // temporary vector for optimization
    
    VEC z(q0.len());
    VEC r(q0.len());
    VEC u(q0.len());
    VEC w(q0.len());
    
    double v = 0;
    
    // the following is based on the algorithm from the lecture
    while (e > tol && k < maxiter) {
        z = temp2; // z=A*q0
        k++;
        q0 = z / Lp_norm(z, 2); // Lp_norm() was defined in VEC.h
        temp1 = q0 * A;
        v = temp1 * q0; // q0*A*q0
        
        temp2 = A * q0;
        r = temp2 - v * q0; // r=A*q0-v*q0
        u = temp1;          // u=q0*A
        w = u / Lp_norm(u, 2);
        
        e = Lp_norm(r, 2) / fabs(w * q0);
        
        lambda = v;
    }
    
    return k; // return iteration number
}

int EViPwr(MAT &A, VEC &q0, double &lambda, double tol, int maxiter) {
    double e = tol + 1;
    q0 = q0 / Lp_norm(q0, 2);
    int k = 0;
    
    VEC temp1(q0.len()); // temporary vector for optimization
    
    VEC z(q0.len());
    VEC r(q0.len());
    VEC u(q0.len());
    VEC w(q0.len());
    
    double v = 0;
    
    // the following is based on the algorithm from the lecture
    while (e > tol && k < maxiter) {
        // LU_Solve() was defined in MAT.h, which solves Ax=q0
        z = LU_Solve(A, q0);
        k++;
        q0 = z / Lp_norm(z, 2);
        temp1 = q0 * A;
        v = temp1 * q0; // v=q0*A*q0
        
        r = A * q0 - v * q0;
        u = temp1; // u=q0*A
        w = u / Lp_norm(u, 2);
        
        e = Lp_norm(r, 2) / fabs(w * q0);
        
        lambda = v;
    }
    
    return k; // return iteration number
}

int EViPwrShft(MAT &A, VEC &q0, double &lambda, double mu, double tol,
               int maxiter) {
    double e = tol + 1;
    q0 = q0 / Lp_norm(q0, 2);
    int k = 0;
    
    VEC temp1(q0.len()); // temporary vector for optimization
    
    VEC z(q0.len());
    VEC r(q0.len());
    VEC u(q0.len());
    VEC w(q0.len());
    
    double v = 0;
    
    MAT A_temp = A;
    
    for (int i = 0; i < A_temp.dim(); i++) {
        A_temp[i][i] -= mu; // A_temp=A-wI
    }
    
    // the following is based on the algorithm from the lecture
    while (e > tol && k < maxiter) {
        z = LU_Solve(A_temp, q0);
        k++;
        q0 = z / Lp_norm(z, 2);
        temp1 = q0 * A;
        v = temp1 * q0; // v=q0*A*q0
        
        r = A * q0 - v * q0;
        u = temp1; // u=q0*A
        w = u / Lp_norm(u, 2);
        
        e = Lp_norm(r, 2) / fabs(w * q0);
        
        lambda = v;
    }
    
    return k; // return iteration number
}

/*
 input Q, R as reference, will do QR factorization, make A=QR
 */
void qrFact(MAT A, MAT &Q, MAT &R) {
    int n = A.dim();
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Q[i][j] = 0;
            R[i][j] = 0;
        }
    }
    
    MAT AT = A.tpose();
    
    R[0][0] = sqrt(AT[0] * AT[0]);
    Q[0] = AT[0] / R[0][0];
    
    for (int j = 1; j < n; j++) {
        Q[j] = AT[j];
        for (int i = 0; i < j; i++) {
            R[i][j] = Q[i] * Q[j];
            Q[j] -= R[i][j] * Q[i];
        }
        R[j][j] = sqrt(Q[j] * Q[j]);
        Q[j] /= R[j][j];
    }
    
    Q = Q.tpose();
}

int EVqr(MAT &A, double tol, int maxiter) {
    int n = A.dim();
    
    MAT A_next = A;
    MAT Q(n);
    MAT R(n);
    
    int k;
    double err = std::numeric_limits<double>::max();
    
    for (k = 0; k < maxiter && tol < err; k++) {
        qrFact(A_next, Q, R); // do QR factorization, QR=A_next
        
        A = Q * R;
        A_next = R * Q;
        
        err = 0;
        for (int i = 1; i < n; i++) {
            if (err < fabs(A[i][i - 1])) {
                err = fabs(A[i][i - 1]);
            }
        }
    }
    return k;
}

int EVqrShifted(MAT &A, double mu, double tol, int maxiter) {
    int n = A.dim();
    
    MAT A_next = A;
    MAT Q(n);
    MAT R(n);
    
    int k;
    double err = std::numeric_limits<double>::max();
    
    // A_next=A_next-uI
    for (int i = 0; i < n; i++) {
        A_next[i][i] -= mu;
    }
    
    for (k = 0; k < maxiter && tol < err; k++) {
        qrFact(A_next, Q, R); // do QR factorization, QR=A_next
        
        A = Q * R;
        A_next = R * Q;
        
        err = 0;
        for (int i = 1; i < n; i++) {
            if (err < fabs(A[i][i - 1])) {
                err = fabs(A[i][i - 1]);
            }
        }
    }
    
    // A=A+uI
    for (int i = 0; i < n; i++) {
        A[i][i] += mu;
    }
    
    return k;
}

double Lagrange(double x, VEC &XDATA, VEC &YDATA) {
    int n = XDATA.len();
    double s = 0;
    
    for (int i = 0; i < n; i++) {
        double L = 1;
        for (int k = 0; k < n; k++) {
            if (k != i) {
                L *= (double(x) - XDATA[k]) / (XDATA[i] - XDATA[k]);
            }
        }
        
        s += YDATA[i] * L; // s = sum(yi*Li(x))
    }
    
    return s;
}

double NEV(double x, VEC &XDATA, VEC &YDATA) {
    int n = XDATA.len();
    VEC NS(n);
    for (int i = 0; i < n; i++) {
        NS[i] = YDATA[i];
    }
    for (int k = 1; k < n; k++) {
        for (int j = 0; j < n - k; j++) {
            NS[j] = ((x - XDATA[j]) * NS[j + 1] - (x - XDATA[k + j]) * NS[j]) /
            (XDATA[j + k] - XDATA[j]);
        }
    }
    return NS[0];
}

void splineM(int N, VEC &X, VEC &Y, VEC &M) {
    
    VEC h(N);
    for (int i = 1; i < N; i++) {
        h[i] = X[i] - X[i - 1];
    }
    
    VEC mu(N);
    VEC lambda(N);
    VEC d(N);
    
    for (int i = 1; i < N - 1; i++) { // construct mu, lambda, and d vector
        mu[i] = h[i] / (h[i] + h[i + 1]);
        lambda[i] = h[i + 1] / (h[i] + h[i + 1]);
        d[i] = (6 / (h[i] + h[i + 1])) *
        (((Y[i + 1] - Y[i]) / h[i + 1]) - ((Y[i] - Y[i - 1]) / h[i]));
    }
    lambda[0] = 0;
    d[0] = 0; // use zero boundary
    mu[N - 1] = 0;
    d[N - 1] = 0; // use zero boundary
    
    MAT A(N); // construct A for A * M = d
    for (int i = 0; i < N; i++) {
        A[i][i] = 2;
        if (i + 1 < N) {
            A[i + 1][i] = mu[i + 1];
            A[i][i + 1] = lambda[i];
        }
    }
    
    M = LU_Solve(A, d); // solve M in A * M = d
}

double spline(double x, int N, VEC &X, VEC &Y, VEC &M) {
    
    VEC h(N); // length of subintervals
    
    for (int i = 1; i < N; i++) {
        h[i] = X[i] - X[i - 1];
    }
    
    VEC A(N);
    VEC B(N);
    
    for (int i = 1; i < N; i++) { // construct A and B
        B[i] = Y[i - 1] - M[i - 1] * h[i] * h[i] / 6;
        A[i] = (Y[i] - Y[i - 1]) / h[i] - h[i] * (M[i] - M[i - 1]) / 6;
    }
    
    int i = 1;
    while (i < N) { // find which subinterval x is in
        if (x >= X[i - 1] && x <= X[i]) {
            break;
        }
        i++;
    }
    
    // caculate s, which is the spline interpolation value at x
    double s = M[i - 1] * pow(X[i] - x, 3) / (6 * h[i]) +
    M[i] * pow(x - X[i - 1], 3) / (6 * h[i]) + A[i] * (x - X[i - 1]) +
    B[i];
    
    return s;
}

/*
 perform parametrix spline interpolation
 X and Y: input sopport points
 N: number of soupport points
 spline_x and spline_y: will store the result of spline interpolation values
                        into them
 subintervalNumber: determine the number of interpolation values will be
                    produced, should be more than the number of support points
 */
void parametric_spline(int subintervalNumber, int N, VEC &X, VEC &Y,
                       VEC &spline_x, VEC &spline_y) {
    
    VEC t(N);
    
    for (int i = 1; i < N; i++) {
        t[i] = t[i - 1] +
        sqrt(pow(X[i] - X[i - 1], 2) + pow(Y[i] - Y[i - 1], 2));
    }

    // perform cubic spline interpolation on t, x
    VEC Mx(N);
    splineM(N, t, X, Mx);
    for (int i = 0; i < subintervalNumber; i++) {
        spline_x[i] = spline(i * t[N - 1] / (subintervalNumber - 1),
                             N, t, X, Mx);
    }
    
    // perform cubic spline interpolation on t, y
    VEC My(N);
    splineM(N, t, Y, My);
    for (int i = 0; i < subintervalNumber; i++) {
        spline_y[i] = spline(i * t[N - 1] / (subintervalNumber - 1),
                             N, t, Y, My);
    }
}
