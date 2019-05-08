
// VEC class functions
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>

#include "VEC.h"

using namespace std;

VEC::VEC(int n) {
    dim = n;
    val = (double *)calloc(n, sizeof(double));
}

VEC::VEC(const VEC &v1) {
    dim = v1.dim;
    val = (double *)calloc(dim, sizeof(double));
    for (int i = 0; i < dim; i++) {
        val[i] = v1.val[i];
    }
}

VEC::VEC(int n, double *v) {
    dim = n;
    val = (double *)calloc(n, sizeof(double));
    for (int i = 0; i < n; i++) {
        val[i] = v[i];
    }
}

VEC::~VEC() { free(val); }

int VEC::len() const { return dim; }

VEC &VEC::operator-() {
    for (int i = 0; i < dim; i++) {
        val[i] = -val[i];
    }
    return *this;
}

VEC &VEC::operator=(const VEC v1) {
    dim = v1.dim;
    for (int i = 0; i < dim; i++) {
        val[i] = v1.val[i];
    }
    return *this;
}

VEC &VEC::operator+=(const VEC v1) {
    for (int i = 0; i < dim; i++) {
        val[i] += v1.val[i];
    }
    return *this;
}

VEC &VEC::operator-=(const VEC v1) {
    for (int i = 0; i < dim; i++) {
        val[i] -= v1.val[i];
    }
    return *this;
}

VEC &VEC::operator*=(double a) {
    for (int i = 0; i < dim; i++) {
        val[i] *= a;
    }
    return *this;
}

VEC &VEC::operator/=(double a) {
    for (int i = 0; i < dim; i++) {
        val[i] /= a;
    }
    return *this;
}

VEC VEC::operator+(const VEC v1) {
    VEC s(*this);
    for (int i = 0; i < dim; i++) {
        s.val[i] += v1.val[i];
    }
    return s;
}

VEC VEC::operator-(const VEC v1) {
    VEC s(*this);
    for (int i = 0; i < dim; i++) {
        s.val[i] -= v1.val[i];
    }
    return s;
}

double VEC::operator*(VEC v1) {
    double s = 0;
    for (int i = 0; i < dim; i++) {
        s += val[i] * v1.val[i];
    }
    return s;
}

VEC VEC::operator*(double a) {
    VEC s(*this);
    for (int i = 0; i < dim; i++) {
        s.val[i] *= a;
    }
    return s;
}

VEC VEC::operator/(double a) {
    VEC s(*this);
    for (int i = 0; i < dim; i++) {
        s.val[i] /= a;
    }
    return s;
}

double &VEC::operator[](int n) {
    if (n < 0)
        n = 0;
    else if (n >= dim)
        n = dim - 1;
    return val[n];
}

VEC *newVEC(int n) {
    VEC *vptr;
    vptr = (VEC *)malloc(sizeof(VEC));
    vptr->dim = n;
    vptr->val = (double *)calloc(n, sizeof(double));
    return vptr;
}

VEC *newVEC(int n, double *v) {
    VEC *vptr;
    vptr = (VEC *)malloc(sizeof(VEC));
    vptr->dim = n;
    vptr->val = (double *)calloc(n, sizeof(double));
    for (int i = 0; i < n; i++) {
        vptr->val[i] = v[i];
    }
    return vptr;
}

VEC operator*(double a, const VEC v1) {
    int n = v1.len();
    VEC v2 = v1;
    for (int i = 0; i < n; i++) {
        v2[i] = a * v2[i];
    }
    return v2;
}

void merge(VEC &s, int l, int m, int r) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
    
    VEC L(n1);
    VEC R(n2);
    
    for (i = 0; i < n1; i++)
        L[i] = s[l + i];
    for (j = 0; j < n2; j++)
        R[j] = s[m + 1 + j];
    
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            s[k] = L[i];
            i++;
        } else {
            s[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        s[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        s[k] = R[j];
        j++;
        k++;
    }
}

void mergeSort(VEC &s, int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;
        
        mergeSort(s, l, m);
        mergeSort(s, m + 1, r);
        
        merge(s, l, m, r);
    }
}

void VEC::sort() {
    mergeSort((*this), 0, (*this).len()-1);
}

void VEC::print() {
    for (int i = 0; i < dim; i++) {
        printf("%g ", val[i]);
    }
    cout << endl;
}

void VEC::print(string s) {
    for (int i = 0; i < dim; i++) {
        printf("%g%s", val[i], s.c_str());
    }
    cout << endl;
}

/*
 return Lp-norm of v
 when n=-1, return infinite norm
 */
double Lp_norm(VEC v, double p) {
    double s = 0;
    if (p == -1) {
        for (int i = 0; i < v.len(); i++) {
            if (s < v[i]) {
                s = v[i];
            }
        }
    } else {
        for (int i = 0; i < v.len(); i++) {
            s += pow(fabs(v[i]), p);
        }
    }
    s = pow(s, 1 / p);
    return s;
}

double p_norm(VEC v1, VEC v2, double p) {
    int n = v1.len();
    double s = 0;
    double temp;
    if (p == -1) {
        for (int i = 0; i < n; i++) {
            temp = fabs(v1[i] - v2[i]);
            if (s < temp) {
                s = temp;
            }
        }
    } else {
        for (int i = 0; i < n; i++) {
            temp = fabs(v1[i] - v2[i]);
            s += pow(temp, p);
        }
    }
    s = pow(s, 1 / p);
    return s;
}
