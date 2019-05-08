
// vector class
#ifndef VEC_H
#define VEC_H

#include "iostream"

class VEC {
public:
    VEC(int n);                    // uninit constructor, val set to 0
    VEC(const VEC &v1);            // copy constructor
    VEC(int n, double *v);         // init constructor
    ~VEC();                        // destructor
    //    int len();
    int len() const;               // dimension of the vector
    VEC &operator-();              // unary operator, nerative value
    VEC &operator=(const VEC v1);  // assignment
    VEC &operator+=(const VEC v1); // V+=v1;
    VEC &operator-=(const VEC v1); // V-=v1;
    VEC &operator*=(double a);     // V*=dbl;
    VEC &operator/=(double a);     // V/+dbl
    VEC operator+(const VEC v1);   // V+v1
    VEC operator-(const VEC v1);   // V-v1
    double operator*(VEC v1);      // inner product
    VEC operator*(double a);       // V*dbl
    VEC operator/(double a);       // V/dbl
    double &operator[](int n);     // indexing
    friend VEC operator*(double a, const VEC v1); // dbl x V
    friend VEC *newVEC(int n);                    // create dynamic VEC
    friend VEC *newVEC(int n, double *v);

    
    void set_element(double x, int index);
    double get_element(int index) const;
    
    void sort();

    
    void print();
    void print(std::string s);
    
private:
    int dim;     // vector length
    double *val; // array to store vector
};

VEC operator*(double a, const VEC v1);
VEC *newVEC(int n); // create dynamic VEC
VEC *newVEC(int n, double *v); // create dynamic VEC


double Lp_norm(VEC v, double p);         // return Lp-norm of v
double p_norm(VEC v1, VEC v2, double p); // return p-norm of v1 and v2

#endif
