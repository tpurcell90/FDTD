#ifndef FDTD_OBJECT
#define FDTD_OBJECT

#include <string>
#include <complex>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <iostream>
#include <memory>

enum Shape {sphere,block,ellipse,cone,cylinder};
typedef std::complex<double> cplx;

// Needs modification for complex I think
class Obj
{
protected:
    double* geoParam; 
    double* material;
    Shape part;
  

public:
    // Constructor
    Obj(Shape s, double* mater, double* geo);

    // Copy Constructor
    Obj(const Obj& o);


    // Dielectric
    cplx dielectric(double freq, double material[]);

    // Geometry Creators Should be in the field section
    //Obj makeSphere(double* mater, double rad);
};

#endif

