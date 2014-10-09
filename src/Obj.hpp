#ifndef FDTD_OBJECT
#define FDTD_OBJECT

#include <string>
#include <complex>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <iostream>
#include <memory>

typedef std::complex<double> cplx;
enum Shape {sphere,block,ellipse,cone,cylinder};

// Needs modification for complex I think
class Obj
{

protected:
    Shape part_;
    std::vector<double> material_;
    std::vector<double> geoParam_;
    std::vector<double> location_;
    std::vector<double> alpha_, zi_, gamma_;
    std::vector<double> upConsts_;

public:

    // Constructor
    Obj(Shape s, std::vector<double> mater, std::vector<double> geo, std::vector<double> loc);
    // Copy Constructor
    Obj(const Obj& o) :  part_(o.part_), material_(o.material_), geoParam_(o.geoParam_), location_(o.location_), alpha_(o.alpha_), zi_(o.zi_), gamma_(o.gamma_){}

    void setUpConsts (double dt);

    // Access Functions
    Shape s() {return part_;}
    std::vector<double> geo() {return geoParam_;}
    std::vector<double> loc() {return location_;}
    std::vector<double> mat() {return material_;}
    std::vector<double> alpha() {return alpha_;}
    std::vector<double> zi() {return zi_;}
    std::vector<double> gamma() {return gamma_;}
    std::vector<double> upConsts() {return upConsts_;}
    // Dielectric Access function
    double dielectric(){return material_[0];}
    bool isObj(std::vector<double> v);
    double dist(std::vector<double> pt1,std::vector<double> pt2);
};

#endif

