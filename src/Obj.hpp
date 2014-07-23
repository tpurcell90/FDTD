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

public:
    // Constructor
    Obj(Shape s, std::vector<double> mater, std::vector<double> geo, std::vector<double> loc) : part_(s), material_(mater), geoParam_(geo), location_(loc) {}
    // Copy Constructor
    Obj(const Obj& o) :  part_(o.part_), material_(o.material_), geoParam_(o.geoParam_), location_(o.location_) {}
    // Access Functions
    Shape s() {return part_;}
    std::vector<double> geo() {return geoParam_;}
    std::vector<double> loc() {return location_;}

    // Dielectric Access function
    cplx dielectric(double freq);
    bool isObj(std::vector<double> v);
    double dist(std::vector<double> pt1,std::vector<double> pt2);
};

#endif

