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
    std::vector<double> geoParam;
    std::vector<double> material;
    Shape part;
    std::vector<double> location;
 
public:
    // Constructor
    Obj(Shape s, std::vector<double> mater, std::vector<double> geo, std::vector<double> loc);

    // Copy Constructor
    Obj(const Obj& o);
    // Access Functions
    Shape s() {return part;}
    std::vector<double> geo();
    std::vector<double> loc();

    // Dielectric Access function
    cplx dielectric(double freq);
    bool isObj(std::vector<double> v);
    double dist(std::vector<double> pt1,std::vector<double> pt2);
};

#endif

