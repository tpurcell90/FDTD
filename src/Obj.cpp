#include "Obj.hpp"
#include <memory>
#include <math.h>

using namespace std;

typedef std::complex<double> cplx;

//Constructor
Obj::Obj(Shape s, double* mater, double* geo)  //: part(s), material(mater), geoParam(geo) {}
{
    part = s;
    material = mater;
    geoParam = geo;
}
Obj::Obj(const Obj& o)  //: part(o.part), material(o.material), geoParam(o.geoParam){}
{
    part = o.part;
    material = o.material;
    geoParam = o.geoParam;
}

//Dielectric 
//
//Not sure how to get the size of the array...
complex<double> dielectric(double freq, double material[])
{
    complex<double> eps(0.0,0.0);
    complex<double> i(0,1);
    for(int ii = 0; ii < 6/3; ii++)
       eps += (pow(material[3*ii],2.0) * material[3*ii+1]) / (pow(material[3*ii],2.0) - pow(freq,2.0) - i*freq*material[3*ii+2]/(2*M_PI));
    return eps;
}

// I know this should be in the FDTDField class I just don't want to change it yet
/*Obj::Obj makeSphere(double* mater, double rad)
{
    double geo[1] = [rad];
    return Obj(sphere, mater, geo);
}*/
