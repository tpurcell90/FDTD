#include "Obj.hpp"
#include <memory>
#include <math.h>

using namespace std;


//Constructor
//Learn how to throw an error is mater does not have a size of 3n+1, where n is an integer
Obj::Obj(Shape s, vector<double> mater, vector<double> geo)  //: part(s), material(mater), geoParam(geo) {}
{
    part = s;
    vector<double> material(mater);
    vector<double> geoParam(geo);
}
Obj::Obj(const Obj& o)  //: part(o.part), material(o.material), geoParam(o.geoParam){}
{
    part = o.part;
    vector<double> material(o.material);
    vector<double> geoParam (o.geoParam);
}

//Access Functions
vector<double> Obj::geo() {return geoParam;}

//Dielectric
cplx Obj::dielectric(double freq) //, double material[])
{
    cplx eps(material[0],0.0);
    cplx i(0,1);
    for(int ii = 0; ii < (material.size()-1)/3; ii++)
       eps += (pow(material[3*ii+1],2.0) * material[3*ii+2]) / (pow(material[3*ii+1],2.0) - pow(freq,2.0) - i*freq*material[3*ii+3]/(2*M_PI));
    return eps;
}

