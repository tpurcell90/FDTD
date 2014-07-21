#include "Obj.hpp"
#include <memory>
#include <math.h>

using namespace std;


//Constructor
//Learn how to throw an error is mater does not have a size of 3n+1, where n is an integer
Obj::Obj(Shape s, vector<double> mater, vector<double> geo, vector<double> loc)  //: part(s), material(mater), geoParam(geo) {}
{
    part = s;
    material = mater;
    geoParam = geo;
    location = loc;
}
Obj::Obj(const Obj& o)  //: part(o.part), material(o.material), geoParam(o.geoParam){}
{
    part = o.part;
    material = o.material;
    geoParam = o.geoParam;
    location = o.location;
}

//Access Functions
vector<double> Obj::geo() {return geoParam;}
vector<double> Obj::loc() {return location;}

//Dielectric
cplx Obj::dielectric(double freq) //, double material[])
{
    cplx eps(material[0],0.0);
    cplx i(0,1);
    for(int ii = 0; ii < (material.size()-1)/3; ii++)
       eps += (pow(material[3*ii+1],2.0) * material[3*ii+2]) / (pow(material[3*ii+1],2.0) - pow(freq,2.0) - i*freq*material[3*ii+3]/(2*M_PI));
    return eps;
}

//Determine if a point is in an object
bool Obj::isObj(std::vector<double> v)
{
    bool isIn = true;
    if(part == block)
    {
        for(int ii = 0; ii < v.size(); ii++)
        {
            if((v[ii] > location[ii] + geoParam[ii]/2.0) || (v[ii] < location[ii] - geoParam[ii]/2.0))
                isIn = false;
        }
    }
    else if(part == sphere)
    {
        if(dist(v,location) > geoParam[0])
            isIn = false;
    }
    else
    {
        isIn = false;
    }
    return isIn;
}
double Obj::dist(vector<double> pt1, vector<double> pt2)
{
    double sum = 0;
    for(int cc = 0; cc < pt1.size(); cc ++)
        sum += pow((pt1[cc]-pt2[cc]),2);
    return sqrt(sum);
}









