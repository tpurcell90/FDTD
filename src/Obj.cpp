#include "Obj.hpp"
#include <memory>
#include <math.h>

using namespace std;

//Dielectric
cplx Obj::dielectric(double freq) //, double material_[])
{
    cplx eps(material_[0],0.0);
    cplx i(0,1);
    for(int ii = 0; ii < (material_.size()-1)/3; ii++)
       eps += (pow(material_[3*ii+1],2.0) * material_[3*ii+2]) / (pow(material_[3*ii+1],2.0) - pow(freq,2.0) - i*freq*material_[3*ii+3]/(2*M_PI));
    return eps;
}

//Determine if a point is in an object
bool Obj::isObj(std::vector<double> v)
{
    bool isIn = true;
    if(part_ == block)
    {
        for(int ii = 0; ii < v.size(); ii++)
        {
            if((v[ii] > location_[ii] + geoParam_[ii]/2.0) || (v[ii] < location_[ii] - geoParam_[ii]/2.0))
                isIn = false;
        }
    }
    else if(part_ == sphere && dist(v,location_) > geoParam_[0])
        isIn = false;
    else
        isIn = false;
    return isIn;
}
double Obj::dist(vector<double> pt1, vector<double> pt2)
{
    double sum = 0;
    for(int cc = 0; cc < pt1.size(); cc ++)
        sum += pow((pt1[cc]-pt2[cc]),2);
    return sqrt(sum);
}









