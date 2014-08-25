#include "Obj.hpp"
#include <memory>
#include <math.h>

using namespace std;

//Dielectric
/**
 * @brief Get the dielectric constant of the Object's Material
 * @details Calculates the complex dielectric function of the object's mateial at a given frequency
 *
 * @param freq Frequency of light
 * @return dielectric function at the frequncy
 */
double Obj::dielectric(double freq) //, double material_[])
{
    //cplx eps(material_[0],0.0);
    //cplx i(0,1);
    //for(int ii = 0; ii < (material_.size()-1)/3; ii++)
    //   eps += (pow(material_[3*ii+1],2.0) * material_[3*ii+2]) / (pow(material_[3*ii+1],2.0) - pow(freq,2.0) - i*freq*material_[3*ii+3]/(2*M_PI));
    return material_[0];
}

/**
 * @brief Determine if a point is in an object
 * @details Deterimes a point in vector form is inside the object or not
 *
 * @param v point to test if it's in an object
 * @return True if inside the object false if not
 */
bool Obj::isObj(std::vector<double> v)
{
    bool isIn = true;
    if(part_ == block)
    {
        for(int ii = 0; ii < v.size(); ii++)
        {
            if((v[ii] > location_[ii] + geoParam_[ii]/2.0 + 1e-5) || (v[ii] < location_[ii] - geoParam_[ii]/2.0 -1e-5))// add a 1/10 res term here
                isIn = false;
        }
    }
    else if(part_ == sphere && dist(v,location_) > geoParam_[0])
        isIn = false;
    else
        isIn = false;
    return isIn;
}
/**
 * @brief Determines teh distance between two points
 *
 * @param pt1 The Second point
 * @param pt2 The Second point
 * @return The distance between points
 */
double Obj::dist(vector<double> pt1, vector<double> pt2)
{
    double sum = 0;
    for(int cc = 0; cc < pt1.size(); cc ++)
        sum += pow((pt1[cc]-pt2[cc]),2);
    return sqrt(sum);
}









